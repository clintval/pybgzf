//! The `pybgzf._pybgzf` extension module.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, PoisonError, TryLockError};
use std::thread;

use pyo3::buffer::PyBuffer;
use pyo3::exceptions::{PyOSError, PyRuntimeError, PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

use crate::Error;
use crate::columns::{Columns, Kind};
use crate::reader::{
    AnyIndex, BgzfReader, IndexedReader as CoreIndexedReader, Query, ends_with_eof_marker,
};
use crate::sniff;
use crate::writer::{IndexFormat, IndexOptions, Writer as CoreWriter, check_options, csi_format};

const FILE_BUFFER: usize = 256 * 1024;

pyo3::create_exception!(
    pybgzf,
    TruncatedWarning,
    PyUserWarning,
    "A BGZF file ends without its end-of-file marker, so it may be truncated at a block boundary."
);

/// Warns that the BGZF data from `path`, or from a stream, may be truncated.
fn warn_truncated(py: Python<'_>, path: Option<&Path>) -> PyResult<()> {
    let data = path.map_or_else(|| "BGZF data".to_owned(), |path| path.display().to_string());
    let message = format!("{data} ends without an end-of-file marker and may be truncated");
    let message =
        std::ffi::CString::new(message).map_err(|e| PyValueError::new_err(e.to_string()))?;
    PyErr::warn(py, &py.get_type::<TruncatedWarning>(), &message, 1)
}

type ColumnsTuple = (i64, i64, Option<i64>, bool, String, i64, Kind);

/// The value of each `pybgzf.IndexFormat` member.
#[pyclass(
    module = "pybgzf._pybgzf",
    rename_all = "UPPERCASE",
    eq,
    frozen,
    hash,
    from_py_object
)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum IndexKind {
    Tbi,
    Csi,
}

enum Sink {
    File(BufWriter<File>),
    Python(Py<PyAny>),
    Closed,
}

fn python_to_io(error: PyErr) -> io::Error {
    io::Error::other(error)
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Sink::File(file) => file.write(buf),
            Sink::Closed => Err(io::Error::other("I/O operation on a closed writer")),
            Sink::Python(object) => Python::attach(|py| {
                let written = object
                    .bind(py)
                    .call_method1("write", (PyBytes::new(py, buf),))
                    .map_err(python_to_io)?;
                if written.is_none() {
                    Ok(buf.len())
                } else {
                    written.extract::<usize>().map_err(python_to_io)
                }
            }),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Sink::File(file) => file.flush(),
            Sink::Closed => Ok(()),
            Sink::Python(object) => Python::attach(|py| {
                let object = object.bind(py);
                if object.hasattr("flush").map_err(python_to_io)? {
                    object.call_method0("flush").map_err(python_to_io)?;
                }
                Ok(())
            }),
        }
    }
}

fn to_python(error: Error) -> PyErr {
    match error {
        Error::Invalid(message) => PyValueError::new_err(message),
        Error::Io(error) => error.downcast::<PyErr>().unwrap_or_else(PyErr::from),
    }
}

fn columns_from_tuple(columns: ColumnsTuple) -> PyResult<Columns> {
    let (refname, start, end, zero_based, meta_char, skip_lines, kind) = columns;
    let column = |value: i64| {
        usize::try_from(value)
            .map_err(|_| PyValueError::new_err("column numbers are 1-based and must be positive"))
    };
    let (refname, start, end) = (
        column(refname)?,
        column(start)?,
        end.map(column).transpose()?,
    );
    let skip_lines = u64::try_from(skip_lines)
        .map_err(|_| PyValueError::new_err("skip_lines must not be negative"))?;
    let meta_char = match meta_char.as_bytes() {
        [byte] if byte.is_ascii() => *byte,
        _ => {
            return Err(PyValueError::new_err(
                "meta_char must be a single ASCII character",
            ));
        }
    };
    let columns = Columns {
        refname,
        start,
        end,
        zero_based,
        meta_char,
        skip_lines,
        kind,
    };
    columns.validate().map_err(PyValueError::new_err)?;
    Ok(columns)
}

fn columns_to_tuple(columns: &Columns) -> ColumnsTuple {
    (
        columns.refname as i64,
        columns.start as i64,
        columns.end.map(|end| end as i64),
        columns.zero_based,
        char::from(columns.meta_char).to_string(),
        columns.skip_lines as i64,
        columns.kind,
    )
}

/// Returns an address that no other running thread shares.
fn this_thread() -> usize {
    thread_local! {
        static MARK: u8 = const { 0 };
    }
    MARK.with(|mark| std::ptr::from_ref(mark).addr())
}

/// A mutex that threads wait for without the GIL, and that raises on a reentrant call.
struct Lock<T> {
    mutex: Mutex<T>,
    owner: AtomicUsize,
}

struct Guard<'a, T> {
    inner: MutexGuard<'a, T>,
    owner: &'a AtomicUsize,
}

impl<T> Deref for Guard<'_, T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.inner
    }
}

impl<T> DerefMut for Guard<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.inner
    }
}

impl<T> Drop for Guard<'_, T> {
    fn drop(&mut self) {
        self.owner.store(0, Ordering::Relaxed);
    }
}

impl<T: Send> Lock<T> {
    fn new(value: T) -> Self {
        Self {
            mutex: Mutex::new(value),
            owner: AtomicUsize::new(0),
        }
    }

    fn guard<'a>(&'a self, inner: MutexGuard<'a, T>, thread: usize) -> Guard<'a, T> {
        self.owner.store(thread, Ordering::Relaxed);
        Guard {
            inner,
            owner: &self.owner,
        }
    }

    fn try_lock(&self) -> Option<Guard<'_, T>> {
        let inner = match self.mutex.try_lock() {
            Ok(inner) => inner,
            Err(TryLockError::Poisoned(poisoned)) => poisoned.into_inner(),
            Err(TryLockError::WouldBlock) => return None,
        };
        Some(self.guard(inner, this_thread()))
    }

    /// Waits for the lock; called without the GIL.
    fn lock(&self) -> PyResult<Guard<'_, T>> {
        let thread = this_thread();
        if self.owner.load(Ordering::Relaxed) == thread {
            return Err(PyRuntimeError::new_err("reentrant call"));
        }
        let inner = self.mutex.lock().unwrap_or_else(PoisonError::into_inner);
        Ok(self.guard(inner, thread))
    }

    /// Calls a quick `f`, releasing the GIL only to wait for another thread.
    fn with<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut T) -> PyResult<R> + Send,
    ) -> PyResult<R> {
        match self.try_lock() {
            Some(mut guard) => f(&mut guard),
            None => py.detach(|| f(&mut *self.lock()?)),
        }
    }

    fn get_mut(&mut self) -> &mut T {
        self.mutex.get_mut().unwrap_or_else(PoisonError::into_inner)
    }
}

fn closed_error() -> PyErr {
    PyValueError::new_err("I/O operation on closed file.")
}

fn open_writer(inner: &mut CoreWriter<Sink>) -> PyResult<&mut CoreWriter<Sink>> {
    if inner.is_finished() {
        return Err(closed_error());
    }
    Ok(inner)
}

/// The Rust half of `pybgzf.BgzfWriter`.
#[pyclass(module = "pybgzf._pybgzf", frozen)]
struct Writer {
    inner: Lock<CoreWriter<Sink>>,
    closed: AtomicBool,
}

impl Writer {
    fn end(
        &self,
        py: Python<'_>,
        end: fn(&mut CoreWriter<Sink>) -> crate::Result<()>,
    ) -> PyResult<bool> {
        py.detach(|| {
            let mut inner = self.inner.lock()?;
            let ended = end(&mut inner);
            *inner.get_mut() = Sink::Closed;
            self.closed.store(true, Ordering::Release);
            ended.map_err(to_python)?;
            Ok(inner.is_complete())
        })
    }
}

#[pymethods]
impl Writer {
    #[new]
    #[pyo3(signature = (dest, *, level, threads, index, index_path, columns, infer, infer_bed, csi_min_shift, csi_depth))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        dest: &Bound<'_, PyAny>,
        level: i64,
        threads: i64,
        index: Option<IndexKind>,
        index_path: Option<PathBuf>,
        columns: Option<ColumnsTuple>,
        infer: bool,
        infer_bed: bool,
        csi_min_shift: i64,
        csi_depth: Option<i64>,
    ) -> PyResult<Self> {
        let columns = columns.map(columns_from_tuple).transpose()?;
        let infer = infer || infer_bed;
        let options = match (index, index_path) {
            (Some(kind), Some(path)) => {
                let format = match kind {
                    IndexKind::Tbi => IndexFormat::Tabix,
                    IndexKind::Csi => csi_format(csi_min_shift, csi_depth).map_err(to_python)?,
                };
                Some(IndexOptions {
                    format,
                    path,
                    columns: if infer { None } else { columns },
                    bed_only: infer_bed,
                })
            }
            _ => None,
        };
        check_options(level, threads, options.as_ref()).map_err(to_python)?;
        let sink = if dest.is_instance_of::<PyString>() {
            let path: PathBuf = dest.extract()?;
            let file = File::create(path).map_err(|e| to_python(Error::Io(e)))?;
            Sink::File(BufWriter::with_capacity(FILE_BUFFER, file))
        } else {
            Sink::Python(dest.clone().unbind())
        };
        let inner = CoreWriter::new(sink, level, threads, options).map_err(to_python)?;
        Ok(Self {
            inner: Lock::new(inner),
            closed: AtomicBool::new(false),
        })
    }

    /// Writes bytes-like `data` and returns the number of bytes written.
    fn write(&self, py: Python<'_>, data: &Bound<'_, PyAny>) -> PyResult<usize> {
        let owned;
        let bytes: &[u8] = if let Ok(bytes) = data.cast::<PyBytes>() {
            bytes.as_bytes()
        } else {
            owned = PyBuffer::<u8>::get(data)?.to_vec(py)?;
            &owned
        };
        let write =
            |inner: &mut CoreWriter<Sink>| open_writer(inner)?.write(bytes).map_err(to_python);
        match self
            .inner
            .try_lock()
            .filter(|inner| !inner.may_block(bytes.len()))
        {
            Some(mut inner) => write(&mut inner)?,
            None => py.detach(|| write(&mut *self.inner.lock()?))?,
        }
        Ok(bytes.len())
    }

    /// Ends the current block and flushes the sink.
    fn flush(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| {
            open_writer(&mut *self.inner.lock()?)?
                .flush()
                .map_err(to_python)
        })
    }

    /// Returns the virtual position of the next byte.
    fn tell(&self, py: Python<'_>) -> PyResult<u64> {
        py.detach(|| {
            open_writer(&mut *self.inner.lock()?)?
                .tell()
                .map_err(to_python)
        })
    }

    /// Writes the end-of-file marker and the index, if any, once any call in progress returns,
    /// and returns whether both were written.
    fn close(&self, py: Python<'_>) -> PyResult<bool> {
        self.end(py, CoreWriter::finish)
    }

    /// Ends the file without the end-of-file marker and removes the index, once any call in
    /// progress returns.
    fn abandon(&self, py: Python<'_>) -> PyResult<()> {
        self.end(py, CoreWriter::abandon).map(drop)
    }

    #[getter]
    fn closed(&self) -> bool {
        self.closed.load(Ordering::Acquire)
    }

    #[getter]
    fn columns(&self, py: Python<'_>) -> PyResult<Option<ColumnsTuple>> {
        self.inner
            .with(py, |inner| Ok(inner.columns().map(columns_to_tuple)))
    }
}

/// Decides a column layout from lines, as `Columns.sniff` does.
#[pyclass(module = "pybgzf._pybgzf")]
#[derive(Default)]
struct Sniffer {
    inner: sniff::Sniffer,
}

#[pymethods]
impl Sniffer {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Looks at one line, without its terminator, and returns the layout once it is decided.
    fn push(&mut self, line: &[u8]) -> PyResult<Option<ColumnsTuple>> {
        let decided = self.inner.push(line).map_err(PyValueError::new_err)?;
        Ok(decided.as_ref().map(columns_to_tuple))
    }
}

/// Checks a column layout, raising `ValueError` if tabix could not record it.
#[pyfunction]
fn validate_columns(columns: ColumnsTuple) -> PyResult<()> {
    columns_from_tuple(columns).map(drop)
}

enum Source {
    File(BufReader<File>),
    Python { object: Py<PyAny>, seekable: bool },
}

impl Read for Source {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Source::File(file) => file.read(buf),
            Source::Python { object, .. } => Python::attach(|py| {
                let data = object
                    .bind(py)
                    .call_method1("read", (buf.len(),))
                    .map_err(python_to_io)?;
                if data.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::WouldBlock,
                        "the source has no data yet",
                    ));
                }
                let data = PyBuffer::<u8>::get(&data)
                    .map_err(python_to_io)?
                    .to_vec(py)
                    .map_err(python_to_io)?;
                if data.len() > buf.len() {
                    return Err(io::Error::other(format!(
                        "read({}) returned {} bytes",
                        buf.len(),
                        data.len()
                    )));
                }
                buf[..data.len()].copy_from_slice(&data);
                Ok(data.len())
            }),
        }
    }
}

impl Seek for Source {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        match self {
            Source::File(file) => file.seek(position),
            Source::Python {
                seekable: false, ..
            } => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "the source is not seekable",
            )),
            Source::Python { object, .. } => Python::attach(|py| {
                let (offset, whence) = match position {
                    SeekFrom::Start(offset) => (offset as i64, 0),
                    SeekFrom::Current(offset) => (offset, 1),
                    SeekFrom::End(offset) => (offset, 2),
                };
                object
                    .bind(py)
                    .call_method1("seek", (offset, whence))
                    .and_then(|at| at.extract::<u64>())
                    .map_err(python_to_io)
            }),
        }
    }
}

fn threads(threads: i64) -> PyResult<NonZero<usize>> {
    crate::check_threads(threads).map_err(PyValueError::new_err)
}

fn io_to_python(error: io::Error) -> PyErr {
    to_python(Error::Io(error))
}

/// Converts an error reading `path`, naming the file as Python's own errors do.
fn path_error(error: io::Error, path: Option<&Path>) -> PyErr {
    let Some(path) = path else {
        return io_to_python(error);
    };
    if error
        .get_ref()
        .is_some_and(<dyn std::error::Error + Send + Sync>::is::<PyErr>)
    {
        return io_to_python(error);
    }
    let name = path.display().to_string();
    match error.raw_os_error() {
        Some(code) => {
            let message = error.to_string();
            let message = message.split(" (os error").next().unwrap_or(&message);
            PyOSError::new_err((code, message.to_string(), name))
        }
        None => io_to_python(io::Error::new(error.kind(), format!("{name}: {error}"))),
    }
}

/// Stops a reader's threads, in the background when `background` is set.
fn stop<T: Send + 'static>(inner: T, background: bool) {
    if background {
        drop(thread::Builder::new().spawn(move || drop(inner)));
    } else {
        drop(inner);
    }
}

/// A reader until it is closed, which stops its threads without the GIL, in the background when
/// `background` is set, such as for a pipe whose writer is idle, so that closing never waits on
/// a read that may not return. Threads take turns with it, and `closed` never waits for a turn.
struct Open<T: Send + 'static> {
    inner: Lock<Option<T>>,
    background: bool,
    closed: AtomicBool,
}

impl<T: Send + 'static> Open<T> {
    fn new(inner: T, background: bool) -> Self {
        Self {
            inner: Lock::new(Some(inner)),
            background,
            closed: AtomicBool::new(false),
        }
    }

    /// Calls a quick `f` with the reader, releasing the GIL only to wait for another thread.
    fn with<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut T) -> PyResult<R> + Send,
    ) -> PyResult<R> {
        self.inner
            .with(py, |inner| f(inner.as_mut().ok_or_else(closed_error)?))
    }

    /// Calls `f` with the reader; called without the GIL.
    fn locked<R>(&self, f: impl FnOnce(&mut T) -> PyResult<R>) -> PyResult<R> {
        f(self.inner.lock()?.as_mut().ok_or_else(closed_error)?)
    }

    /// Calls `f` with the reader without the GIL.
    fn detached<R: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut T) -> PyResult<R> + Send,
    ) -> PyResult<R> {
        py.detach(|| self.locked(f))
    }

    /// Stops the reader once any call in progress returns.
    fn close(&self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| {
            let inner = {
                let mut inner = self.inner.lock()?;
                self.closed.store(true, Ordering::Release);
                inner.take()
            };
            if let Some(inner) = inner {
                stop(inner, self.background);
            }
            Ok(())
        })
    }

    fn closed(&self) -> bool {
        self.closed.load(Ordering::Acquire)
    }
}

impl<T: Send + 'static> Drop for Open<T> {
    fn drop(&mut self) {
        if let Some(inner) = self.inner.get_mut().take() {
            let background = self.background;
            Python::attach(|py| py.detach(|| stop(inner, background)));
        }
    }
}

/// The Rust half of `pybgzf.BgzfReader`.
#[pyclass(module = "pybgzf._pybgzf", frozen)]
struct Reader {
    inner: Open<BgzfReader<Source>>,
    path: Option<PathBuf>,
}

impl Reader {
    /// Reads with `f` without the GIL, then warns if the data ended without an end-of-file marker.
    fn read_with<T: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut BgzfReader<Source>) -> io::Result<T> + Send,
    ) -> PyResult<T> {
        let path = self.path.as_deref();
        let (value, truncated) = self.inner.detached(py, |inner| {
            let value = f(inner).map_err(|e| path_error(e, path))?;
            Ok((value, inner.take_missing_eof_marker()))
        })?;
        if truncated {
            warn_truncated(py, path)?;
        }
        Ok(value)
    }
}

#[pymethods]
impl Reader {
    #[new]
    #[pyo3(signature = (src, *, threads))]
    fn new(py: Python<'_>, src: &Bound<'_, PyAny>, threads: i64) -> PyResult<Self> {
        let threads = self::threads(threads)?;
        let (source, path, may_block) = if src.is_instance_of::<PyString>() {
            let path: PathBuf = src.extract()?;
            let file = File::open(&path).map_err(|e| path_error(e, Some(&path)))?;
            let regular = file
                .metadata()
                .map_err(|e| path_error(e, Some(&path)))?
                .is_file();
            let file = BufReader::with_capacity(FILE_BUFFER, file);
            (Source::File(file), Some(path), !regular)
        } else {
            let seekable = src.hasattr("seekable")? && src.call_method0("seekable")?.is_truthy()?;
            let object = src.clone().unbind();
            (Source::Python { object, seekable }, None, true)
        };
        let inner = py
            .detach(|| BgzfReader::new(source, threads))
            .map_err(|e| path_error(e, path.as_deref()))?;
        let inner = Open::new(inner, may_block && threads.get() > 1);
        Ok(Self { inner, path })
    }

    /// Reads up to `size` bytes, fewer only at the end of the stream.
    fn read<'py>(&self, py: Python<'py>, size: usize) -> PyResult<Bound<'py, PyBytes>> {
        let buf = self.read_with(py, |inner| {
            let mut buf = Vec::with_capacity(size);
            inner.take(size as u64).read_to_end(&mut buf)?;
            Ok(buf)
        })?;
        Ok(PyBytes::new(py, &buf))
    }

    /// Reads everything that is left.
    fn readall<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let buf = self.read_with(py, |inner| {
            let mut buf = Vec::new();
            inner.read_to_end(&mut buf)?;
            Ok(buf)
        })?;
        Ok(PyBytes::new(py, &buf))
    }

    /// Reads through the next newline, or at most `size` bytes when `size` is not negative.
    fn readline<'py>(&self, py: Python<'py>, size: i64) -> PyResult<Bound<'py, PyBytes>> {
        let limit = u64::try_from(size).unwrap_or(u64::MAX);
        let line = self.read_with(py, |inner| {
            let mut line = Vec::new();
            inner.take(limit).read_until(b'\n', &mut line)?;
            Ok(line)
        })?;
        Ok(PyBytes::new(py, &line))
    }

    /// Returns the virtual position of the next byte.
    fn tell(&self, py: Python<'_>) -> PyResult<u64> {
        self.inner.with(py, |inner| Ok(inner.virtual_position()))
    }

    /// Moves to a virtual position.
    fn seek(&self, py: Python<'_>, position: u64) -> PyResult<u64> {
        let path = self.path.as_deref();
        self.inner.detached(py, |inner| {
            inner.seek(position).map_err(|error| match error.kind() {
                io::ErrorKind::InvalidInput => PyValueError::new_err(error.to_string()),
                _ => path_error(error, path),
            })?;
            Ok(inner.virtual_position())
        })
    }

    /// Stops any worker threads.
    fn close(&self, py: Python<'_>) -> PyResult<()> {
        self.inner.close(py)
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.closed()
    }
}

/// The Rust half of `pybgzf.IndexedReader`.
#[pyclass(module = "pybgzf._pybgzf", frozen)]
struct IndexedReader {
    inner: Open<CoreIndexedReader<BufReader<File>>>,
    path: PathBuf,
}

#[pymethods]
impl IndexedReader {
    #[new]
    #[pyo3(signature = (path, index_path, *, threads))]
    #[allow(clippy::needless_pass_by_value)]
    fn new(py: Python<'_>, path: PathBuf, index_path: PathBuf, threads: i64) -> PyResult<Self> {
        let threads = self::threads(threads)?;
        let error = |error: io::Error, read: &Path| match error.kind() {
            io::ErrorKind::InvalidData => {
                PyValueError::new_err(format!("{}: {error}", index_path.display()))
            }
            _ => path_error(error, Some(read)),
        };
        let index = py
            .detach(|| AnyIndex::read(&index_path))
            .map_err(|e| error(e, &index_path))?;
        let inner = py
            .detach(|| CoreIndexedReader::new(BgzfReader::from_path(&path, threads)?, index))
            .map_err(|e| error(e, &path))?;
        if !py
            .detach(|| ends_with_eof_marker(&path))
            .map_err(|error| path_error(error, Some(&path)))?
        {
            warn_truncated(py, Some(&path))?;
        }
        let inner = Open::new(inner, false);
        Ok(Self { inner, path })
    }

    /// Starts a query for lines overlapping the 0-based, half-open `[start, end)` on `refname`.
    fn query(slf: Bound<'_, Self>, refname: &str, start: i64, end: i64) -> PyResult<QueryIterator> {
        let query = slf.get().inner.with(slf.py(), |inner| {
            inner
                .query(refname.as_bytes(), start, end)
                .map_err(io_to_python)
        })?;
        let state = Lock::new(QueryState {
            query,
            lines: VecDeque::new(),
        });
        Ok(QueryIterator {
            reader: slf.unbind(),
            state,
        })
    }

    #[getter]
    fn refnames(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.inner.with(py, |inner| {
            inner
                .names()
                .map(|name| {
                    String::from_utf8(name.to_vec()).map_err(|_| {
                        PyValueError::new_err(format!(
                            "the reference name {:?} in the index is not UTF-8",
                            String::from_utf8_lossy(name)
                        ))
                    })
                })
                .collect()
        })
    }

    #[getter]
    fn columns(&self, py: Python<'_>) -> PyResult<ColumnsTuple> {
        self.inner
            .with(py, |inner| Ok(columns_to_tuple(inner.columns())))
    }

    /// Stops any worker threads.
    fn close(&self, py: Python<'_>) -> PyResult<()> {
        self.inner.close(py)
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.closed()
    }
}

const QUERY_BATCH: usize = 256 * 1024;

struct QueryState {
    query: Query,
    lines: VecDeque<String>,
}

/// Lines overlapping a region, read in batches with the GIL released.
#[pyclass(module = "pybgzf._pybgzf", frozen)]
struct QueryIterator {
    reader: Py<IndexedReader>,
    state: Lock<QueryState>,
}

impl QueryIterator {
    /// Returns the next line, reading a batch if none is left; called without the GIL.
    fn next_line(&self) -> PyResult<Option<String>> {
        let mut state = self.state.lock()?;
        let state = &mut *state;
        let reader = self.reader.get();
        loop {
            if let Some(line) = state.lines.pop_front() {
                return Ok(Some(line));
            }
            if state.query.done {
                return Ok(None);
            }
            let mut batch = Vec::new();
            reader.inner.locked(|inner| {
                inner
                    .next_lines(&mut state.query, &mut batch, QUERY_BATCH)
                    .map_err(|error| match error {
                        Error::Io(error) => path_error(error, Some(&reader.path)),
                        Error::Invalid(message) => PyValueError::new_err(message),
                    })
            })?;
            state.lines.extend(batch);
        }
    }
}

#[pymethods]
impl QueryIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyString>>> {
        let line = match self
            .state
            .try_lock()
            .and_then(|mut state| state.lines.pop_front())
        {
            Some(line) => Some(line),
            None => py.detach(|| self.next_line())?,
        };
        Ok(line.map(|line| PyString::new(py, &line)))
    }
}

#[pymodule(gil_used = false)]
fn _pybgzf(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<Writer>()?;
    module.add_class::<Reader>()?;
    module.add_class::<IndexedReader>()?;
    module.add_class::<QueryIterator>()?;
    module.add_class::<Sniffer>()?;
    module.add_class::<IndexKind>()?;
    module.add_class::<Kind>()?;
    module.add_function(wrap_pyfunction!(validate_columns, module)?)?;
    module.add("BLOCK_SIZE", crate::block::BLOCK_SIZE)?;
    module.add(
        "TruncatedWarning",
        module.py().get_type::<TruncatedWarning>(),
    )?;
    Ok(())
}
