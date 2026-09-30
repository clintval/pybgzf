//! The `pybgzf._pybgzf` extension module.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::thread;

use pyo3::buffer::PyBuffer;
use pyo3::exceptions::{PyOSError, PyUserWarning, PyValueError};
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

type ColumnsTuple = (i64, i64, Option<i64>, bool, String, i64, LineKind);

/// The value of each `pybgzf.LineFormat` member.
#[pyclass(module = "pybgzf._pybgzf", eq, frozen, hash, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum LineKind {
    #[pyo3(name = "GENERIC")]
    Generic,
    #[pyo3(name = "SAM")]
    Sam,
    #[pyo3(name = "VCF")]
    Vcf,
}

impl From<LineKind> for Kind {
    fn from(kind: LineKind) -> Self {
        match kind {
            LineKind::Generic => Kind::Generic,
            LineKind::Sam => Kind::Sam,
            LineKind::Vcf => Kind::Vcf,
        }
    }
}

impl From<Kind> for LineKind {
    fn from(kind: Kind) -> Self {
        match kind {
            Kind::Generic => LineKind::Generic,
            Kind::Sam => LineKind::Sam,
            Kind::Vcf => LineKind::Vcf,
        }
    }
}

/// The value of each `pybgzf.IndexFormat` member.
#[pyclass(module = "pybgzf._pybgzf", eq, frozen, hash, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum IndexKind {
    #[pyo3(name = "TBI")]
    Tabix,
    #[pyo3(name = "CSI")]
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
    let kind = Kind::from(kind);
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
        LineKind::from(columns.kind),
    )
}

/// The Rust half of `pybgzf.BgzfWriter`.
#[pyclass(module = "pybgzf._pybgzf")]
struct Writer {
    inner: CoreWriter<Sink>,
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
        let options = match index {
            None => {
                if index_path.is_some() {
                    return Err(PyValueError::new_err(
                        "index_path is only used when index is set",
                    ));
                }
                None
            }
            Some(kind) => {
                let format = match kind {
                    IndexKind::Tabix => IndexFormat::Tabix,
                    IndexKind::Csi => csi_format(csi_min_shift, csi_depth).map_err(to_python)?,
                };
                if columns.is_none() && !infer {
                    return Err(PyValueError::new_err(
                        "columns is required when index is set",
                    ));
                }
                let path = index_path.ok_or_else(|| {
                    PyValueError::new_err("index_path is required when index is set")
                })?;
                Some(IndexOptions {
                    format,
                    path,
                    columns: if infer { None } else { columns },
                    bed_only: infer_bed,
                })
            }
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
        Ok(Self { inner })
    }

    /// Writes bytes-like `data` and returns the number of bytes written.
    fn write(&mut self, py: Python<'_>, data: &Bound<'_, PyAny>) -> PyResult<usize> {
        let owned;
        let bytes: &[u8] = if let Ok(bytes) = data.cast::<PyBytes>() {
            bytes.as_bytes()
        } else {
            owned = PyBuffer::<u8>::get(data)?.to_vec(py)?;
            &owned
        };
        let inner = &mut self.inner;
        if inner.may_block(bytes.len()) {
            py.detach(|| inner.write(bytes)).map_err(to_python)?;
        } else {
            inner.write(bytes).map_err(to_python)?;
        }
        Ok(bytes.len())
    }

    /// Ends the current block and flushes the sink.
    fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        let inner = &mut self.inner;
        py.detach(|| inner.flush()).map_err(to_python)
    }

    /// Returns the virtual position of the next byte.
    fn tell(&mut self, py: Python<'_>) -> PyResult<u64> {
        let inner = &mut self.inner;
        py.detach(|| inner.tell()).map_err(to_python)
    }

    /// Writes the end-of-file marker and the index, if any.
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        let inner = &mut self.inner;
        let finished = py.detach(|| {
            let finished = inner.finish();
            *inner.get_mut() = Sink::Closed;
            finished
        });
        finished.map_err(to_python)
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.is_finished()
    }

    #[getter]
    fn columns(&self) -> Option<ColumnsTuple> {
        self.inner.columns().map(columns_to_tuple)
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

fn closed_error() -> PyErr {
    PyValueError::new_err("I/O operation on closed file.")
}

fn drop_detached<T: Send>(py: Python<'_>, value: T) {
    py.detach(move || drop(value));
}

/// Stops a reader's threads, in the background when its source may block, such as a pipe whose
/// writer is idle, so that closing never waits on a read that may not return.
fn release(py: Python<'_>, reader: BgzfReader<Source>, may_block: bool) {
    if may_block && reader.reads_ahead() {
        py.detach(|| drop(thread::Builder::new().spawn(move || drop(reader))));
    } else {
        drop_detached(py, reader);
    }
}

/// The Rust half of `pybgzf.BgzfReader`.
#[pyclass(module = "pybgzf._pybgzf")]
struct Reader {
    inner: Option<BgzfReader<Source>>,
    path: Option<PathBuf>,
    may_block: bool,
}

impl Reader {
    fn inner(&mut self) -> PyResult<(&mut BgzfReader<Source>, Option<&Path>)> {
        let inner = self.inner.as_mut().ok_or_else(closed_error)?;
        Ok((inner, self.path.as_deref()))
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
        Ok(Self {
            inner: Some(inner),
            path,
            may_block,
        })
    }

    /// Reads up to `size` bytes, fewer only at the end of the stream.
    fn read<'py>(&mut self, py: Python<'py>, size: usize) -> PyResult<Bound<'py, PyBytes>> {
        let (inner, path) = self.inner()?;
        let mut buf = Vec::with_capacity(size);
        py.detach(|| (&mut *inner).take(size as u64).read_to_end(&mut buf))
            .map_err(|e| path_error(e, path))?;
        if inner.take_missing_eof_marker() {
            warn_truncated(py, path)?;
        }
        Ok(PyBytes::new(py, &buf))
    }

    /// Reads everything that is left.
    fn readall<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let (inner, path) = self.inner()?;
        let mut buf = Vec::new();
        py.detach(|| inner.read_to_end(&mut buf))
            .map_err(|e| path_error(e, path))?;
        if inner.take_missing_eof_marker() {
            warn_truncated(py, path)?;
        }
        Ok(PyBytes::new(py, &buf))
    }

    /// Reads through the next newline.
    fn readline<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let (inner, path) = self.inner()?;
        let mut line = Vec::new();
        py.detach(|| inner.read_until(b'\n', &mut line))
            .map_err(|e| path_error(e, path))?;
        if inner.take_missing_eof_marker() {
            warn_truncated(py, path)?;
        }
        Ok(PyBytes::new(py, &line))
    }

    /// Returns the virtual position of the next byte.
    fn tell(&mut self) -> PyResult<u64> {
        Ok(self.inner()?.0.virtual_position())
    }

    /// Moves to a virtual position.
    fn seek(&mut self, py: Python<'_>, position: u64) -> PyResult<u64> {
        let (inner, path) = self.inner()?;
        py.detach(|| inner.seek(position))
            .map_err(|error| match error.kind() {
                io::ErrorKind::InvalidInput => PyValueError::new_err(error.to_string()),
                _ => path_error(error, path),
            })?;
        Ok(inner.virtual_position())
    }

    /// Stops any worker threads.
    fn close(&mut self, py: Python<'_>) {
        if let Some(inner) = self.inner.take() {
            release(py, inner, self.may_block);
        }
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.is_none()
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            Python::attach(|py| release(py, inner, self.may_block));
        }
    }
}

/// The Rust half of `pybgzf.IndexedReader`.
#[pyclass(module = "pybgzf._pybgzf")]
struct IndexedReader {
    inner: Option<CoreIndexedReader<BufReader<File>>>,
    path: PathBuf,
}

impl IndexedReader {
    fn inner(&mut self) -> PyResult<&mut CoreIndexedReader<BufReader<File>>> {
        self.inner.as_mut().ok_or_else(closed_error)
    }
}

#[pymethods]
impl IndexedReader {
    #[new]
    #[pyo3(signature = (path, index_path, *, threads))]
    #[allow(clippy::needless_pass_by_value)]
    fn new(py: Python<'_>, path: PathBuf, index_path: PathBuf, threads: i64) -> PyResult<Self> {
        let threads = self::threads(threads)?;
        let index =
            py.detach(|| AnyIndex::read(&index_path))
                .map_err(|error| match error.kind() {
                    io::ErrorKind::InvalidData => {
                        PyValueError::new_err(format!("{}: {error}", index_path.display()))
                    }
                    _ => path_error(error, Some(&index_path)),
                })?;
        let inner = py
            .detach(|| CoreIndexedReader::new(BgzfReader::from_path(&path, threads)?, index))
            .map_err(|error| match error.kind() {
                io::ErrorKind::InvalidData => {
                    PyValueError::new_err(format!("{}: {error}", index_path.display()))
                }
                _ => path_error(error, Some(&path)),
            })?;
        if !py
            .detach(|| ends_with_eof_marker(&path))
            .map_err(|error| path_error(error, Some(&path)))?
        {
            warn_truncated(py, Some(&path))?;
        }
        Ok(Self {
            inner: Some(inner),
            path,
        })
    }

    /// Starts a query for lines overlapping the 0-based, half-open `[start, end)` on `refname`.
    fn query(slf: Bound<'_, Self>, refname: &str, start: i64, end: i64) -> PyResult<QueryIterator> {
        if start < 0 || end < start {
            return Err(PyValueError::new_err(format!(
                "start must be at least 0 and end at least start, not {start} and {end}"
            )));
        }
        let query = slf
            .borrow_mut()
            .inner()?
            .query(refname.as_bytes(), start, end)
            .map_err(io_to_python)?;
        Ok(QueryIterator {
            reader: slf.unbind(),
            query,
            lines: VecDeque::new(),
        })
    }

    #[getter]
    fn refnames(&mut self) -> PyResult<Vec<String>> {
        self.inner()?
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
    }

    #[getter]
    fn columns(&mut self) -> PyResult<ColumnsTuple> {
        Ok(columns_to_tuple(self.inner()?.columns()))
    }

    /// Stops any worker threads.
    fn close(&mut self, py: Python<'_>) {
        if let Some(inner) = self.inner.take() {
            drop_detached(py, inner);
        }
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.is_none()
    }
}

impl Drop for IndexedReader {
    fn drop(&mut self) {
        if let Some(inner) = self.inner.take() {
            Python::attach(|py| drop_detached(py, inner));
        }
    }
}

const QUERY_BATCH: usize = 256 * 1024;

/// Lines overlapping a region, read in batches with the GIL released.
#[pyclass(module = "pybgzf._pybgzf")]
struct QueryIterator {
    reader: Py<IndexedReader>,
    query: Query,
    lines: VecDeque<String>,
}

#[pymethods]
impl QueryIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyString>>> {
        loop {
            if let Some(line) = self.lines.pop_front() {
                return Ok(Some(PyString::new(py, &line)));
            }
            if self.query.done {
                return Ok(None);
            }
            let mut reader = self.reader.bind(py).borrow_mut();
            let inner = reader.inner()?;
            let query = &mut self.query;
            let mut batch = Vec::new();
            py.detach(|| inner.next_lines(query, &mut batch, QUERY_BATCH))
                .map_err(|error| match error {
                    Error::Io(error) => path_error(error, Some(&reader.path)),
                    Error::Invalid(message) => PyValueError::new_err(message),
                })?;
            self.lines.extend(batch);
        }
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
    module.add_class::<LineKind>()?;
    module.add_function(wrap_pyfunction!(validate_columns, module)?)?;
    module.add("BLOCK_SIZE", crate::block::BLOCK_SIZE)?;
    module.add(
        "TruncatedWarning",
        module.py().get_type::<TruncatedWarning>(),
    )?;
    Ok(())
}
