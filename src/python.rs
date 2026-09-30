//! The `pybgzf._pybgzf` extension module.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use pyo3::buffer::PyBuffer;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

use crate::columns::{Columns, Kind};
use crate::sniff;
use crate::writer::{Error, IndexFormat, IndexOptions, Writer as CoreWriter, check_options};

const FILE_BUFFER: usize = 256 * 1024;

type ColumnsTuple = (i64, i64, Option<i64>, bool, String, i64, String);

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
    let kind = Kind::from_name(&kind)
        .ok_or_else(|| PyValueError::new_err(format!("unknown line format: {kind}")))?;
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
        columns.kind.name().to_string(),
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
    #[pyo3(signature = (dest, *, level, threads, index, index_path, columns, infer, csi_min_shift, csi_depth))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        dest: &Bound<'_, PyAny>,
        level: u8,
        threads: usize,
        index: Option<&str>,
        index_path: Option<PathBuf>,
        columns: Option<ColumnsTuple>,
        infer: bool,
        csi_min_shift: u32,
        csi_depth: Option<u32>,
    ) -> PyResult<Self> {
        let columns = columns.map(columns_from_tuple).transpose()?;
        let options = match index {
            None => {
                if index_path.is_some() {
                    return Err(PyValueError::new_err(
                        "index_path is only used when index is set",
                    ));
                }
                None
            }
            Some(format) => {
                let format = match format {
                    "tbi" => IndexFormat::Tabix,
                    "csi" => IndexFormat::Csi {
                        min_shift: csi_min_shift,
                        depth: csi_depth,
                    },
                    other => {
                        return Err(PyValueError::new_err(format!(
                            "unknown index format: {other}"
                        )));
                    }
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

#[pymodule]
fn _pybgzf(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<Writer>()?;
    module.add_class::<Sniffer>()?;
    module.add_function(wrap_pyfunction!(validate_columns, module)?)?;
    module.add("BLOCK_SIZE", crate::block::BLOCK_SIZE)?;
    Ok(())
}
