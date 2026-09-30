//! Streaming BGZF compression with on-the-fly tabix and CSI indexing.

pub mod block;
pub mod columns;
pub mod index;
pub mod khash;
pub mod reader;
pub mod sniff;
pub mod writer;

#[cfg(feature = "python")]
mod python;

use std::io;
use std::num::NonZero;

/// An error: either I/O failed or the data is invalid, such as a line that cannot be indexed.
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    Invalid(String),
}

impl From<io::Error> for Error {
    fn from(error: io::Error) -> Self {
        Error::Io(error)
    }
}

pub type Result<T> = std::result::Result<T, Error>;

/// The most threads a reader or writer may use.
pub const MAX_THREADS: usize = 1024;

/// Checks that a thread count is between 1 and [`MAX_THREADS`].
pub fn check_threads(threads: i64) -> std::result::Result<NonZero<usize>, String> {
    usize::try_from(threads)
        .ok()
        .and_then(NonZero::new)
        .filter(|threads| threads.get() <= MAX_THREADS)
        .ok_or_else(|| format!("threads must be between 1 and {MAX_THREADS}, not {threads}"))
}
