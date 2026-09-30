//! Streaming BGZF compression with on-the-fly tabix and CSI indexing.

pub mod block;
pub mod columns;
pub mod index;
pub mod khash;
pub mod reader;
pub mod sniff;
pub mod writer;

/// The names of the references in a file, in the order they are first seen.
pub type Names = indexmap::IndexSet<Vec<u8>, ahash::RandomState>;

#[cfg(feature = "python")]
mod python;
