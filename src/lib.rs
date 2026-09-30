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
