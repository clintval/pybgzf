//! Streaming BGZF compression with on-the-fly tabix and CSI indexing.

pub mod block;
pub mod columns;
pub mod index;
pub mod sniff;
pub mod writer;

#[cfg(feature = "python")]
mod python;
