//! Reading BGZF sequentially or by region through a tabix or CSI index.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::num::NonZero;
use std::panic::{self, AssertUnwindSafe};
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use indexmap::IndexSet;
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::{MultithreadedReader, Reader as SerialReader, Seek as _};
use noodles_core::{Position, region::Interval};
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::Header;
use noodles_csi::binning_index::index::header::Format;
use noodles_csi::binning_index::index::header::format::CoordinateSystem;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::Index as ReferenceIndex;

use crate::Error;
use crate::columns::{Columns, Kind};
use crate::index::max_position;

const FILE_BUFFER: usize = 128 * 1024;

/// The empty block that ends every complete BGZF file.
const EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Returns whether a file ends with the BGZF end-of-file marker.
pub fn ends_with_eof_marker(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let Some(start) = file.metadata()?.len().checked_sub(EOF_MARKER.len() as u64) else {
        return Ok(false);
    };
    let mut tail = [0_u8; EOF_MARKER.len()];
    let _ = file.seek(SeekFrom::Start(start))?;
    file.read_exact(&mut tail)?;
    Ok(tail == EOF_MARKER)
}

/// Keeps reads off a noodles fast path that repeats the last block at the end of a stream
/// without an end-of-file marker.
const SERIAL_READ_LIMIT: usize = u16::MAX as usize;

/// Where a source is and its last bytes, shared with the thread reading it, so that bytes after
/// the last complete block, such as a truncated block header, and a missing end-of-file marker
/// are noticed.
#[derive(Default)]
struct Progress {
    offset: u64,
    tail: Vec<u8>,
}

fn lock(progress: &Mutex<Progress>) -> MutexGuard<'_, Progress> {
    progress.lock().unwrap_or_else(PoisonError::into_inner)
}

struct Tracked<R> {
    inner: R,
    progress: Arc<Mutex<Progress>>,
}

impl<R: Read> Read for Tracked<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        let mut progress = lock(&self.progress);
        progress.offset += n as u64;
        progress
            .tail
            .extend_from_slice(&buf[n.saturating_sub(EOF_MARKER.len())..n]);
        let excess = progress.tail.len().saturating_sub(EOF_MARKER.len());
        progress.tail.drain(..excess);
        Ok(n)
    }
}

impl<R: Seek> Seek for Tracked<R> {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        let offset = self.inner.seek(position)?;
        *lock(&self.progress) = Progress {
            offset,
            tail: Vec::new(),
        };
        Ok(offset)
    }
}

enum Inner<R: Read + Send + 'static> {
    Serial(SerialReader<Tracked<R>>),
    Parallel(MultithreadedReader<Tracked<R>>),
    Exhausted { source: Tracked<R>, position: u64 },
    Failed { position: u64 },
}

/// Decompresses BGZF in the calling thread or on worker threads, tracking virtual positions.
pub struct BgzfReader<R: Read + Seek + Send + 'static> {
    inner: Inner<R>,
    threads: NonZero<usize>,
    progress: Arc<Mutex<Progress>>,
    missing_eof_marker: Option<bool>,
}

fn parallel<R: Read + Send + 'static>(
    threads: NonZero<usize>,
    source: R,
) -> io::Result<MultithreadedReader<R>> {
    panic::catch_unwind(AssertUnwindSafe(|| {
        MultithreadedReader::with_worker_count(threads, source)
    }))
    .map_err(|_| io::Error::other("could not start the decompression threads"))
}

fn failed_earlier() -> io::Error {
    io::Error::other("the reader failed earlier")
}

fn not_in_file(position: u64) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("virtual offset {position} is not in this file"),
    )
}

impl BgzfReader<BufReader<File>> {
    /// Opens a BGZF file.
    pub fn from_path<P: AsRef<Path>>(path: P, threads: NonZero<usize>) -> io::Result<Self> {
        let file = BufReader::with_capacity(FILE_BUFFER, File::open(path)?);
        Self::new(file, threads)
    }
}

impl<R: Read + Seek + Send + 'static> BgzfReader<R> {
    /// Reads BGZF from `source`, decompressing on `threads` worker threads when more than one.
    pub fn new(source: R, threads: NonZero<usize>) -> io::Result<Self> {
        let progress = Arc::default();
        let source = Tracked {
            inner: source,
            progress: Arc::clone(&progress),
        };
        let mut reader = Self {
            inner: Inner::Failed { position: 0 },
            threads,
            progress,
            missing_eof_marker: None,
        };
        reader.start(source)?;
        Ok(reader)
    }

    fn start(&mut self, source: Tracked<R>) -> io::Result<()> {
        self.inner = if self.threads.get() == 1 {
            Inner::Serial(SerialReader::new(source))
        } else {
            Inner::Parallel(parallel(self.threads, source)?)
        };
        Ok(())
    }

    /// Returns the virtual position of the next byte to be read.
    pub fn virtual_position(&self) -> u64 {
        match &self.inner {
            Inner::Serial(reader) => u64::from(reader.virtual_position()),
            Inner::Parallel(reader) => u64::from(reader.virtual_position()),
            Inner::Exhausted { position, .. } | Inner::Failed { position } => *position,
        }
    }

    fn compressed_position(&self) -> u64 {
        match &self.inner {
            Inner::Serial(reader) => reader.position(),
            Inner::Parallel(reader) => reader.position(),
            Inner::Exhausted { position, .. } | Inner::Failed { position } => position >> 16,
        }
    }

    /// Moves to a virtual position, which must be the start of a line or record for the reads
    /// that follow to make sense, or the end of the file.
    ///
    /// Errors with [`io::ErrorKind::InvalidInput`] if the position is in no block of the file.
    pub fn seek(&mut self, position: u64) -> io::Result<()> {
        if let Inner::Exhausted { .. } = self.inner
            && let Inner::Exhausted { source, .. } =
                std::mem::replace(&mut self.inner, Inner::Failed { position: 0 })
        {
            self.start(source)?;
        }
        let target = VirtualPosition::from(position);
        let (block, offset) = target.into();
        let sought = match &mut self.inner {
            Inner::Serial(reader) => reader.seek(target).map(drop),
            Inner::Parallel(reader) => reader.seek_to_virtual_position(target).map(drop),
            Inner::Exhausted { .. } | Inner::Failed { .. } => return Err(failed_earlier()),
        };
        // Some systems refuse to seek a file past their largest offset rather than past its end.
        sought.map_err(|error| match error.kind() {
            io::ErrorKind::InvalidInput => not_in_file(position),
            _ => error,
        })?;
        if self.compressed_position() == block {
            let mut source = self.stop()?;
            let end = source.seek(SeekFrom::End(0))?;
            let at_end = block == end && offset == 0;
            let position = if at_end { position } else { end << 16 };
            self.inner = Inner::Exhausted { source, position };
            return if at_end {
                Ok(())
            } else {
                Err(not_in_file(position))
            };
        }
        if self.virtual_position() == position {
            return Ok(());
        }
        if offset == 0 && self.fill_buf()?.is_empty() {
            let source = self.stop()?;
            self.inner = Inner::Exhausted { source, position };
            return Ok(());
        }
        Err(not_in_file(position))
    }

    /// Stops any worker threads and returns the source, reporting any error they met.
    fn stop(&mut self) -> io::Result<Tracked<R>> {
        let position = self.virtual_position();
        match std::mem::replace(&mut self.inner, Inner::Failed { position }) {
            Inner::Serial(reader) => Ok(reader.into_inner()),
            Inner::Parallel(mut reader) => reader.finish(),
            Inner::Exhausted { source, .. } => Ok(source),
            Inner::Failed { .. } => Err(failed_earlier()),
        }
    }

    /// Returns true the first time the stream is found to end without the BGZF end-of-file
    /// marker, which means it may be truncated at a block boundary.
    pub fn take_missing_eof_marker(&mut self) -> bool {
        if self.missing_eof_marker == Some(true) {
            self.missing_eof_marker = Some(false);
            return true;
        }
        false
    }

    /// Checks at the end of the stream that no bytes follow the last complete block, then, when
    /// reading on threads, stops them, which reports any error they met, such as a corrupt block.
    fn exhaust(&mut self) -> io::Result<()> {
        let consumed = self.compressed_position();
        let (read, marked) = {
            let progress = lock(&self.progress);
            let marked = progress.tail.is_empty() || EOF_MARKER.ends_with(&progress.tail);
            (progress.offset, marked)
        };
        if read > consumed {
            let position = self.virtual_position();
            let parallel = matches!(self.inner, Inner::Parallel(_));
            if parallel {
                self.stop()?;
            }
            self.inner = Inner::Failed { position };
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "the BGZF file is truncated: {} bytes follow the last complete block",
                    read - consumed
                ),
            ));
        }
        if let Inner::Parallel(_) = self.inner {
            let position = self.virtual_position();
            let source = self.stop()?;
            self.inner = Inner::Exhausted { source, position };
        }
        if self.missing_eof_marker.is_none() && !marked {
            self.missing_eof_marker = Some(true);
        }
        Ok(())
    }
}

impl<R: Read + Seek + Send + 'static> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if let Inner::Serial(reader) = &mut self.inner {
            let n = buf.len().min(SERIAL_READ_LIMIT);
            let n = reader.read(&mut buf[..n])?;
            if n > 0 || buf.is_empty() {
                return Ok(n);
            }
            self.exhaust()?;
            return Ok(0);
        }
        let data = self.fill_buf()?;
        let n = data.len().min(buf.len());
        buf[..n].copy_from_slice(&data[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl<R: Read + Seek + Send + 'static> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let empty = match &mut self.inner {
            Inner::Serial(reader) => reader.fill_buf()?.is_empty(),
            Inner::Parallel(reader) => reader.fill_buf()?.is_empty(),
            Inner::Exhausted { .. } => false,
            Inner::Failed { .. } => return Err(failed_earlier()),
        };
        if empty {
            self.exhaust()?;
        }
        match &mut self.inner {
            Inner::Serial(reader) => reader.fill_buf(),
            Inner::Parallel(reader) => reader.fill_buf(),
            Inner::Exhausted { .. } => Ok(&[]),
            Inner::Failed { .. } => Err(failed_earlier()),
        }
    }

    fn consume(&mut self, amount: usize) {
        match &mut self.inner {
            Inner::Serial(reader) => reader.consume(amount),
            Inner::Parallel(reader) => reader.consume(amount),
            Inner::Exhausted { .. } | Inner::Failed { .. } => {}
        }
    }
}

/// A tabix or CSI index read from disk.
pub enum AnyIndex {
    Tabix(noodles_tabix::Index),
    Csi(noodles_csi::Index),
}

impl AnyIndex {
    /// Reads a BGZF-compressed tabix or CSI index, recognized by its magic number.
    pub fn read<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let invalid = || io::Error::new(io::ErrorKind::InvalidData, "not a tabix or CSI index");
        let mut magic = [0_u8; 4];
        noodles_bgzf::io::Reader::new(File::open(&path)?)
            .read_exact(&mut magic)
            .map_err(|_| invalid())?;
        match &magic {
            b"TBI\x01" => Ok(AnyIndex::Tabix(noodles_tabix::fs::read(path)?)),
            b"CSI\x01" => Ok(AnyIndex::Csi(noodles_csi::fs::read(path)?)),
            _ => Err(invalid()),
        }
    }

    fn binning(&self) -> &dyn BinningIndex {
        match self {
            AnyIndex::Tabix(index) => index,
            AnyIndex::Csi(index) => index,
        }
    }

    fn chunks(&self, tid: usize, interval: Interval) -> io::Result<Vec<Chunk>> {
        match self {
            AnyIndex::Tabix(index) => clipped_chunks(index, tid, interval),
            AnyIndex::Csi(index) => clipped_chunks(index, tid, interval),
        }
    }
}

/// Returns the chunks overlapping `interval`, skipping what precedes the first record that could
/// overlap it, as recorded in the linear or binned index, as htslib does.
fn clipped_chunks<I: ReferenceIndex>(
    index: &noodles_csi::binning_index::Index<I>,
    tid: usize,
    interval: Interval,
) -> io::Result<Vec<Chunk>> {
    let chunks = index.query(tid, interval)?;
    let Some(reference) = index.reference_sequences().get(tid) else {
        return Ok(chunks);
    };
    let start = interval.start().unwrap_or(Position::MIN);
    let min = reference.min_offset(index.min_shift(), index.depth(), start);
    Ok(chunks
        .into_iter()
        .filter(|chunk| chunk.end() > min)
        .map(|chunk| Chunk::new(chunk.start().max(min), chunk.end()))
        .collect())
}

/// Returns the columns an index header describes.
pub fn header_columns(header: &Header) -> Columns {
    let (kind, zero_based) = match header.format() {
        Format::Generic(CoordinateSystem::Bed) => (Kind::Generic, true),
        Format::Generic(CoordinateSystem::Gff) => (Kind::Generic, false),
        Format::Sam => (Kind::Sam, false),
        Format::Vcf => (Kind::Vcf, false),
    };
    let start = header.start_position_index() + 1;
    let end = match kind {
        Kind::Generic => header
            .end_position_index()
            .map(|end| end + 1)
            .filter(|&end| end != start),
        Kind::Sam | Kind::Vcf => None,
    };
    Columns {
        refname: header.reference_sequence_name_index() + 1,
        start,
        end,
        zero_based,
        meta_char: header.line_comment_prefix(),
        skip_lines: u64::from(header.line_skip_count()),
        kind,
    }
}

/// Where a region query is in its list of chunks.
pub struct Query {
    tid: usize,
    beg: i64,
    end: i64,
    chunks: Vec<Chunk>,
    next_chunk: usize,
    chunk_end: Option<u64>,
    position: Option<u64>,
    failure: Option<String>,
    /// True once every line has been returned.
    pub done: bool,
}

/// A BGZF file with its tabix or CSI index, for reading the lines that overlap a region.
pub struct IndexedReader<R: Read + Seek + Send + 'static> {
    reader: BgzfReader<R>,
    index: AnyIndex,
    columns: Columns,
    names: IndexSet<Vec<u8>>,
    line: Vec<u8>,
}

impl<R: Read + Seek + Send + 'static> IndexedReader<R> {
    /// Pairs a reader with an index that has a tabix header.
    pub fn new(reader: BgzfReader<R>, index: AnyIndex) -> io::Result<Self> {
        let header = index.binning().header().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "the index has no tabix header describing its columns",
            )
        })?;
        let columns = header_columns(header);
        let names = header
            .reference_sequence_names()
            .iter()
            .map(|name| name.to_vec())
            .collect();
        Ok(Self {
            reader,
            index,
            columns,
            names,
            line: Vec::new(),
        })
    }

    /// Returns the reference names in the index, in order.
    pub fn names(&self) -> impl Iterator<Item = &[u8]> {
        self.names.iter().map(Vec::as_slice)
    }

    /// Returns the columns recorded in the index header.
    pub fn columns(&self) -> &Columns {
        &self.columns
    }

    /// Starts a query for lines overlapping the 0-based, half-open `[beg, end)` on `name`.
    pub fn query(&self, name: &[u8], beg: i64, end: i64) -> io::Result<Query> {
        let index = self.index.binning();
        let limit = max_position(u32::from(index.min_shift()), u32::from(index.depth()));
        let mut query = Query {
            tid: 0,
            beg,
            end,
            chunks: Vec::new(),
            next_chunk: 0,
            chunk_end: None,
            position: None,
            failure: None,
            done: true,
        };
        let Some(tid) = self.names.get_index_of(name) else {
            return Ok(query);
        };
        if beg >= end {
            return Ok(query);
        }
        let to_position = |value: i64| {
            Position::new(usize::try_from(value).map_err(io::Error::other)?)
                .ok_or_else(|| io::Error::other("positions are 1-based"))
        };
        let (first, last) = (beg.min(limit - 2) + 1, end.min(limit - 1));
        let interval = Interval::from(to_position(first)?..=to_position(last)?);
        query.tid = tid;
        query.chunks = self.index.chunks(tid, interval)?;
        query.done = query.chunks.is_empty();
        Ok(query)
    }

    /// Reads the next lines of `query`, without their line terminators, appending up to about
    /// `budget` bytes of them to `lines`.
    ///
    /// A line that cannot be parsed, or that overlaps the region but is not UTF-8, is an error
    /// returned once the lines before it have been.
    pub fn next_lines(
        &mut self,
        query: &mut Query,
        lines: &mut Vec<String>,
        budget: usize,
    ) -> crate::Result<()> {
        if let Some(message) = query.failure.take() {
            query.done = true;
            return Err(Error::Invalid(message));
        }
        let mut used = 0;
        while !query.done && used < budget {
            if let Some(position) = query.position
                && position != self.reader.virtual_position()
            {
                self.reader.seek(position)?;
            }
            let at_chunk_end = match (query.position, query.chunk_end) {
                (Some(position), Some(end)) => position >= end,
                _ => true,
            };
            if at_chunk_end {
                let Some(chunk) = query.chunks.get(query.next_chunk) else {
                    query.done = true;
                    break;
                };
                let start = u64::from(chunk.start());
                if query.position != Some(start) {
                    self.reader.seek(start)?;
                }
                query.chunk_end = Some(u64::from(chunk.end()));
                query.next_chunk += 1;
            }
            let offset = self.reader.virtual_position();
            self.line.clear();
            if self.reader.read_until(b'\n', &mut self.line)? == 0 {
                query.done = true;
                break;
            }
            query.position = Some(self.reader.virtual_position());
            let content = trim(&self.line);
            if content.first() == Some(&self.columns.meta_char) {
                continue;
            }
            let interval = match self.columns.parse(content) {
                Ok(interval) => interval,
                Err(message) => {
                    let line = String::from_utf8_lossy(content);
                    query.failure = Some(format!(
                        "the line at virtual offset {offset} cannot be parsed: {message}: {line:?}"
                    ));
                    break;
                }
            };
            if self.names.get_index(query.tid).map(Vec::as_slice) != Some(interval.name)
                || interval.beg >= query.end
            {
                query.done = true;
                break;
            }
            if interval.end > query.beg && query.end > interval.beg {
                let Ok(text) = std::str::from_utf8(content) else {
                    query.failure = Some(format!(
                        "the line on {:?} at virtual offset {offset} is not UTF-8",
                        String::from_utf8_lossy(interval.name)
                    ));
                    break;
                };
                used += text.len();
                lines.push(text.to_owned());
            }
        }
        Ok(())
    }
}

fn trim(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;
    use std::io::Cursor;

    use super::*;
    use crate::writer::{IndexFormat, IndexOptions, Writer};

    fn threads(n: usize) -> NonZero<usize> {
        NonZero::new(n).unwrap()
    }

    fn bed(lines: usize) -> Vec<u8> {
        let mut text = String::new();
        for name in ["chr1", "chr2"] {
            for i in 0..lines {
                let start = i * 37;
                writeln!(
                    text,
                    "{name}\t{start}\t{}\tfeature{i}",
                    start + 1 + (i % 500)
                )
                .unwrap();
            }
        }
        text.into_bytes()
    }

    fn written(
        data: &[u8],
        format: IndexFormat,
    ) -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.bed.gz");
        let index_path = dir.path().join("x.bed.gz.idx");
        let options = IndexOptions {
            format,
            path: index_path.clone(),
            columns: Some(Columns::bed()),
            bed_only: false,
        };
        let mut writer = Writer::new(File::create(&path).unwrap(), 6, 2, Some(options)).unwrap();
        writer.write(data).unwrap();
        writer.finish().unwrap();
        (dir, path, index_path)
    }

    #[test]
    fn reads_back_what_was_written_with_any_thread_count() {
        let data = bed(20_000);
        let (_dir, path, _) = written(&data, IndexFormat::Tabix);
        for n in [1, 3] {
            let mut reader = BgzfReader::from_path(&path, threads(n)).unwrap();
            let mut out = Vec::new();
            reader.read_to_end(&mut out).unwrap();
            assert_eq!(out, data);
        }
    }

    #[test]
    fn a_missing_end_of_file_marker_is_reported_once_with_any_thread_count() {
        let data = bed(20_000);
        let (_dir, path, _) = written(&data, IndexFormat::Tabix);
        let compressed = std::fs::read(&path).unwrap();
        let cut = compressed[..compressed.len() - EOF_MARKER.len()].to_vec();
        for n in [1, 3] {
            let mut complete =
                BgzfReader::new(Cursor::new(compressed.clone()), threads(n)).unwrap();
            let _ = complete.read_to_end(&mut Vec::new()).unwrap();
            assert!(!complete.take_missing_eof_marker(), "threads={n}");
            let mut reader = BgzfReader::new(Cursor::new(cut.clone()), threads(n)).unwrap();
            let _ = reader.read_to_end(&mut Vec::new()).unwrap();
            assert!(reader.take_missing_eof_marker(), "threads={n}");
            assert!(!reader.take_missing_eof_marker(), "threads={n}");
        }
        assert!(ends_with_eof_marker(&path).unwrap());
    }

    #[test]
    fn truncated_input_is_an_error_with_any_thread_count() {
        let data = bed(20_000);
        let (_dir, path, _) = written(&data, IndexFormat::Tabix);
        let compressed = std::fs::read(&path).unwrap();
        let truncated = compressed[..compressed.len() / 2].to_vec();
        for n in [1, 3] {
            let mut reader = BgzfReader::new(Cursor::new(truncated.clone()), threads(n)).unwrap();
            assert!(reader.read_to_end(&mut Vec::new()).is_err(), "threads={n}");
        }
    }

    #[test]
    fn seeks_to_virtual_positions() {
        let data = bed(20_000);
        let (_dir, path, _) = written(&data, IndexFormat::Tabix);
        for n in [1, 2] {
            let mut reader = BgzfReader::from_path(&path, threads(n)).unwrap();
            let mut marks = Vec::new();
            let mut line = Vec::new();
            for _ in 0..30_000 {
                marks.push((reader.virtual_position(), {
                    line.clear();
                    reader.read_until(b'\n', &mut line).unwrap();
                    line.clone()
                }));
            }
            for (position, expected) in marks.iter().rev().step_by(997) {
                reader.seek(*position).unwrap();
                line.clear();
                reader.read_until(b'\n', &mut line).unwrap();
                assert_eq!(&line, expected);
            }
        }
    }

    fn brute_force(data: &[u8], name: &[u8], beg: i64, end: i64) -> Vec<String> {
        data.split(|&b| b == b'\n')
            .filter(|line| !line.is_empty())
            .filter(|line| {
                let parsed = Columns::bed().parse(line).unwrap();
                parsed.name == name && parsed.beg < end && parsed.end > beg && beg < end
            })
            .map(|line| String::from_utf8(line.to_vec()).unwrap())
            .collect()
    }

    #[test]
    fn queries_match_brute_force() {
        let data = bed(20_000);
        for format in [
            IndexFormat::Tabix,
            IndexFormat::Csi {
                min_shift: 14,
                depth: None,
            },
        ] {
            let (_dir, path, index_path) = written(&data, format);
            for n in [1, 2] {
                let reader = BgzfReader::from_path(&path, threads(n)).unwrap();
                let mut indexed =
                    IndexedReader::new(reader, AnyIndex::read(&index_path).unwrap()).unwrap();
                for (name, beg, end) in [
                    (&b"chr1"[..], 0, 1),
                    (b"chr1", 5000, 90_000),
                    (b"chr2", 739_000, 740_000),
                    (b"chr2", 100, 100),
                    (b"chr2", 10_000_000, 10_000_100),
                    (b"chr3", 0, 100),
                ] {
                    let mut query = indexed.query(name, beg, end).unwrap();
                    let mut lines = Vec::new();
                    while !query.done {
                        indexed.next_lines(&mut query, &mut lines, 1000).unwrap();
                    }
                    assert_eq!(
                        lines,
                        brute_force(&data, name, beg, end),
                        "{format:?} {n} {beg}-{end}"
                    );
                }
            }
        }
    }

    #[test]
    fn queries_start_at_the_linear_index() {
        let data = bed(20_000);
        let (_dir, path, index_path) = written(&data, IndexFormat::Tabix);
        let reader = BgzfReader::from_path(&path, threads(1)).unwrap();
        let indexed = IndexedReader::new(reader, AnyIndex::read(&index_path).unwrap()).unwrap();
        let query = indexed.query(b"chr1", 700_000, 700_100).unwrap();
        let first = u64::from(query.chunks[0].start());
        assert!(first >> 16 > 0, "the query starts at {first}");
    }

    #[test]
    fn interleaved_queries_keep_their_place() {
        let data = bed(5000);
        let (_dir, path, index_path) = written(&data, IndexFormat::Tabix);
        let reader = BgzfReader::from_path(&path, threads(1)).unwrap();
        let mut indexed = IndexedReader::new(reader, AnyIndex::read(&index_path).unwrap()).unwrap();
        let mut first = indexed.query(b"chr1", 0, 100_000).unwrap();
        let mut second = indexed.query(b"chr2", 0, 100_000).unwrap();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        loop {
            indexed.next_lines(&mut first, &mut a, 10).unwrap();
            indexed.next_lines(&mut second, &mut b, 10).unwrap();
            if first.done && second.done {
                break;
            }
        }
        assert_eq!(a, brute_force(&data, b"chr1", 0, 100_000));
        assert_eq!(b, brute_force(&data, b"chr2", 0, 100_000));
    }

    #[test]
    fn header_columns_round_trip() {
        let data = bed(10);
        let (_dir, _path, index_path) = written(&data, IndexFormat::Tabix);
        let index = AnyIndex::read(&index_path).unwrap();
        assert_eq!(
            header_columns(index.binning().header().unwrap()),
            Columns::bed()
        );
    }

    #[test]
    fn rejects_files_that_are_not_indexes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.gz");
        let mut writer = Writer::new(File::create(&path).unwrap(), 6, 1, None).unwrap();
        writer.write(b"not an index").unwrap();
        writer.finish().unwrap();
        assert!(AnyIndex::read(&path).is_err());
        assert!(
            BgzfReader::new(Cursor::new(Vec::new()), threads(1))
                .unwrap()
                .read_until(b'\n', &mut Vec::new())
                .is_ok()
        );
    }
}
