//! Reading BGZF sequentially or by region through a tabix or CSI index.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};
use std::num::NonZero;
use std::path::Path;

use indexmap::IndexSet;
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::{MultithreadedReader, Reader as SerialReader, Seek as _};
use noodles_core::{Position, region::Interval};
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::Header;
use noodles_csi::binning_index::index::header::Format;
use noodles_csi::binning_index::index::header::format::CoordinateSystem;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::columns::{Columns, Kind};
use crate::index::max_position;

const FILE_BUFFER: usize = 128 * 1024;

enum Inner<R: Read + Send + 'static> {
    Serial(SerialReader<R>),
    Parallel(MultithreadedReader<R>),
}

/// Decompresses BGZF in the calling thread or on worker threads, tracking virtual positions.
pub struct BgzfReader<R: Read + Seek + Send + 'static> {
    inner: Inner<R>,
}

impl BgzfReader<BufReader<File>> {
    /// Opens a BGZF file.
    pub fn from_path<P: AsRef<Path>>(path: P, threads: NonZero<usize>) -> io::Result<Self> {
        let file = BufReader::with_capacity(FILE_BUFFER, File::open(path)?);
        Ok(Self::new(file, threads))
    }
}

impl<R: Read + Seek + Send + 'static> BgzfReader<R> {
    /// Reads BGZF from `source`, decompressing on `threads` worker threads when more than one.
    pub fn new(source: R, threads: NonZero<usize>) -> Self {
        let inner = if threads.get() == 1 {
            Inner::Serial(SerialReader::new(source))
        } else {
            Inner::Parallel(MultithreadedReader::with_worker_count(threads, source))
        };
        Self { inner }
    }

    /// Returns the virtual position of the next byte to be read.
    pub fn virtual_position(&self) -> u64 {
        u64::from(match &self.inner {
            Inner::Serial(reader) => reader.virtual_position(),
            Inner::Parallel(reader) => reader.virtual_position(),
        })
    }

    /// Moves to a virtual position, which must be the start of a line or record for the reads
    /// that follow to make sense.
    pub fn seek(&mut self, position: u64) -> io::Result<()> {
        let position = VirtualPosition::from(position);
        match &mut self.inner {
            Inner::Serial(reader) => reader.seek(position).map(drop),
            Inner::Parallel(reader) => reader.seek_to_virtual_position(position).map(drop),
        }
    }

    /// Reads up to and including the next newline into `line`, returning the bytes read.
    pub fn read_line(&mut self, line: &mut Vec<u8>) -> io::Result<usize> {
        self.read_until(b'\n', line)
    }

    /// Fills `buf` as far as possible, returning fewer bytes only at the end of the stream.
    pub fn read_full(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut filled = 0;
        while filled < buf.len() {
            let n = self.read(&mut buf[filled..])?;
            if n == 0 {
                break;
            }
            filled += n;
        }
        Ok(filled)
    }

    /// Stops any worker threads and returns the source.
    pub fn finish(self) -> io::Result<R> {
        match self.inner {
            Inner::Serial(reader) => Ok(reader.into_inner()),
            Inner::Parallel(mut reader) => reader.finish(),
        }
    }
}

impl<R: Read + Seek + Send + 'static> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match &mut self.inner {
            Inner::Serial(reader) => reader.read(buf),
            Inner::Parallel(reader) => reader.read(buf),
        }
    }
}

impl<R: Read + Seek + Send + 'static> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match &mut self.inner {
            Inner::Serial(reader) => reader.fill_buf(),
            Inner::Parallel(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amount: usize) {
        match &mut self.inner {
            Inner::Serial(reader) => reader.consume(amount),
            Inner::Parallel(reader) => reader.consume(amount),
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

    fn header(&self) -> Option<&Header> {
        match self {
            AnyIndex::Tabix(index) => index.header(),
            AnyIndex::Csi(index) => index.header(),
        }
    }

    fn bins(&self) -> (u32, u32) {
        match self {
            AnyIndex::Tabix(index) => (u32::from(index.min_shift()), u32::from(index.depth())),
            AnyIndex::Csi(index) => (u32::from(index.min_shift()), u32::from(index.depth())),
        }
    }

    fn chunks(&self, tid: usize, interval: Interval) -> io::Result<Vec<Chunk>> {
        match self {
            AnyIndex::Tabix(index) => index.query(tid, interval),
            AnyIndex::Csi(index) => index.query(tid, interval),
        }
    }
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
    done: bool,
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
        let header = index.header().ok_or_else(|| {
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
        let (min_shift, depth) = self.index.bins();
        let limit = max_position(min_shift, depth);
        let mut query = Query {
            tid: 0,
            beg,
            end,
            chunks: Vec::new(),
            next_chunk: 0,
            chunk_end: None,
            position: None,
            done: true,
        };
        let Some(tid) = self.names.get_index_of(name) else {
            return Ok(query);
        };
        let clamped_end = end.min(limit - 1);
        if beg >= clamped_end {
            return Ok(query);
        }
        let to_position = |value: i64| {
            Position::new(usize::try_from(value).map_err(io::Error::other)?)
                .ok_or_else(|| io::Error::other("positions are 1-based"))
        };
        let interval = Interval::from(to_position(beg + 1)?..=to_position(clamped_end)?);
        query.tid = tid;
        query.chunks = self.index.chunks(tid, interval)?;
        query.done = query.chunks.is_empty();
        Ok(query)
    }

    /// Reads the next lines of `query`, without their line terminators, appending up to about
    /// `budget` bytes of them to `lines`. Returns false once the query is exhausted.
    pub fn next_lines(
        &mut self,
        query: &mut Query,
        lines: &mut Vec<Vec<u8>>,
        budget: usize,
    ) -> Result<bool, QueryError> {
        let mut used = 0;
        while !query.done && used < budget {
            if query
                .position
                .is_some_and(|position| position != self.reader.virtual_position())
            {
                self.reader.seek(query.position.expect("checked above"))?;
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
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                query.done = true;
                break;
            }
            query.position = Some(self.reader.virtual_position());
            let content = trim(&self.line);
            if content.first() == Some(&self.columns.meta_char) {
                continue;
            }
            let interval = self.columns.parse(content).map_err(|message| {
                QueryError::Invalid(format!("{message}: {:?}", String::from_utf8_lossy(content)))
            })?;
            if self.names.get_index_of(interval.name) != Some(query.tid)
                || interval.beg >= query.end
            {
                query.done = true;
                break;
            }
            if interval.end > query.beg && query.end > interval.beg {
                used += content.len();
                lines.push(content.to_vec());
            }
        }
        Ok(!query.done)
    }

    /// Stops any worker threads.
    pub fn finish(self) -> io::Result<R> {
        self.reader.finish()
    }
}

/// An error from a query: either I/O failed or a line could not be parsed.
#[derive(Debug)]
pub enum QueryError {
    Io(io::Error),
    Invalid(String),
}

impl From<io::Error> for QueryError {
    fn from(error: io::Error) -> Self {
        QueryError::Io(error)
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
                    reader.read_line(&mut line).unwrap();
                    line.clone()
                }));
            }
            for (position, expected) in marks.iter().rev().step_by(997) {
                reader.seek(*position).unwrap();
                line.clear();
                reader.read_line(&mut line).unwrap();
                assert_eq!(&line, expected);
            }
        }
    }

    fn brute_force(data: &[u8], name: &[u8], beg: i64, end: i64) -> Vec<Vec<u8>> {
        data.split(|&b| b == b'\n')
            .filter(|line| !line.is_empty())
            .filter(|line| {
                let parsed = Columns::bed().parse(line).unwrap();
                parsed.name == name && parsed.beg < end && parsed.end > beg && beg < end
            })
            .map(<[u8]>::to_vec)
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
                    while indexed.next_lines(&mut query, &mut lines, 1000).unwrap() {}
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
    fn interleaved_queries_keep_their_place() {
        let data = bed(5000);
        let (_dir, path, index_path) = written(&data, IndexFormat::Tabix);
        let reader = BgzfReader::from_path(&path, threads(1)).unwrap();
        let mut indexed = IndexedReader::new(reader, AnyIndex::read(&index_path).unwrap()).unwrap();
        let mut first = indexed.query(b"chr1", 0, 100_000).unwrap();
        let mut second = indexed.query(b"chr2", 0, 100_000).unwrap();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        loop {
            let more_a = indexed.next_lines(&mut first, &mut a, 10).unwrap();
            let more_b = indexed.next_lines(&mut second, &mut b, 10).unwrap();
            if !more_a && !more_b {
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
        assert_eq!(header_columns(index.header().unwrap()), Columns::bed());
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
                .read_line(&mut Vec::new())
                .is_ok()
        );
    }
}
