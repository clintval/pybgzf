//! A streaming BGZF writer that indexes lines as they are written.

use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::num::NonZero;
use std::path::PathBuf;

use bgzf::CompressionLevel;
use indexmap::IndexSet;
use memchr::{memchr, memmem};

use crate::block::{BLOCK_SIZE, BlockWriter, LogicalPosition};
use crate::columns::{Columns, Kind};
use crate::index::{self, IndexBuilder, Layout, max_position};
use crate::sniff::Sniffer;

const TABIX_MIN_SHIFT: u32 = 14;
const TABIX_DEPTH: u32 = 5;
const TABIX_MAX_SHIFT: u32 = 31;
const NO_COLUMNS: &str =
    "could not infer columns because no data lines were written; pass columns explicitly";

/// An error from writing: either I/O failed or the data cannot be indexed.
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

/// Which index to build, and how.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndexFormat {
    Tabix,
    Csi { min_shift: u32, depth: Option<u32> },
}

/// Everything needed to index lines while writing.
pub struct IndexOptions {
    pub format: IndexFormat,
    pub path: PathBuf,
    pub columns: Option<Columns>,
    pub bed_only: bool,
}

struct Record {
    tid: usize,
    beg: i64,
    end: i64,
    end_position: LogicalPosition,
}

struct Indexer {
    format: IndexFormat,
    path: PathBuf,
    file: Option<File>,
    columns: Option<Columns>,
    sniffer: Sniffer,
    partial: Vec<u8>,
    line_number: u64,
    names: IndexSet<Vec<u8>>,
    last: Option<(usize, i64, u64)>,
    bins: Option<(u32, u32)>,
    longest_vcf_contig: i64,
    longest_sam_reference: i64,
    first_offset: LogicalPosition,
    pending: VecDeque<Record>,
    builder: Option<IndexBuilder>,
    failure: Option<String>,
}

fn longest(line: &[u8], prefix: &[u8], key: &[u8], skip_padding: bool) -> Option<i64> {
    let rest = line.strip_prefix(prefix)?;
    let at = memmem::find(rest, key)?;
    let mut value = &rest[at + key.len()..];
    if skip_padding {
        let padding = value
            .iter()
            .take_while(|&&b| b == b' ' || b == b'=')
            .count();
        value = &value[padding..];
    }
    let digits = value.iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    Some(value[..digits].iter().fold(0_i64, |total, digit| {
        total
            .saturating_mul(10)
            .saturating_add(i64::from(digit - b'0'))
    }))
}

fn csi_bins(mut min_shift: u32, longest_reference: i64) -> std::result::Result<(u32, u32), String> {
    const MAX_DEPTH: u32 = 9;
    const MAX_SHIFT: u32 = 62;
    if longest_reference <= 0 {
        let depth = match min_shift {
            0..10 => MAX_DEPTH,
            10..25 => MAX_DEPTH - (min_shift - 10) / 3,
            _ => 4,
        };
        return Ok((min_shift, depth));
    }
    let needed = longest_reference.saturating_add(256);
    if needed > 1_i64 << MAX_SHIFT {
        return Err(format!(
            "a reference of length {longest_reference} in the header is too long for a CSI index, which holds up to {}",
            (1_i64 << MAX_SHIFT) - 256
        ));
    }
    let mut depth = (TABIX_MAX_SHIFT + 2).saturating_sub(min_shift) / 3;
    if needed <= max_position(min_shift, MAX_DEPTH) {
        while needed > max_position(min_shift, depth) {
            depth += 1;
        }
    } else {
        depth = MAX_DEPTH;
        while needed > max_position(min_shift, depth) {
            min_shift += 1;
        }
    }
    Ok((min_shift, depth))
}

fn trim_line(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

impl Indexer {
    fn create(options: IndexOptions) -> io::Result<Self> {
        let file = File::create(&options.path).map_err(|error| {
            let path = options.path.display();
            io::Error::new(
                error.kind(),
                format!("cannot create the index {path}: {error}"),
            )
        })?;
        Ok(Self {
            format: options.format,
            path: options.path,
            file: Some(file),
            columns: options.columns,
            sniffer: if options.bed_only {
                Sniffer::bed()
            } else {
                Sniffer::default()
            },
            partial: Vec::new(),
            line_number: 0,
            names: IndexSet::new(),
            last: None,
            bins: None,
            longest_vcf_contig: 0,
            longest_sam_reference: 0,
            first_offset: LogicalPosition {
                block: 0,
                offset: 0,
            },
            pending: VecDeque::new(),
            builder: None,
            failure: None,
        })
    }

    fn remove_file(&mut self) -> io::Result<()> {
        self.file.take();
        match fs::remove_file(&self.path) {
            Err(error) if error.kind() != io::ErrorKind::NotFound => Err(error),
            _ => Ok(()),
        }
    }

    fn note_meta(&mut self, line: &[u8]) {
        if let Some(length) = longest(line, b"##contig", b"length", true) {
            self.longest_vcf_contig = self.longest_vcf_contig.max(length);
        }
        if let Some(length) = longest(line, b"@SQ", b"\tLN:", false) {
            self.longest_sam_reference = self.longest_sam_reference.max(length);
        }
    }

    fn decide_bins(&self, columns: &Columns) -> std::result::Result<(u32, u32), String> {
        match self.format {
            IndexFormat::Tabix => Ok((TABIX_MIN_SHIFT, TABIX_DEPTH)),
            IndexFormat::Csi {
                min_shift,
                depth: Some(depth),
            } => Ok((min_shift, depth)),
            IndexFormat::Csi {
                min_shift,
                depth: None,
            } => {
                let longest_reference = match columns.kind {
                    Kind::Vcf => self.longest_vcf_contig,
                    Kind::Sam => self.longest_sam_reference,
                    Kind::Generic => 0,
                };
                csi_bins(min_shift, longest_reference)
            }
        }
    }

    fn empty_bins(&self) -> (u32, u32) {
        match self.format {
            IndexFormat::Tabix => (TABIX_MIN_SHIFT, TABIX_DEPTH),
            IndexFormat::Csi { min_shift, depth } => (
                min_shift,
                depth.unwrap_or((TABIX_MAX_SHIFT + 2).saturating_sub(min_shift) / 3),
            ),
        }
    }

    /// Classifies a line, without its terminator, as a header line or a record to index.
    fn classify(&mut self, line: &[u8]) -> std::result::Result<Option<(usize, i64, i64)>, String> {
        let number = self.line_number + 1;
        let columns = match &self.columns {
            Some(columns) => columns.clone(),
            None => {
                let Some(columns) = self.sniffer.push(line)? else {
                    self.line_number = number;
                    self.note_meta(line);
                    return Ok(None);
                };
                self.columns = Some(columns.clone());
                columns
            }
        };
        if columns.is_meta(number, line) {
            self.line_number = number;
            self.note_meta(line);
            return Ok(None);
        }
        let interval = columns
            .parse(line)
            .map_err(|e| format!("line {number}: {e}"))?;
        if memchr(0, interval.name).is_some() {
            return Err(format!(
                "line {number}: the reference name {:?} contains a NUL byte",
                String::from_utf8_lossy(interval.name)
            ));
        }
        let (min_shift, depth) = match self.bins {
            Some(bins) => bins,
            None => self
                .decide_bins(&columns)
                .map_err(|e| format!("line {number}: {e}"))?,
        };
        let name = interval.name;
        let limit = max_position(min_shift, depth);
        if interval.beg > limit || interval.end > limit {
            return Err(match self.format {
                IndexFormat::Tabix => format!(
                    "line {number}: position {} on {:?} is beyond the tabix limit of 2^29 ({limit}); use a CSI index instead",
                    interval.end.max(interval.beg),
                    String::from_utf8_lossy(name),
                ),
                IndexFormat::Csi { .. } => format!(
                    "line {number}: position {} on {:?} is beyond the limit of {limit} for a CSI index with min_shift={min_shift} and depth={depth}; increase csi_depth",
                    interval.end.max(interval.beg),
                    String::from_utf8_lossy(name),
                ),
            });
        }
        let same_reference = self
            .last
            .is_some_and(|(tid, _, _)| self.names.get_index(tid).is_some_and(|n| n == name));
        let tid = if same_reference {
            let (tid, previous_beg, line) = self.last.expect("a previous record exists");
            if interval.beg < previous_beg {
                return Err(format!(
                    "line {number}: records are not sorted: this record on {:?} starts before the one on line {line}",
                    String::from_utf8_lossy(name),
                ));
            }
            tid
        } else if self.names.contains(name) {
            return Err(format!(
                "line {number}: records for {:?} are not contiguous; the lines of each reference must form one contiguous run",
                String::from_utf8_lossy(name),
            ));
        } else {
            self.names.insert_full(name.to_vec()).0
        };
        if interval.end < interval.beg {
            let start = interval.beg + i64::from(!columns.zero_based);
            return Err(format!(
                "line {number}: the end {} is before the start {start}",
                interval.end
            ));
        }
        self.bins = Some((min_shift, depth));
        self.line_number = number;
        self.last = Some((tid, interval.beg, number));
        Ok(Some((tid, interval.beg, interval.end)))
    }

    fn commit(&mut self, record: Option<(usize, i64, i64)>, end_position: LogicalPosition) {
        match record {
            Some((tid, beg, end)) => self.pending.push_back(Record {
                tid,
                beg,
                end,
                end_position,
            }),
            None if self.bins.is_none() => self.first_offset = end_position,
            None => {}
        }
    }

    fn resolve<W: Write>(&mut self, blocks: &mut BlockWriter<W>) {
        while let Some(record) = self.pending.front() {
            let Some(end_offset) = blocks.resolve(record.end_position) else {
                break;
            };
            if self.builder.is_none() {
                let (min_shift, depth) = self.bins.expect("bins are decided at the first record");
                let first = blocks
                    .resolve(self.first_offset)
                    .expect("the first offset precedes the first record");
                self.builder = Some(IndexBuilder::new(min_shift, depth, first));
            }
            let builder = self.builder.as_mut().expect("the builder exists");
            builder.push(record.tid, record.beg, record.end, end_offset);
            self.pending.pop_front();
        }
        let needed = match (&self.builder, self.pending.front()) {
            (None, _) => self.first_offset.block,
            (Some(_), Some(record)) => record.end_position.block,
            (Some(_), None) => blocks.position().block,
        };
        blocks.forget_before(needed);
    }

    fn write_index(&mut self, final_offset: u64, first_offset: u64) -> Result<()> {
        let columns = self
            .columns
            .take()
            .ok_or_else(|| Error::Invalid(NO_COLUMNS.into()))?;
        let (min_shift, depth) = self.bins.unwrap_or_else(|| self.empty_bins());
        let builder = self
            .builder
            .take()
            .unwrap_or_else(|| IndexBuilder::new(min_shift, depth, first_offset));
        let layout = Layout {
            format: index::format(columns.kind, columns.zero_based),
            refname: columns.refname,
            start: columns.start,
            end: columns.end,
            meta_char: columns.meta_char,
            skip_lines: columns.skip_lines as u32,
        };
        let file = self
            .file
            .take()
            .ok_or_else(|| io::Error::other("the index was already written"))?;
        let file = BufWriter::new(file);
        let mut file = match self.format {
            IndexFormat::Tabix => index::write_tabix(
                file,
                &index::tabix(builder, final_offset, &layout, &self.names),
            )?,
            IndexFormat::Csi { .. } => index::write_csi(
                file,
                &index::csi(builder, final_offset, &layout, &self.names),
            )?,
        };
        file.flush()?;
        Ok(())
    }
}

/// Checks a compression level, which must be between 0 and 12.
pub fn check_level(level: i64) -> Result<CompressionLevel> {
    u8::try_from(level)
        .ok()
        .and_then(|level| CompressionLevel::new(level).ok())
        .ok_or_else(|| Error::Invalid(format!("level must be between 0 and 12, not {level}")))
}

/// Checks the CSI bin parameters and returns the index format they describe.
pub fn csi_format(min_shift: i64, depth: Option<i64>) -> Result<IndexFormat> {
    let min_shift = u32::try_from(min_shift)
        .ok()
        .filter(|shift| (1..=TABIX_MAX_SHIFT).contains(shift))
        .ok_or_else(|| {
            Error::Invalid(format!(
                "csi_min_shift must be between 1 and {TABIX_MAX_SHIFT}, not {min_shift}"
            ))
        })?;
    let depth = depth
        .map(|depth| {
            u32::try_from(depth)
                .ok()
                .filter(|depth| (1..=9).contains(depth))
                .ok_or_else(|| {
                    Error::Invalid(format!("csi_depth must be between 1 and 9, not {depth}"))
                })
        })
        .transpose()?;
    Ok(IndexFormat::Csi { min_shift, depth })
}

/// Checks writer options before anything is created.
pub fn check_options(
    level: i64,
    threads: i64,
    index: Option<&IndexOptions>,
) -> Result<(CompressionLevel, NonZero<usize>)> {
    let compression = check_level(level)?;
    let threads = crate::check_threads(threads).map_err(Error::Invalid)?;
    if let Some(options) = index {
        if let Some(columns) = &options.columns {
            columns.validate().map_err(Error::Invalid)?;
        }
        if let IndexFormat::Csi { min_shift, depth } = options.format {
            csi_format(i64::from(min_shift), depth.map(i64::from))?;
        }
    }
    Ok((compression, threads))
}

/// A BGZF writer that optionally builds a tabix or CSI index from the lines it writes.
pub struct Writer<W: Write> {
    blocks: BlockWriter<W>,
    indexer: Option<Indexer>,
    final_columns: Option<Columns>,
    io_failure: Option<String>,
    finished: bool,
}

impl<W: Write> Writer<W> {
    /// Creates a writer compressing at `level` (0-12) on `threads` threads.
    pub fn new(sink: W, level: u8, threads: usize, index: Option<IndexOptions>) -> Result<Self> {
        let (level, threads) = check_options(i64::from(level), threads as i64, index.as_ref())?;
        Ok(Self {
            blocks: BlockWriter::new(sink, level, threads)?,
            indexer: index.map(Indexer::create).transpose()?,
            final_columns: None,
            io_failure: None,
            finished: false,
        })
    }

    /// Returns the columns in use, which may not be decided yet when they are being inferred.
    pub fn columns(&self) -> Option<&Columns> {
        match &self.indexer {
            Some(indexer) => indexer.columns.as_ref(),
            None => self.final_columns.as_ref(),
        }
    }

    /// Returns true if writing `len` more bytes may compress a block or wait on workers.
    pub fn may_block(&self, len: usize) -> bool {
        self.blocks.buffered() + len >= BLOCK_SIZE
    }

    fn check_stream(&self) -> Result<()> {
        if self.finished {
            return Err(Error::Io(io::Error::other(
                "I/O operation on a closed writer",
            )));
        }
        self.check_io()
    }

    fn check(&self) -> Result<()> {
        self.check_stream()?;
        if let Some(message) = self
            .indexer
            .as_ref()
            .and_then(|indexer| indexer.failure.as_ref())
        {
            return Err(Error::Invalid(format!(
                "indexing failed earlier: {message}"
            )));
        }
        Ok(())
    }

    fn guard<T>(&mut self, result: Result<T>) -> Result<T> {
        match &result {
            Err(Error::Io(error)) => self.io_failure = Some(error.to_string()),
            Err(Error::Invalid(message)) => {
                if let Some(indexer) = &mut self.indexer {
                    indexer.failure = Some(message.clone());
                }
            }
            Ok(_) => {}
        }
        result
    }

    /// Writes `data`, indexing every line it completes.
    pub fn write(&mut self, data: &[u8]) -> Result<()> {
        self.check()?;
        let result = self.write_unchecked(data);
        self.guard(result)
    }

    fn write_unchecked(&mut self, data: &[u8]) -> Result<()> {
        let Some(indexer) = &mut self.indexer else {
            self.blocks.write(data)?;
            self.blocks.forget_before(u64::MAX);
            return Ok(());
        };
        let mut rest = data;
        let mut block = self.blocks.position().block;
        while let Some(newline) = memchr(b'\n', rest) {
            let mut line = std::mem::take(&mut indexer.partial);
            let classified = if line.is_empty() {
                indexer.classify(trim_line(&rest[..newline]))
            } else {
                line.extend_from_slice(&rest[..newline]);
                indexer.classify(trim_line(&line))
            };
            line.clear();
            indexer.partial = line;
            let record = classified.map_err(Error::Invalid)?;
            self.blocks.write(&rest[..=newline])?;
            indexer.commit(record, self.blocks.position());
            rest = &rest[newline + 1..];
            if self.blocks.position().block != block {
                block = self.blocks.position().block;
                indexer.resolve(&mut self.blocks);
            }
        }
        self.blocks.write(rest)?;
        indexer.partial.extend_from_slice(rest);
        indexer.resolve(&mut self.blocks);
        Ok(())
    }

    /// Ends the current block and flushes everything written so far to the sink.
    pub fn flush(&mut self) -> Result<()> {
        self.check_stream()?;
        let result = self.blocks.flush().map_err(Error::Io);
        if let Some(indexer) = &mut self.indexer {
            indexer.resolve(&mut self.blocks);
        }
        self.guard(result)
    }

    /// Returns the virtual position of the next byte to be written.
    pub fn tell(&mut self) -> Result<u64> {
        self.check_stream()?;
        let result = self.blocks.tell().map_err(Error::Io);
        self.guard(result)
    }

    /// Returns true once the writer has been finished.
    pub fn is_finished(&self) -> bool {
        self.finished
    }

    /// Writes the end-of-file marker and then the index, if one was requested.
    ///
    /// After an indexing error, the data is written without the end-of-file marker, so that
    /// readers see it as truncated. After any error, the index file, which was created empty
    /// with the writer, is removed. An indexing error already reported by [`Writer::write`] is
    /// not reported again.
    pub fn finish(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        let Some(mut indexer) = self.indexer.take() else {
            self.check_io()?;
            return self.blocks.finish().map(drop).map_err(Error::Io);
        };
        let reported = indexer.failure.is_some();
        match self.finish_indexed(&mut indexer) {
            Ok(()) => Ok(()),
            Err(error) => {
                indexer.remove_file()?;
                match error {
                    Error::Invalid(_) if reported => Ok(()),
                    error => Err(error),
                }
            }
        }
    }

    fn finish_indexed(&mut self, indexer: &mut Indexer) -> Result<()> {
        self.check_io()?;
        if indexer.failure.is_none() && !indexer.partial.is_empty() {
            let line = std::mem::take(&mut indexer.partial);
            match indexer.classify(trim_line(&line)) {
                Ok(record) => indexer.commit(record, self.blocks.position()),
                Err(message) => indexer.failure = Some(message),
            }
        }
        self.final_columns.clone_from(&indexer.columns);
        if indexer.failure.is_none() && indexer.columns.is_none() {
            indexer.failure = Some(NO_COLUMNS.into());
        }
        if let Some(message) = indexer.failure.take() {
            self.blocks.abandon()?;
            return Err(Error::Invalid(message));
        }
        let final_offset = self.blocks.finish()?;
        indexer.resolve(&mut self.blocks);
        let first_offset = self
            .blocks
            .resolve(indexer.first_offset)
            .unwrap_or(final_offset);
        indexer.write_index(final_offset, first_offset)
    }

    fn check_io(&self) -> Result<()> {
        match &self.io_failure {
            Some(message) => Err(Error::Io(io::Error::other(format!(
                "the writer failed earlier: {message}"
            )))),
            None => Ok(()),
        }
    }

    /// Returns the sink.
    pub fn get_ref(&self) -> &W {
        self.blocks.get_ref()
    }

    /// Returns the sink mutably.
    pub fn get_mut(&mut self) -> &mut W {
        self.blocks.get_mut()
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;

    fn options(
        dir: &tempfile::TempDir,
        format: IndexFormat,
        columns: Option<Columns>,
    ) -> IndexOptions {
        IndexOptions {
            format,
            path: dir.path().join("out.bed.gz.idx"),
            columns,
            bed_only: false,
        }
    }

    fn decompress(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        bgzf::Reader::new(data).read_to_end(&mut out).unwrap();
        out
    }

    fn error(result: Result<()>) -> String {
        match result {
            Err(Error::Invalid(message)) => message,
            Err(Error::Io(error)) => panic!("unexpected I/O error: {error}"),
            Ok(()) => panic!("expected an error"),
        }
    }

    #[test]
    fn writes_without_an_index() {
        let mut writer = Writer::new(Vec::new(), 6, 1, None).unwrap();
        writer.write(b"hello\nworld").unwrap();
        writer.finish().unwrap();
        assert_eq!(decompress(writer.get_ref()), b"hello\nworld");
    }

    #[test]
    fn unsorted_records_name_the_line() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, Some(Columns::bed()))),
        )
        .unwrap();
        writer.write(b"#header\nchr1\t10\t20\n").unwrap();
        let message = error(writer.write(b"chr1\t5\t20\n"));
        assert!(message.starts_with("line 3:"), "{message}");
        assert!(error(writer.write(b"chr1\t50\t60\n")).contains("failed earlier"));
        writer.finish().unwrap();
        assert!(!dir.path().join("out.bed.gz.idx").exists());
    }

    #[test]
    fn references_must_be_contiguous() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, Some(Columns::bed()))),
        )
        .unwrap();
        writer.write(b"chr1\t1\t2\nchr2\t1\t2\n").unwrap();
        assert!(error(writer.write(b"chr1\t3\t4\n")).contains("not contiguous"));
    }

    #[test]
    fn tabix_positions_are_limited() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, Some(Columns::bed()))),
        )
        .unwrap();
        assert!(error(writer.write(b"chr1\t1\t600000000\n")).contains("CSI"));
        let csi = IndexFormat::Csi {
            min_shift: 14,
            depth: None,
        };
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, csi, Some(Columns::bed()))),
        )
        .unwrap();
        writer.write(b"chr1\t1\t600000000\n").unwrap();
        writer.finish().unwrap();
    }

    #[test]
    fn a_final_line_without_a_newline_is_checked_at_finish() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, Some(Columns::bed()))),
        )
        .unwrap();
        writer.write(b"chr1\t10\t20\nchr1\t5").unwrap();
        assert!(error(writer.finish()).starts_with("line 2:"));
        assert_eq!(decompress(writer.get_ref()), b"chr1\t10\t20\nchr1\t5");
    }

    #[test]
    fn inferred_columns_are_decided_while_streaming() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, None)),
        )
        .unwrap();
        writer.write(b"track name=x\n#c").unwrap();
        assert!(writer.columns().is_none());
        writer.write(b"omment\nchr1\t0\t5\n").unwrap();
        assert_eq!(
            writer.columns(),
            Some(&Columns {
                skip_lines: 1,
                ..Columns::bed()
            })
        );
        writer.finish().unwrap();
        assert!(dir.path().join("out.bed.gz.idx").exists());
    }

    #[test]
    fn inference_without_data_fails_at_finish() {
        let dir = tempfile::tempdir().unwrap();
        let mut writer = Writer::new(
            Vec::new(),
            6,
            1,
            Some(options(&dir, IndexFormat::Tabix, None)),
        )
        .unwrap();
        writer.write(b"# only a comment\n").unwrap();
        assert!(error(writer.finish()).contains("could not infer"));
    }

    #[test]
    fn csi_bins_follow_htslib() {
        assert_eq!(csi_bins(14, 0), Ok((14, 8)));
        assert_eq!(csi_bins(10, 0), Ok((10, 9)));
        assert_eq!(csi_bins(14, 248_956_422), Ok((14, 6)));
        assert_eq!(csi_bins(14, 1 << 33), Ok((14, 7)));
        assert_eq!(csi_bins(4, 1 << 40), Ok((14, 9)));
        assert_eq!(csi_bins(14, (1 << 62) - 256), Ok((35, 9)));
        assert!(csi_bins(14, (1 << 62) - 255).is_err());
        assert!(csi_bins(14, i64::MAX).is_err());
    }

    #[test]
    fn a_first_record_ending_blocks_after_the_header_is_indexed() {
        let dir = tempfile::tempdir().unwrap();
        let options = options(&dir, IndexFormat::Tabix, Some(Columns::bed()));
        let mut writer = Writer::new(Vec::new(), 1, 1, Some(options)).unwrap();
        writer.write(b"#header\n").unwrap();
        let mut line = b"chr1\t1\t2\t".to_vec();
        line.resize(3 * BLOCK_SIZE, b'x');
        line.push(b'\n');
        writer.write(&line).unwrap();
        writer.write(b"chr1\t5\t6\n").unwrap();
        writer.finish().unwrap();
        assert!(dir.path().join("out.bed.gz.idx").exists());
    }

    #[test]
    fn memory_stays_bounded_within_one_large_write() {
        let data: Vec<u8> = (0..200_000)
            .flat_map(|i| format!("chr1\t{i}\t{}\n", i + 5).into_bytes())
            .collect();
        for threads in [1, 3] {
            let dir = tempfile::tempdir().unwrap();
            let options = options(&dir, IndexFormat::Tabix, Some(Columns::bed()));
            let mut writer = Writer::new(Vec::new(), 1, threads, Some(options)).unwrap();
            writer.write(&data).unwrap();
            let pending = writer.indexer.as_ref().unwrap().pending.len();
            assert!(pending < 50_000, "threads={threads} pending={pending}");
            writer.finish().unwrap();
            let mut writer = Writer::new(Vec::new(), 1, threads, None).unwrap();
            writer.write(&data).unwrap();
            assert!(writer.blocks.remembered_blocks() <= 1, "threads={threads}");
        }
    }

    #[test]
    fn contig_lengths_are_read_from_headers() {
        assert_eq!(
            longest(
                b"##contig=<ID=1,length=249250621>",
                b"##contig",
                b"length",
                true
            ),
            Some(249_250_621)
        );
        assert_eq!(
            longest(b"@SQ\tSN:1\tLN:1000", b"@SQ", b"\tLN:", false),
            Some(1000)
        );
        assert_eq!(
            longest(b"##INFO=<ID=x>", b"##contig", b"length", true),
            None
        );
        assert_eq!(
            longest(b"@SQ\tLN:99999999999999999999", b"@SQ", b"\tLN:", false),
            Some(i64::MAX)
        );
    }
}
