//! A BGZF block writer that knows the virtual position of every byte it has written.
//!
//! Blocks are compressed with the public [`bgzf::Compressor`], either in the calling thread or on
//! a pool of worker threads, and are always written in order. Each block gets a number when it is
//! handed off for compression, so a byte's position is first known as a block number and an
//! offset ([`LogicalPosition`]) and is resolved to a virtual position once its block has been
//! written. This design is adapted from the multithreaded writer in fgumi
//! (<https://github.com/fulcrumgenomics/fgumi>, MIT licensed, Copyright Fulcrum Genomics LLC).

use std::collections::VecDeque;
use std::io::{self, Write};
use std::num::NonZero;
use std::thread::{self, JoinHandle};

use bgzf::{BgzfError, CompressionLevel, Compressor};
use crossbeam_channel::{Receiver, Sender, TryRecvError, bounded};

/// The number of uncompressed bytes in a full block, as in `bgzip`; flushing ends a block early.
pub const BLOCK_SIZE: usize = bgzf::BGZF_BLOCK_SIZE;

/// The position of a byte as a block number and an offset into that block's uncompressed data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct LogicalPosition {
    pub block: u64,
    pub offset: u32,
}

type Compressed = io::Result<(Vec<u8>, Vec<u8>)>;
type Job = (Vec<u8>, Vec<u8>, Sender<Compressed>);

enum Engine {
    Serial(Compressor),
    Parallel {
        jobs: Option<Sender<Job>>,
        pending: VecDeque<Receiver<Compressed>>,
        workers: Vec<JoinHandle<()>>,
        max_pending: usize,
    },
}

/// Writes BGZF blocks to a sink while tracking where each block starts.
pub struct BlockWriter<W: Write> {
    sink: W,
    buf: Vec<u8>,
    engine: Engine,
    next_block: u64,
    blocks_written: u64,
    bytes_written: u64,
    window: VecDeque<(u64, u32)>,
    window_base: u64,
    spare_inputs: Vec<Vec<u8>>,
    spare_outputs: Vec<Vec<u8>>,
    finished: bool,
}

fn to_io(error: BgzfError) -> io::Error {
    match error {
        BgzfError::Io(e) => e,
        other => io::Error::other(other.to_string()),
    }
}

fn compress(compressor: &mut Compressor, input: &[u8], output: &mut Vec<u8>) -> io::Result<()> {
    output.clear();
    compressor.compress(input, output).map_err(to_io)
}

fn spawn_workers(
    level: CompressionLevel,
    count: usize,
    jobs: &Receiver<Job>,
) -> io::Result<Vec<JoinHandle<()>>> {
    (0..count)
        .map(|_| {
            let jobs = jobs.clone();
            thread::Builder::new().spawn(move || {
                let mut compressor = Compressor::new(level);
                while let Ok((input, mut output, reply)) = jobs.recv() {
                    let result =
                        compress(&mut compressor, &input, &mut output).map(|()| (input, output));
                    let _ = reply.send(result);
                }
            })
        })
        .collect()
}

impl<W: Write> BlockWriter<W> {
    /// Creates a writer that compresses in the calling thread when `threads` is one, and on
    /// `threads` worker threads otherwise.
    pub fn new(sink: W, level: CompressionLevel, threads: NonZero<usize>) -> io::Result<Self> {
        let engine = if threads.get() == 1 {
            Engine::Serial(Compressor::new(level))
        } else {
            let (sender, receiver) = bounded(threads.get() * 2);
            let workers = spawn_workers(level, threads.get(), &receiver)?;
            Engine::Parallel {
                jobs: Some(sender),
                pending: VecDeque::new(),
                workers,
                max_pending: threads.get() * 2,
            }
        };
        Ok(Self {
            sink,
            buf: Vec::with_capacity(BLOCK_SIZE),
            engine,
            next_block: 0,
            blocks_written: 0,
            bytes_written: 0,
            window: VecDeque::new(),
            window_base: 0,
            spare_inputs: Vec::new(),
            spare_outputs: Vec::new(),
            finished: false,
        })
    }

    /// Returns the position of the next byte to be written.
    pub fn position(&self) -> LogicalPosition {
        LogicalPosition {
            block: self.next_block,
            offset: self.buf.len() as u32,
        }
    }

    /// Returns the number of uncompressed bytes waiting in the current block.
    pub fn buffered(&self) -> usize {
        self.buf.len()
    }

    /// Buffers `data`, compressing each block as soon as it is full.
    pub fn write(&mut self, mut data: &[u8]) -> io::Result<()> {
        self.check_open()?;
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.buf.len()).min(data.len());
            self.buf.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.buf.len() == BLOCK_SIZE {
                self.send()?;
            }
        }
        Ok(())
    }

    /// Ends the current block, if it holds any data, and hands it off for compression.
    pub fn end_block(&mut self) -> io::Result<()> {
        self.check_open()?;
        if self.buf.is_empty() {
            Ok(())
        } else {
            self.send()
        }
    }

    /// Ends the current block, waits for every block to be written, and flushes the sink.
    pub fn flush(&mut self) -> io::Result<()> {
        self.end_block()?;
        self.drain(0)?;
        self.sink.flush()
    }

    /// Returns the virtual position of the next byte, waiting for earlier blocks to be written.
    pub fn tell(&mut self) -> io::Result<u64> {
        self.check_open()?;
        self.drain(0)?;
        Ok((self.bytes_written << 16) | self.buf.len() as u64)
    }

    /// Writes all remaining data and, if `eof`, the BGZF end-of-file marker, returning the
    /// virtual position of the marker. Without it, readers see the stream as truncated.
    pub fn finish(&mut self, eof: bool) -> io::Result<u64> {
        self.end_block()?;
        self.drain(0)?;
        self.stop_workers();
        let eof_start = self.bytes_written;
        if eof {
            let mut marker = Vec::new();
            Compressor::append_eof(&mut marker);
            self.sink.write_all(&marker)?;
        }
        self.sink.flush()?;
        self.finished = true;
        Ok(eof_start << 16)
    }

    #[cfg(test)]
    pub fn get_ref(&self) -> &W {
        &self.sink
    }

    /// Returns the sink mutably.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.sink
    }

    /// Resolves a position to a virtual position, if its block has been written.
    ///
    /// A position at the very end of a block resolves to the start of the next block, which is
    /// where a BGZF reader such as htslib reports it to be.
    pub fn resolve(&self, position: LogicalPosition) -> Option<u64> {
        let LogicalPosition { block, offset } = position;
        if block >= self.blocks_written {
            let at_end = self.finished && block == self.blocks_written && offset == 0;
            return at_end.then_some(self.bytes_written << 16);
        }
        let index = usize::try_from(block.checked_sub(self.window_base)?).ok()?;
        let (start, size) = *self.window.get(index)?;
        if offset < size {
            return Some((start << 16) | u64::from(offset));
        }
        let next = if block + 1 < self.blocks_written {
            self.window.get(index + 1)?.0
        } else {
            self.bytes_written
        };
        Some(next << 16)
    }

    /// Forgets where blocks before `block` start, keeping at least the last written block.
    pub fn forget_before(&mut self, block: u64) {
        while self.window_base < block && self.window.len() > 1 {
            self.window.pop_front();
            self.window_base += 1;
        }
    }

    #[cfg(test)]
    pub(crate) fn remembered_blocks(&self) -> usize {
        self.window.len()
    }

    fn check_open(&self) -> io::Result<()> {
        if self.finished {
            return Err(io::Error::other("the BGZF stream is already finished"));
        }
        Ok(())
    }

    fn send(&mut self) -> io::Result<()> {
        let spare = self
            .spare_inputs
            .pop()
            .unwrap_or_else(|| Vec::with_capacity(BLOCK_SIZE));
        let input = std::mem::replace(&mut self.buf, spare);
        let mut output = self.spare_outputs.pop().unwrap_or_default();
        self.next_block += 1;
        match &mut self.engine {
            Engine::Serial(compressor) => {
                compress(compressor, &input, &mut output)?;
                self.emit(input, output)
            }
            Engine::Parallel {
                jobs,
                pending,
                max_pending,
                ..
            } => {
                let (reply, receiver) = bounded(1);
                if jobs
                    .as_ref()
                    .is_none_or(|jobs| jobs.send((input, output, reply)).is_err())
                {
                    return Err(io::Error::other("compression workers have stopped"));
                }
                pending.push_back(receiver);
                let max_pending = *max_pending;
                self.drain(max_pending)
            }
        }
    }

    fn emit(&mut self, mut input: Vec<u8>, output: Vec<u8>) -> io::Result<()> {
        self.sink.write_all(&output)?;
        self.window
            .push_back((self.bytes_written, input.len() as u32));
        self.bytes_written += output.len() as u64;
        self.blocks_written += 1;
        input.clear();
        self.spare_inputs.push(input);
        self.spare_outputs.push(output);
        Ok(())
    }

    /// Writes compressed blocks in order: those already compressed and, while more than
    /// `max_pending` are pending, the next one once it is compressed.
    fn drain(&mut self, max_pending: usize) -> io::Result<()> {
        let exited = || io::Error::other("a compression worker exited unexpectedly");
        loop {
            let Engine::Parallel { pending, .. } = &mut self.engine else {
                return Ok(());
            };
            let Some(receiver) = pending.front() else {
                return Ok(());
            };
            let result = if pending.len() > max_pending {
                receiver.recv().map_err(|_| exited())?
            } else {
                match receiver.try_recv() {
                    Ok(result) => result,
                    Err(TryRecvError::Empty) => return Ok(()),
                    Err(TryRecvError::Disconnected) => return Err(exited()),
                }
            };
            pending.pop_front();
            let (input, output) = result?;
            self.emit(input, output)?;
        }
    }

    fn stop_workers(&mut self) {
        if let Engine::Parallel {
            jobs,
            workers,
            pending,
            ..
        } = &mut self.engine
        {
            pending.clear();
            jobs.take();
            for worker in workers.drain(..) {
                let _ = worker.join();
            }
        }
    }
}

impl<W: Write> Drop for BlockWriter<W> {
    fn drop(&mut self) {
        self.stop_workers();
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;

    fn level(n: u8) -> CompressionLevel {
        CompressionLevel::new(n).unwrap()
    }

    fn threads(n: usize) -> NonZero<usize> {
        NonZero::new(n).unwrap()
    }

    fn decompress(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        bgzf::Reader::new(data).read_to_end(&mut out).unwrap();
        out
    }

    fn block_starts(data: &[u8]) -> Vec<u64> {
        let mut starts = Vec::new();
        let mut at = 0usize;
        while at < data.len() {
            starts.push(at as u64);
            let bsize = u16::from_le_bytes([data[at + 16], data[at + 17]]) as usize;
            at += bsize + 1;
        }
        starts
    }

    fn sample(len: usize) -> Vec<u8> {
        (0..len)
            .map(|i| b"ACGT\tchr1\n"[(i * 7 + i / 13) % 10])
            .collect()
    }

    #[test]
    fn empty_stream_is_only_the_eof_marker() {
        let mut writer = BlockWriter::new(Vec::new(), level(6), threads(1)).unwrap();
        assert_eq!(writer.finish(true).unwrap(), 0);
        assert_eq!(writer.get_ref().len(), 28);
        assert!(decompress(writer.get_ref()).is_empty());
    }

    #[test]
    fn round_trips_and_is_identical_across_thread_counts() {
        let data = sample(BLOCK_SIZE * 7 + 1234);
        let mut outputs = Vec::new();
        for n in [1, 2, 3, 8] {
            let mut writer = BlockWriter::new(Vec::new(), level(6), threads(n)).unwrap();
            for chunk in data.chunks(10_007) {
                writer.write(chunk).unwrap();
            }
            writer.finish(true).unwrap();
            assert_eq!(decompress(writer.get_ref()), data);
            outputs.push(writer.get_ref().clone());
        }
        assert!(outputs.windows(2).all(|pair| pair[0] == pair[1]));
    }

    #[test]
    fn positions_resolve_to_block_starts() {
        for n in [1, 4] {
            let mut writer = BlockWriter::new(Vec::new(), level(6), threads(n)).unwrap();
            let mut marks = Vec::new();
            for chunk in sample(BLOCK_SIZE * 3 + 100).chunks(5000) {
                writer.write(chunk).unwrap();
                marks.push(writer.position());
            }
            writer.write(b"x").unwrap();
            writer.finish(true).unwrap();
            let starts = block_starts(writer.get_ref());
            for mark in marks {
                let expected = (starts[mark.block as usize] << 16) | u64::from(mark.offset);
                assert_eq!(
                    writer.resolve(mark),
                    Some(expected),
                    "threads={n} mark={mark:?}"
                );
            }
        }
    }

    #[test]
    fn a_full_block_moves_the_position_to_the_next_block() {
        let mut writer = BlockWriter::new(Vec::new(), level(6), threads(1)).unwrap();
        writer.write(&sample(BLOCK_SIZE)).unwrap();
        assert_eq!(
            writer.position(),
            LogicalPosition {
                block: 1,
                offset: 0
            }
        );
    }

    #[test]
    fn the_end_of_a_flushed_block_resolves_to_the_next_block() {
        for n in [1, 3] {
            let mut writer = BlockWriter::new(Vec::new(), level(6), threads(n)).unwrap();
            writer.write(b"chr1\t1\t2\n").unwrap();
            let end_of_first = writer.position();
            writer.flush().unwrap();
            writer.write(b"chr1\t2\t3\n").unwrap();
            let end_of_second = writer.position();
            let eof = writer.finish(true).unwrap();
            let starts = block_starts(writer.get_ref());
            assert_eq!(writer.resolve(end_of_first), Some(starts[1] << 16));
            assert_eq!(writer.resolve(end_of_second), Some(eof));
            assert_eq!(eof, starts[2] << 16);
        }
    }

    #[test]
    fn unwritten_positions_do_not_resolve() {
        let mut writer = BlockWriter::new(Vec::new(), level(6), threads(1)).unwrap();
        writer.write(b"abc").unwrap();
        assert_eq!(writer.resolve(writer.position()), None);
        writer.end_block().unwrap();
        assert_eq!(
            writer.resolve(LogicalPosition {
                block: 0,
                offset: 1
            }),
            Some(1)
        );
    }

    #[test]
    fn tell_reports_the_virtual_position_of_the_next_byte() {
        for n in [1, 2] {
            let mut writer = BlockWriter::new(Vec::new(), level(6), threads(n)).unwrap();
            assert_eq!(writer.tell().unwrap(), 0);
            writer.write(&sample(BLOCK_SIZE + 10)).unwrap();
            let told = writer.tell().unwrap();
            writer.finish(true).unwrap();
            let starts = block_starts(writer.get_ref());
            assert_eq!(told, (starts[1] << 16) | 10);
        }
    }

    #[test]
    fn forgetting_old_blocks_keeps_recent_positions() {
        let mut writer = BlockWriter::new(Vec::new(), level(1), threads(1)).unwrap();
        writer.write(&sample(BLOCK_SIZE * 4)).unwrap();
        writer.forget_before(3);
        let last = LogicalPosition {
            block: 3,
            offset: 5,
        };
        let starts = block_starts(writer.get_ref());
        assert_eq!(writer.resolve(last), Some((starts[3] << 16) | 5));
    }

    #[test]
    fn level_zero_stores_blocks() {
        let data = sample(BLOCK_SIZE * 2);
        let mut writer = BlockWriter::new(Vec::new(), level(0), threads(2)).unwrap();
        writer.write(&data).unwrap();
        writer.finish(true).unwrap();
        assert_eq!(decompress(writer.get_ref()), data);
    }

    struct Broken;

    impl Write for Broken {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn sink_errors_surface() {
        for n in [1, 2] {
            let mut writer = BlockWriter::new(Broken, level(6), threads(n)).unwrap();
            let result = writer
                .write(&sample(BLOCK_SIZE * 8))
                .and_then(|()| writer.finish(true).map(drop));
            assert_eq!(result.unwrap_err().kind(), io::ErrorKind::BrokenPipe);
        }
    }
}
