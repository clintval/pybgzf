"""Threads sharing one reader, writer, or query take turns, as with the standard library's files."""

import bisect
import faulthandler
import functools
import gzip
import io
import random
import shutil
import threading
from collections import Counter
from collections.abc import Callable
from collections.abc import Iterator
from pathlib import Path

import pytest
from pybgzf import BgzfReader
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexedReader
from pybgzf import IndexFormat
from typing_extensions import override

from tests.helpers import HAS_HTSLIB
from tests.helpers import tabix
from tests.indexes import read_index

BLOCK_SIZE = 65280
CLOSED = "I/O operation on closed file."
EOF_MARKER = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
HEADER = b"#chrom\tstart\tend\tname\tpayload\n"
REFERENCES = ("chr1", "chr10", "chr2")
THREADS = 8
TIMEOUT = 60.0

RECORD = 12
RECORDS = 400_000


@pytest.fixture(autouse=True)
def watchdog() -> Iterator[None]:
    """Exit with every thread's traceback if a test hangs, even while a thread holds the GIL."""
    faulthandler.dump_traceback_later(2 * TIMEOUT, exit=True)
    yield
    faulthandler.cancel_dump_traceback_later()


class Worker(threading.Thread):
    """A daemon thread that keeps what its call returns or raises."""

    def __init__(self, call: Callable[[], object]) -> None:
        super().__init__(daemon=True)
        self.call: Callable[[], object] = call
        self.result: object = None
        self.error: BaseException | None = None
        self.start()

    @override
    def run(self) -> None:
        try:
            self.result = self.call()
        except BaseException as error:
            self.error = error

    def finish(self) -> object:
        """Wait for the call and return its result, or raise what it raised."""
        self.join(TIMEOUT)
        assert not self.is_alive(), "the thread did not finish"
        if self.error is not None:
            raise self.error
        return self.result


def concurrently(call: Callable[[int], None]) -> None:
    """Call `call` with each thread number, from that many threads started at once."""
    barrier = threading.Barrier(THREADS, timeout=TIMEOUT)

    def run(number: int) -> None:
        barrier.wait()
        call(number)

    workers = [Worker(functools.partial(run, number)) for number in range(THREADS)]
    for worker in workers:
        worker.finish()


def until_closed(step: Callable[[], bool]) -> None:
    """Call `step` until it returns False or raises because another thread closed the file."""
    try:
        while step():
            pass
    except ValueError as error:
        assert str(error) == CLOSED


class Gate:
    """A file-like object whose first read or write waits until the gate is opened."""

    def __init__(self, data: bytes = b"") -> None:
        self.buffer: io.BytesIO = io.BytesIO(data)
        self.entered: threading.Event = threading.Event()
        self.opened: threading.Event = threading.Event()

    def _wait(self) -> None:
        if not self.entered.is_set():
            self.entered.set()
            assert self.opened.wait(TIMEOUT)

    def read(self, size: int) -> bytes:
        self._wait()
        return self.buffer.read(size)

    def write(self, data: bytes) -> int:
        self._wait()
        return self.buffer.write(data)


def record(number: int) -> bytes:
    return f"{number:0{RECORD - 1}d}\n".encode()


@pytest.fixture(scope="module")
def records_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("records") / "records.gz"
    with BgzfWriter(path, threads=4) as writer:
        writer.write(b"".join(record(number) for number in range(RECORDS)))
    return path


def as_records(chunks: list[bytes]) -> list[bytes]:
    """Split whole records out of chunks, which must each hold only whole records."""
    for chunk in chunks:
        assert len(chunk) % RECORD == 0, chunk
    return [chunk[at : at + RECORD] for chunk in chunks for at in range(0, len(chunk), RECORD)]


def read_some(reader: BgzfReader, rng: random.Random, most: int) -> bytes:
    """Read a line or up to `most` records."""
    if rng.random() < 0.5:
        return reader.readline()
    return reader.read(RECORD * rng.randint(1, most))


class Features:
    """Sorted BED lines on every reference, spanning many blocks and query batches."""

    def __init__(self, path: Path) -> None:
        rng = random.Random(3)
        self.lines: dict[str, list[str]] = {}
        self.starts: dict[str, list[int]] = {}
        for name in REFERENCES:
            position = 0
            lines: list[str] = []
            for number in range(30_000):
                position += rng.randint(0, 300)
                end = position + rng.choice([1, 50, 5_000, 200_000])
                lines.append(f"{name}\t{position}\t{end}\tfeature{number}")
            self.lines[name] = lines
            self.starts[name] = [int(line.split("\t")[1]) for line in lines]
        self.path: Path = path
        with BgzfWriter(path, threads=2, index=IndexFormat.TBI, columns=Columns.BED) as writer:
            for name in REFERENCES:
                writer.write("".join(f"{line}\n" for line in self.lines[name]).encode())

    def overlapping(self, refname: str, start: int, end: int) -> list[str]:
        """The lines overlapping `[start, end)`, found without an index."""
        first = bisect.bisect_left(self.starts[refname], start - 200_000)
        last = bisect.bisect_left(self.starts[refname], end)
        selected = self.lines[refname][first:last]
        return [line for line in selected if int(line.split("\t")[2]) > start and start < end]


@pytest.fixture(scope="module")
def features(tmp_path_factory: pytest.TempPathFactory) -> Features:
    return Features(tmp_path_factory.mktemp("features") / "features.bed.gz")


def feature(name: str, start: int, label: str, rng: random.Random, payload: int) -> bytes:
    """A BED line with a random end and up to `payload` bytes after its name."""
    end = start + rng.choice([1, 40, 700, 20_000, 300_000, 5_000_000])
    return f"{name}\t{start}\t{end}\t{label}\t{'x' * rng.randint(0, payload)}\n".encode()


def assert_complete_and_indexed(path: Path, index: IndexFormat, expected: list[bytes]) -> None:
    """Assert the file holds each expected line once and whole, and its index finds them."""
    compressed = path.read_bytes()
    assert compressed.endswith(EOF_MARKER)
    lines = gzip.decompress(compressed).splitlines(keepends=True)
    assert lines[0] == HEADER
    assert Counter(lines[1:]) == Counter(expected)
    fields = [line.decode().removesuffix("\n").split("\t") for line in lines[1:]]
    with IndexedReader(path) as reader:
        for name in REFERENCES:
            for start, end in [(0, 10**10), (0, 1), (30_000, 60_000), (100_000, 10**7)]:
                assert list(reader.query(name, start, end)) == [
                    "\t".join(field)
                    for field in fields
                    if field[0] == name and int(field[1]) < end and int(field[2]) > start
                ]
    if HAS_HTSLIB:
        copy = path.parent / "htslib" / path.name
        copy.parent.mkdir()
        shutil.copyfile(path, copy)
        tabix("-f", *(["-C"] if index is IndexFormat.CSI else []), "-p", "bed", copy)
        suffix = f".{index.name.lower()}"
        ours, theirs = Path(f"{path}{suffix}"), Path(f"{copy}{suffix}")
        assert read_index(ours) == read_index(theirs)
        assert gzip.decompress(ours.read_bytes()) == gzip.decompress(theirs.read_bytes())


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
@pytest.mark.parametrize("threads", [1, 4])
def test_threads_share_one_writer(tmp_path: Path, index: IndexFormat, threads: int) -> None:
    path = tmp_path / "shared.bed.gz"
    rounds = 30
    turn = threading.Barrier(THREADS, timeout=TIMEOUT)
    written: list[list[bytes]] = [[] for _ in range(THREADS)]
    with BgzfWriter(path, threads=threads, index=index, columns=Columns.BED) as writer:
        writer.write(HEADER)

        def write(number: int) -> None:
            rng = random.Random(number)
            try:
                for round_number in range(rounds):
                    name = REFERENCES[round_number * len(REFERENCES) // rounds]
                    start = round_number % (rounds // len(REFERENCES)) * 50_000
                    lines = [
                        feature(name, start, f"thread{number}-{round_number}-{n}", rng, 2_000)
                        for n in range(25)
                    ]
                    while lines:
                        size = rng.randint(1, 3)
                        batch, lines = lines[:size], lines[size:]
                        assert writer.write(b"".join(batch)) == sum(map(len, batch))
                        written[number].extend(batch)
                    turn.wait()
            except BaseException:
                turn.abort()
                raise

        concurrently(write)
    assert_complete_and_indexed(path, index, [line for lines in written for line in lines])


@pytest.mark.parametrize("threads", [1, 3])
def test_threads_share_one_reader(records_path: Path, threads: int) -> None:
    chunks: list[list[bytes]] = [[] for _ in range(THREADS)]
    with BgzfReader(records_path, threads=threads) as reader:

        def read(number: int) -> None:
            rng = random.Random(number)
            while chunk := read_some(reader, rng, 400):
                chunks[number].append(chunk)

        concurrently(read)
    for got in chunks:
        assert as_records(got) == sorted(as_records(got))
    assert sorted(as_records([chunk for got in chunks for chunk in got])) == [
        record(number) for number in range(RECORDS)
    ]


@pytest.mark.parametrize("threads", [1, 2])
def test_threads_share_one_query(features: Features, threads: int) -> None:
    expected = features.overlapping("chr10", 0, 10**10)
    received: list[list[str]] = [[] for _ in range(THREADS)]
    with IndexedReader(features.path, threads=threads) as reader:
        lines = reader.query("chr10", 0, 10**10)

        def iterate(number: int) -> None:
            received[number].extend(lines)

        concurrently(iterate)
    assert Counter(line for got in received for line in got) == Counter(expected)
    order = {line: at for at, line in enumerate(expected)}
    for got in received:
        positions = [order[line] for line in got]
        assert positions == sorted(positions)


def test_threads_share_one_indexed_reader(features: Features) -> None:
    with IndexedReader(features.path, threads=2) as reader:

        def query(number: int) -> None:
            rng = random.Random(number)
            for _ in range(20):
                name = rng.choice(REFERENCES)
                start = rng.randint(0, 9_000_000)
                end = start + rng.choice([0, 1, 1_000, 100_000, 3_000_000])
                assert list(reader.query(name, start, end)) == features.overlapping(
                    name, start, end
                )
                assert reader.refnames == list(REFERENCES)
                assert reader.columns == Columns.BED

        concurrently(query)


def test_a_waiting_thread_releases_the_interpreter(records_path: Path) -> None:
    gate = Gate(records_path.read_bytes())
    reader = BgzfReader(gate)
    first = Worker(lambda: reader.read(RECORD))
    try:
        assert gate.entered.wait(TIMEOUT)
        second = Worker(lambda: reader.read(RECORD))
        second.join(0.2)
        assert second.is_alive()
    finally:
        gate.opened.set()
    assert first.finish() == record(0)
    assert second.finish() == record(1)
    reader.close()


def test_closed_does_not_wait_for_a_read(records_path: Path) -> None:
    gate = Gate(records_path.read_bytes())
    reader = BgzfReader(gate)
    inner = getattr(reader, "_inner")  # noqa: B009
    reading = Worker(lambda: inner.read(RECORD))
    try:
        assert gate.entered.wait(TIMEOUT)
        assert Worker(lambda: inner.closed).finish() is False
    finally:
        gate.opened.set()
    assert reading.finish() == record(0)
    reader.close()
    assert inner.closed


def test_close_waits_for_a_read(records_path: Path) -> None:
    gate = Gate(records_path.read_bytes())
    reader = BgzfReader(gate)
    reading = Worker(lambda: reader.read(RECORD))
    try:
        assert gate.entered.wait(TIMEOUT)
        closing = Worker(reader.close)
        closing.join(0.2)
        assert closing.is_alive()
    finally:
        gate.opened.set()
    assert reading.finish() == record(0)
    closing.finish()
    assert reader.closed
    with pytest.raises(ValueError, match=CLOSED):
        reader.read(1)


def test_close_waits_for_a_write() -> None:
    gate = Gate()
    data = bytes(range(256)) * 1_000
    writer = BgzfWriter(gate)
    writing = Worker(lambda: writer.write(data))
    try:
        assert gate.entered.wait(TIMEOUT)
        closing = Worker(writer.close)
        closing.join(0.2)
        assert closing.is_alive()
    finally:
        gate.opened.set()
    assert writing.finish() == len(data)
    closing.finish()
    assert gzip.decompress(gate.buffer.getvalue()) == data
    assert gate.buffer.getvalue().endswith(EOF_MARKER)


@pytest.mark.parametrize("threads", [1, 3])
def test_close_races_reads(records_path: Path, threads: int) -> None:
    chunks: list[list[bytes]] = [[] for _ in range(THREADS)]
    started = threading.Event()
    reader = BgzfReader(records_path, threads=threads)

    def read(number: int) -> None:
        if number == 0:
            assert started.wait(TIMEOUT)
            reader.close()
            return
        rng = random.Random(number)

        def step() -> bool:
            chunk = read_some(reader, rng, 40)
            chunks[number].append(chunk)
            if len(chunks[number]) == 100:
                started.set()
            return bool(chunk)

        until_closed(step)

    concurrently(read)
    assert reader.closed
    got = sorted(as_records([chunk for got in chunks for chunk in got]))
    assert got == [record(number) for number in range(len(got))]


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_close_races_writes(tmp_path: Path, index: IndexFormat) -> None:
    path = tmp_path / "raced.bed.gz"
    written: list[list[bytes]] = [[] for _ in range(THREADS)]
    started = threading.Event()
    writer = BgzfWriter(path, threads=2, index=index, columns=Columns.BED)
    writer.write(HEADER)

    def write(number: int) -> None:
        if number == 0:
            assert started.wait(TIMEOUT)
            writer.close()
            return
        rng = random.Random(number)

        def step() -> bool:
            line = feature("chr1", 1_000, f"thread{number}-{len(written[number])}", rng, 200)
            writer.write(line)
            written[number].append(line)
            if len(written[number]) % 50 == 0:
                writer.flush()
                assert writer.tell() > 0
            if len(written[number]) == 200:
                started.set()
            return len(written[number]) < 20_000

        until_closed(step)

    concurrently(write)
    assert writer.closed
    assert_complete_and_indexed(path, index, [line for lines in written for line in lines])


def test_close_races_queries(features: Features) -> None:
    received: list[list[str]] = [[] for _ in range(THREADS)]
    started = threading.Event()
    reader = IndexedReader(features.path, threads=2)
    regions = [
        (REFERENCES[number % len(REFERENCES)], number * 1_000, 10**10) for number in range(8)
    ]

    def query(number: int) -> None:
        if number == 0:
            assert started.wait(TIMEOUT)
            reader.close()
            return

        def step() -> bool:
            received[number].clear()
            for line in reader.query(*regions[number]):
                received[number].append(line)
                if len(received[number]) == 1_000:
                    started.set()
            return True

        until_closed(step)

    concurrently(query)
    assert reader.closed
    for number, got in enumerate(received[1:], start=1):
        assert got == features.overlapping(*regions[number])[: len(got)]


def close_one(files: list[io.RawIOBase], number: int) -> None:
    files[number % len(files)].close()


def test_threads_may_close_at_once(records_path: Path) -> None:
    for _ in range(200):
        files: list[io.RawIOBase] = [BgzfWriter(io.BytesIO()), BgzfReader(records_path)]
        concurrently(functools.partial(close_one, files))
        assert all(file.closed for file in files)


class Reentrant:
    """A source and sink that reads from or writes to the reader or writer calling it."""

    def __init__(self) -> None:
        self.reader: BgzfReader | None = None
        self.writer: BgzfWriter | None = None

    def read(self, size: int) -> bytes:
        assert self.reader is not None
        return self.reader.read(size)

    def write(self, data: bytes) -> int:
        assert self.writer is not None
        return self.writer.write(data)


def test_reentrant_calls_raise_instead_of_waiting_forever() -> None:
    sink = Reentrant()
    sink.writer = BgzfWriter(sink)
    with pytest.raises(RuntimeError, match="reentrant call"):
        sink.writer.write(b"x" * (2 * BLOCK_SIZE))
    with pytest.raises(OSError, match="failed earlier"):
        sink.writer.close()
    source = Reentrant()
    source.reader = BgzfReader(source)
    with pytest.raises(RuntimeError, match="reentrant call"):
        source.reader.read(10)
    source.reader.close()
