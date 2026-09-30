import gzip
import io
import os
import random
import shutil
import struct
import threading
import warnings
from pathlib import Path
from typing import cast

import pybgzf
import pytest
from pybgzf import BgzfReader
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexedReader
from pybgzf import IndexFormat
from pybgzf import ReadableBinary

from tests.helpers import bed_lines
from tests.helpers import requires_htslib
from tests.helpers import tabix

REFERENCES = ("chr1", "chr10", "chr2")


@pytest.fixture(scope="module")
def lines() -> list[str]:
    return bed_lines(seed=5, references=REFERENCES)


@pytest.fixture(scope="module")
def data_dir(tmp_path_factory: pytest.TempPathFactory, lines: list[str]) -> Path:
    directory = tmp_path_factory.mktemp("data")
    for index in (IndexFormat.TBI, IndexFormat.CSI):
        path = directory / index.name.lower() / "features.bed.gz"
        path.parent.mkdir()
        with pybgzf.writer(path, threads=3, index=index, columns=Columns.BED) as handle:
            handle.write("#chrom\tstart\tend\tname\tscore\n")
            handle.writelines(lines)
    return directory


def bed_path(data_dir: Path, index: IndexFormat) -> Path:
    return data_dir / index.name.lower() / "features.bed.gz"


def test_round_trips_with_every_thread_count(data_dir: Path, lines: list[str]) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    expected = gzip.decompress(path.read_bytes())
    assert expected.decode().endswith("".join(lines))
    for threads in (1, 2, 5):
        with BgzfReader(path, threads=threads) as reader:
            assert reader.readall() == expected
        with pybgzf.reader(path, threads=threads) as handle:
            assert handle.read() == expected.decode()


def test_matches_stdlib_gzip_line_by_line(data_dir: Path) -> None:
    path = bed_path(data_dir, IndexFormat.CSI)
    with gzip.open(path, "rt") as expected, pybgzf.reader(path, threads=4) as actual:
        assert list(actual) == list(expected)


def test_small_reads_and_readinto(data_dir: Path) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    expected = gzip.decompress(path.read_bytes())
    chunks: list[bytes] = []
    with BgzfReader(path, threads=2) as reader:
        buffer = bytearray(7777)
        while n := reader.readinto(buffer):
            chunks.append(bytes(buffer[:n]))
        assert reader.read(10) == b""
    assert b"".join(chunks) == expected


@pytest.mark.parametrize("threads", [1, 3])
def test_readline_and_iteration(data_dir: Path, lines: list[str], threads: int) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    with BgzfReader(path, threads=threads) as reader:
        assert reader.readline() == b"#chrom\tstart\tend\tname\tscore\n"
        assert reader.readline(0) == b""
        assert reader.readline(4) == b"chr1"
        assert reader.readline() == lines[0].encode()[4:]
        assert next(iter(reader)) == lines[1].encode()
        pieces = [lines[0].encode(), lines[1].encode()]
        while piece := reader.readline(1000):
            assert len(piece) <= 1000 and b"\n" not in piece[:-1]
            pieces.append(piece)
    assert b"".join(pieces) == "".join(lines).encode()


def test_tell_and_seek_use_virtual_offsets(data_dir: Path) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    for threads in (1, 3):
        with BgzfReader(path, threads=threads) as reader:
            assert not reader.seekable()
            marks: list[tuple[int, bytes]] = []
            while True:
                offset = reader.tell()
                line = reader.readline()
                if not line:
                    break
                marks.append((offset, line))
            rng = random.Random(threads)
            for offset, line in rng.sample(marks, 200):
                assert reader.seek(offset) == offset
                assert reader.tell() == offset
                assert reader.readline() == line
            with pytest.raises(io.UnsupportedOperation):
                reader.seek(0, io.SEEK_END)
            with pytest.raises(ValueError, match="negative"):
                reader.seek(-1)


def test_offsets_agree_with_the_writer(tmp_path: Path) -> None:
    path = tmp_path / "out.gz"
    offsets: list[int] = []
    with BgzfWriter(path, threads=2) as writer:
        for number in range(50_000):
            offsets.append(writer.tell())
            writer.write(f"line {number}\n".encode())
    with BgzfReader(path, threads=2) as reader:
        for number, offset in enumerate(offsets[::997]):
            reader.seek(offset)
            assert reader.readline() == f"line {number * 997}\n".encode()


def test_reads_file_like_objects_and_pipes(data_dir: Path) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    compressed = path.read_bytes()
    expected = gzip.decompress(compressed)
    with BgzfReader(io.BytesIO(compressed), threads=3) as reader:
        assert reader.readall() == expected
    read_fd, write_fd = os.pipe()

    def feed() -> None:
        with os.fdopen(write_fd, "wb") as pipe:
            pipe.write(compressed)

    writer = threading.Thread(target=feed)
    writer.start()
    with os.fdopen(read_fd, "rb") as pipe, BgzfReader(pipe, threads=2) as reader:
        if not pipe.seekable():
            with pytest.raises(OSError, match="not seekable"):
                reader.seek(0)
        assert reader.readall() == expected
    writer.join()


def test_closed_readers_reject_io(data_dir: Path) -> None:
    reader = BgzfReader(bed_path(data_dir, IndexFormat.TBI))
    assert reader.readable()
    reader.close()
    reader.close()
    for call in (reader.readall, reader.readline, reader.tell, lambda: reader.seek(0)):
        with pytest.raises(ValueError, match="closed"):
            call()


def test_invalid_reader_options(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="threads"):
        BgzfReader(tmp_path / "x.gz", threads=0)
    with pytest.raises(ValueError, match="threads must be between 1 and 1024"):
        BgzfReader(tmp_path / "x.gz", threads=100_000)
    with pytest.raises(ValueError, match="threads must be between 1 and 1024, not -1"):
        BgzfReader(tmp_path / "x.gz", threads=-1)
    with pytest.raises(FileNotFoundError):
        BgzfReader(tmp_path / "missing.gz")


def overlapping(lines: list[str], refname: str, start: int, end: int) -> list[str]:
    selected: list[str] = []
    for line in lines:
        name, first, last = line.split("\t")[:3]
        if name == refname and int(first) < end and int(last) > start and start < end:
            selected.append(line.removesuffix("\n"))
    return selected


def random_regions(seed: int, count: int) -> list[tuple[str, int, int]]:
    rng = random.Random(seed)
    regions = [
        ("chr1", 0, 0),
        ("chr1", 0, 1),
        ("chr2", 5_000, 5_000),
        ("chr10", 10**9, 10**9 + 10),
        ("chr10", 0, 10**12),
        ("chrUn", 0, 100),
    ]
    for _ in range(count):
        start = rng.randint(0, 25_000_000)
        width = rng.choice([0, 1, 50, 2_000, 100_000, 5_000_000])
        regions.append((rng.choice(REFERENCES), start, start + width))
    return regions


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
@pytest.mark.parametrize("threads", [1, 3])
def test_queries_match_brute_force(
    data_dir: Path, lines: list[str], index: IndexFormat, threads: int
) -> None:
    with IndexedReader(bed_path(data_dir, index), threads=threads) as reader:
        assert reader.refnames == list(REFERENCES)
        assert reader.columns == Columns.BED
        for refname, start, end in random_regions(seed=threads, count=150):
            assert list(reader.query(refname, start, end)) == overlapping(
                lines, refname, start, end
            )


def test_interleaved_queries(data_dir: Path, lines: list[str]) -> None:
    with IndexedReader(bed_path(data_dir, IndexFormat.TBI), threads=2) as reader:
        first = reader.query("chr1", 0, 10**9)
        second = reader.query("chr2", 0, 10**9)
        firsts: list[str] = []
        seconds: list[str] = []
        while True:
            a, b = next(first, None), next(second, None)
            if a is None and b is None:
                break
            firsts.extend([] if a is None else [a])
            seconds.extend([] if b is None else [b])
    assert firsts == overlapping(lines, "chr1", 0, 10**9)
    assert seconds == overlapping(lines, "chr2", 0, 10**9)


def test_query_arguments_are_checked(data_dir: Path) -> None:
    with IndexedReader(bed_path(data_dir, IndexFormat.TBI)) as reader:
        with pytest.raises(ValueError, match="start must be at least 0"):
            reader.query("chr1", -1, 5)
        with pytest.raises(ValueError, match="end at least start"):
            reader.query("chr1", 5, 4)
    assert reader.closed
    with pytest.raises(ValueError, match="closed"):
        reader.query("chr1", 0, 1)
    with pytest.raises(ValueError, match="threads must be between 1 and 1024"):
        IndexedReader(bed_path(data_dir, IndexFormat.TBI), threads=100_000)


def test_finds_the_index_next_to_the_file(tmp_path: Path, data_dir: Path) -> None:
    path = tmp_path / "features.bed.gz"
    shutil.copyfile(bed_path(data_dir, IndexFormat.TBI), path)
    with pytest.raises(FileNotFoundError, match="pass index_path"):
        IndexedReader(path)
    shutil.copyfile(f"{bed_path(data_dir, IndexFormat.TBI)}.tbi", f"{path}.tbi")
    with IndexedReader(path) as reader:
        assert reader.refnames == list(REFERENCES)
    Path(f"{path}.csi").write_bytes(b"not an index")
    with pytest.raises(ValueError, match="not a tabix or CSI index"):
        IndexedReader(path)
    with IndexedReader(path, index_path=f"{path}.tbi") as reader:
        assert list(reader.query("chr1", 0, 1)) == list(reader.query("chr1", 0, 1))


def test_queries_follow_the_columns_in_the_index(tmp_path: Path) -> None:
    path = tmp_path / "genes.gff.gz"
    text = (
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=a\n"
        "chr1\tsrc\tgene\t150\t400\t.\t+\t.\tID=b\n"
        "#comment\n"
        "chr1\tsrc\tgene\t401\t500\t.\t+\t.\tID=c\n"
    )
    with pybgzf.writer(path, index=IndexFormat.CSI, columns=pybgzf.INFER) as handle:
        handle.write(text)
    with IndexedReader(path) as reader:
        assert reader.columns == Columns.GFF
        assert [line.split("\t")[8] for line in reader.query("chr1", 199, 200)] == ["ID=a", "ID=b"]
        assert [line.split("\t")[8] for line in reader.query("chr1", 400, 401)] == ["ID=c"]
        assert list(reader.query("chr1", 0, 99)) == []


@requires_htslib
@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_queries_match_tabix(data_dir: Path, index: IndexFormat) -> None:
    path = bed_path(data_dir, index)
    with IndexedReader(path, threads=2) as reader:
        for refname, start, end in random_regions(seed=11, count=60):
            if start == end or start >= 2**29:
                continue
            expected = tabix(path, f"{refname}:{start + 1}-{end}").splitlines()
            assert list(reader.query(refname, start, end)) == expected


@requires_htslib
def test_reads_files_written_by_bgzip(tmp_path: Path, lines: list[str]) -> None:
    path = tmp_path / "bgzip.bed"
    path.write_text("".join(lines))
    import subprocess

    subprocess.run(["bgzip", "-@", "2", str(path)], check=True)
    subprocess.run(["tabix", "-p", "bed", f"{path}.gz"], check=True)
    with pybgzf.reader(f"{path}.gz", threads=2) as handle:
        assert handle.read() == "".join(lines)
    with IndexedReader(f"{path}.gz") as reader:
        assert list(reader.query("chr2", 1_000, 90_000)) == overlapping(
            lines, "chr2", 1_000, 90_000
        )


def corrupted(data_dir: Path, kind: str) -> bytes:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    if kind == "truncated":
        return compressed[: len(compressed) // 2]
    if kind == "garbage":
        return random.Random(3).randbytes(len(compressed))
    return gzip.compress(gzip.decompress(compressed))


CORRUPTIONS = ["truncated", "garbage", "plain gzip"]


@pytest.mark.parametrize("kind", CORRUPTIONS)
@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("from_path", [True, False], ids=["path", "file-like"])
def test_corrupt_input_raises(
    tmp_path: Path, data_dir: Path, kind: str, threads: int, from_path: bool
) -> None:
    data = corrupted(data_dir, kind)
    path = tmp_path / "corrupt.gz"
    path.write_bytes(data)

    def source() -> Path | io.BytesIO:
        return path if from_path else io.BytesIO(data)

    with pytest.raises(OSError), BgzfReader(source(), threads=threads) as reader:
        reader.readall()
    with pytest.raises(OSError), BgzfReader(source(), threads=threads) as reader:
        while reader.readline():
            pass
    with pytest.raises(OSError), pybgzf.reader(source(), threads=threads) as handle:
        handle.read()


@pytest.mark.filterwarnings("ignore::pybgzf.TruncatedWarning")
@pytest.mark.parametrize("kind", CORRUPTIONS)
@pytest.mark.parametrize("threads", [1, 4])
def test_corrupt_input_raises_from_queries(
    tmp_path: Path, data_dir: Path, kind: str, threads: int
) -> None:
    path = tmp_path / "corrupt.bed.gz"
    path.write_bytes(corrupted(data_dir, kind))
    shutil.copyfile(f"{bed_path(data_dir, IndexFormat.TBI)}.tbi", f"{path}.tbi")
    with IndexedReader(path, threads=threads) as reader, pytest.raises(OSError):
        for refname in REFERENCES:
            list(reader.query(refname, 0, 10**9))


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_queries_reach_the_last_position_an_index_holds(tmp_path: Path, index: IndexFormat) -> None:
    path = tmp_path / "edge.bed.gz"
    last = "chr1\t536870911\t536870912"
    with BgzfWriter(path, index=index, columns=Columns.BED, csi_depth=5) as writer:
        writer.write(f"chr1\t0\t1\n{last}\n".encode())
    with IndexedReader(path) as reader:
        assert list(reader.query("chr1", 536870911, 536870912)) == [last]
        assert list(reader.query("chr1", 536870900, 10**12)) == [last]
        assert list(reader.query("chr1", 536870912, 10**12)) == []


@pytest.mark.parametrize("csi_min_shift", [1, 2, 3])
def test_vcf_with_the_smallest_csi_bins_can_be_queried(tmp_path: Path, csi_min_shift: int) -> None:
    path = tmp_path / "small.vcf.gz"
    lines = [f"1\t{position}\t.\tA\tG\t.\t.\t." for position in range(1, 5_000, 7)]
    header = "##contig=<ID=1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    with pybgzf.writer(
        path, index=IndexFormat.CSI, columns=Columns.VCF, csi_min_shift=csi_min_shift
    ) as handle:
        handle.write(header + "".join(f"{line}\n" for line in lines))
    with IndexedReader(path) as reader:
        assert list(reader.query("1", 0, 10_000)) == lines
        assert list(reader.query("1", 98, 99)) == ["1\t99\t.\tA\tG\t.\t.\t."]


def test_csi_indexes_deeper_than_nine_levels_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "deep.vcf.gz"
    with pybgzf.writer(path) as handle:
        handle.write("1\t1\t.\tA\tG\t.\t.\t.\n")
    header = struct.pack("<7i", 2, 1, 2, 0, ord("#"), 0, 2) + b"1\0"
    index = struct.pack("<4s3i", b"CSI\x01", 2, 10, len(header)) + header + struct.pack("<2i", 1, 0)
    with BgzfWriter(f"{path}.csi") as writer:
        writer.write(index)
    with pytest.raises(ValueError, match="has 10 bin levels"):
        IndexedReader(path)


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="named pipes are not available")
@pytest.mark.parametrize("from_path", [True, False], ids=["path", "file-like"])
def test_close_returns_while_a_pipe_writer_is_idle(
    tmp_path: Path, data_dir: Path, from_path: bool
) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()[:40_000]
    fifo = tmp_path / "fifo"
    os.mkfifo(fifo)
    ready = threading.Event()

    def feed() -> None:
        with fifo.open("wb") as pipe:
            pipe.write(compressed)
            pipe.flush()
            ready.wait()

    feeder = threading.Thread(target=feed)
    feeder.start()
    with fifo.open("rb") as pipe:
        reader = BgzfReader(fifo if from_path else pipe, threads=2)
        assert len(reader.read(10)) == 10
        closer = threading.Thread(target=reader.close)
        closer.start()
        closer.join(timeout=5)
        closed = not closer.is_alive()
        ready.set()
        closer.join()
    feeder.join()
    assert closed


def test_lines_that_are_not_utf8_name_their_reference_and_offset(tmp_path: Path) -> None:
    path = tmp_path / "latin1.bed.gz"
    bad = "chr2\t10\t20\tcafé\n".encode("latin-1")
    with BgzfWriter(path, index=IndexFormat.TBI, columns=Columns.BED) as writer:
        writer.write(b"chr1\t1\t2\tok\nchr2\t5\t6\tok\n" + bad)
    with BgzfReader(path) as reader:
        offset = reader.tell()
        while reader.readline() != bad:
            offset = reader.tell()
    with IndexedReader(path) as reader:
        assert list(reader.query("chr1", 0, 10)) == ["chr1\t1\t2\tok"]
        lines = reader.query("chr2", 0, 100)
        assert next(lines) == "chr2\t5\t6\tok"
        with pytest.raises(ValueError, match=f'"chr2" at virtual offset {offset} is not UTF-8'):
            next(lines)


def test_reference_names_that_are_not_utf8_raise(tmp_path: Path) -> None:
    path = tmp_path / "latin1.bed.gz"
    with BgzfWriter(path, index=IndexFormat.TBI, columns=Columns.BED) as writer:
        writer.write("chré\t1\t2\n".encode("latin-1"))
    with IndexedReader(path) as reader, pytest.raises(ValueError, match="not UTF-8"):
        _ = reader.refnames


@pytest.mark.parametrize("threads", [1, 3])
def test_reads_files_without_an_end_of_file_marker(data_dir: Path, threads: int) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    expected = gzip.decompress(compressed)
    chunks: list[bytes] = []
    with (
        pytest.warns(pybgzf.TruncatedWarning, match="without an end-of-file marker") as caught,
        BgzfReader(io.BytesIO(compressed[:-28]), threads=threads) as reader,
    ):
        for _ in range(len(expected) // 65536 + 2):
            chunks.append(reader.read(65536))
    assert b"".join(chunks) == expected
    assert chunks[-1] == b""
    assert len(caught) == 1


@pytest.mark.parametrize("threads", [1, 3])
def test_complete_files_do_not_warn(data_dir: Path, threads: int) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with BgzfReader(path, threads=threads) as reader:
            end = len(reader.readall())
            assert reader.read(1) == b""
            assert end > 0
        with pybgzf.reader(path, threads=threads) as handle:
            assert handle.read()
        with BgzfReader(io.BytesIO(b""), threads=threads) as reader:
            assert reader.readall() == b""


@pytest.mark.parametrize("threads", [1, 3])
def test_a_path_without_an_end_of_file_marker_warns_naming_it(
    data_dir: Path, threads: int, tmp_path: Path
) -> None:
    path = tmp_path / "cut.bed.gz"
    _ = path.write_bytes(bed_path(data_dir, IndexFormat.TBI).read_bytes()[:-28])
    with (
        pytest.warns(pybgzf.TruncatedWarning, match="cut.bed.gz"),
        pybgzf.reader(path, threads=threads) as handle,
    ):
        assert handle.read()


@pytest.mark.parametrize("threads", [1, 3])
def test_a_pipe_without_an_end_of_file_marker_warns(data_dir: Path, threads: int) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    read_end, write_end = os.pipe()

    def write() -> None:
        with os.fdopen(write_end, "wb") as sink:
            _ = sink.write(compressed[:-28])

    writer = threading.Thread(target=write, daemon=True)
    writer.start()
    with (
        pytest.warns(pybgzf.TruncatedWarning),
        os.fdopen(read_end, "rb") as source,
        BgzfReader(source, threads=threads) as reader,
    ):
        assert reader.readall() == gzip.decompress(compressed)
    writer.join(30)
    assert not writer.is_alive()


def test_a_missing_end_of_file_marker_raises_when_warnings_are_errors(data_dir: Path) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    with warnings.catch_warnings():
        warnings.simplefilter("error", pybgzf.TruncatedWarning)
        with (
            BgzfReader(io.BytesIO(compressed[:-28])) as reader,
            pytest.raises(pybgzf.TruncatedWarning),
        ):
            _ = reader.readall()


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_an_indexed_file_without_an_end_of_file_marker_warns_when_opened(
    data_dir: Path, index: IndexFormat, tmp_path: Path
) -> None:
    source = bed_path(data_dir, index)
    suffix = ".tbi" if index is IndexFormat.TBI else ".csi"
    path = tmp_path / "cut.bed.gz"
    _ = path.write_bytes(source.read_bytes()[:-28])
    _ = shutil.copyfile(f"{source}{suffix}", f"{path}{suffix}")
    with pytest.warns(pybgzf.TruncatedWarning, match="cut.bed.gz"):
        reader = IndexedReader(path)
    with reader:
        assert list(reader.query("chr1", 0, 1 << 20))


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_a_query_of_a_file_its_index_does_not_describe_says_so(
    data_dir: Path, index: IndexFormat, tmp_path: Path
) -> None:
    suffix = ".tbi" if index is IndexFormat.TBI else ".csi"
    path = tmp_path / "other.bed.gz"
    with pybgzf.writer(path, columns=Columns.BED) as handle:
        _ = handle.write("chr1\t0\t1\tonly\n")
    with (
        IndexedReader(path, index_path=f"{bed_path(data_dir, index)}{suffix}") as reader,
        pytest.raises(OSError, match="changed since it was opened"),
    ):
        _ = list(reader.query("chr2", 0, 1 << 29))


def block_starts(data: bytes) -> list[int]:
    starts: list[int] = []
    at = 0
    while at < len(data):
        starts.append(at)
        at += int.from_bytes(data[at + 16 : at + 18], "little") + 1
    return starts


@pytest.mark.parametrize("threads", [1, 3])
@pytest.mark.parametrize("cut", [5, 20], ids=["in a header", "in a body"])
def test_truncation_inside_a_block_raises(data_dir: Path, threads: int, cut: int) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    start = block_starts(compressed)[3]
    with (
        pytest.raises(OSError, match="truncated|unexpected end of file|failed to fill"),
        BgzfReader(io.BytesIO(compressed[: start + cut]), threads=threads) as reader,
    ):
        reader.readall()


@pytest.mark.parametrize("threads", [1, 3])
def test_seeking_to_the_end(data_dir: Path, threads: int) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    size = path.stat().st_size
    with BgzfReader(path, threads=threads) as reader:
        first = reader.readline()
        rest = reader.readall()
        assert reader.read(10) == b""
        end = reader.tell()
        for position in (end, (size - 28) << 16, end):
            reader.seek(0)
            reader.readline()
            assert reader.seek(position) == position
            assert reader.tell() == position
            assert reader.readline() == b""
        reader.seek(0)
        assert reader.readline() == first
        assert reader.readall() == rest


@pytest.mark.parametrize("threads", [1, 3])
def test_seeking_outside_the_file_raises(data_dir: Path, threads: int) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    size = path.stat().st_size
    with BgzfReader(path, threads=threads) as reader:
        first = reader.readline()
        for position in (65535, (size + 100) << 16, 2**64 - 1, 2**64, -1):
            with pytest.raises(ValueError, match="not in this file|negative|too large"):
                reader.seek(position)
        reader.seek(0)
        assert reader.readline() == first


class _Source:
    def __init__(self, data: bytes, kind: str) -> None:
        self.data: io.BytesIO = io.BytesIO(data)
        self.kind: str = kind

    def read(self, size: int, /) -> object:
        if self.kind == "raises" and self.data.tell() > 100_000:
            raise PermissionError("the source failed")
        if self.kind == "str":
            return "text"
        if self.kind == "none":
            return None
        if self.kind == "too much":
            return self.data.read(size + 1)
        if self.kind == "short":
            return self.data.read(min(size, 7))
        return self.data.read(size)


def readable(data: bytes, kind: str) -> ReadableBinary:
    return cast(ReadableBinary, cast(object, _Source(data, kind)))


@pytest.mark.parametrize("threads", [1, 3])
@pytest.mark.parametrize(
    ("kind", "error"),
    [("raises", PermissionError), ("str", TypeError), ("none", OSError), ("too much", OSError)],
)
def test_source_errors_raise(
    data_dir: Path, threads: int, kind: str, error: type[Exception]
) -> None:
    compressed = bed_path(data_dir, IndexFormat.TBI).read_bytes()
    with pytest.raises(error), BgzfReader(readable(compressed, kind), threads=threads) as reader:
        reader.readall()


@pytest.mark.parametrize("threads", [1, 3])
def test_short_reads_from_a_source(threads: int) -> None:
    stream = io.BytesIO()
    with BgzfWriter(stream) as writer:
        writer.write(b"".join(b"line %d\n" % number for number in range(20_000)))
    compressed = stream.getvalue()
    with BgzfReader(readable(compressed, "short"), threads=threads) as reader:
        assert reader.readall() == gzip.decompress(compressed)


def test_sources_without_read_raise_at_once() -> None:
    with pytest.raises(TypeError, match="read"):
        BgzfReader(cast(str, cast(object, b"features.bed.gz")))


def test_query_bounds_beyond_64_bits(data_dir: Path, lines: list[str]) -> None:
    with IndexedReader(bed_path(data_dir, IndexFormat.TBI)) as reader:
        assert list(reader.query("chr2", 0, 10**30)) == overlapping(lines, "chr2", 0, 10**12)
        assert list(reader.query("chr2", 10**30, 10**31)) == []


def test_open_reader_checks_the_encoding_first(data_dir: Path) -> None:
    with pytest.raises(LookupError):
        pybgzf.reader(bed_path(data_dir, IndexFormat.TBI), threads=2, encoding="no-such-codec")


@pytest.mark.filterwarnings("ignore::pybgzf.TruncatedWarning")
@pytest.mark.parametrize("threads", [1, 3])
def test_errors_name_the_file(tmp_path: Path, data_dir: Path, threads: int) -> None:
    missing = tmp_path / "missing.gz"
    with pytest.raises(FileNotFoundError) as error:
        BgzfReader(missing, threads=threads)
    assert error.value.filename == str(missing)
    path = tmp_path / "truncated.gz"
    path.write_bytes(corrupted(data_dir, "truncated"))
    with pytest.raises(OSError, match="truncated.gz"), BgzfReader(path, threads=threads) as reader:
        reader.readall()
    shutil.copyfile(f"{bed_path(data_dir, IndexFormat.TBI)}.tbi", f"{path}.tbi")
    with (
        pytest.raises(OSError, match="truncated.gz"),
        IndexedReader(path, threads=threads) as indexed,
    ):
        for refname in REFERENCES:
            list(indexed.query(refname, 0, 10**9))
