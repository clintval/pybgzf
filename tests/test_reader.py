import gzip
import io
import os
import random
import shutil
import threading
from pathlib import Path

import pybgzf
import pytest
from pybgzf import BgzfReader
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexedReader
from pybgzf import IndexFormat

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
        path = directory / index.value / "features.bed.gz"
        path.parent.mkdir()
        with pybgzf.open(path, threads=3, index=index, columns=Columns.BED) as handle:
            handle.write("#chrom\tstart\tend\tname\tscore\n")
            handle.writelines(lines)
    return directory


def bed_path(data_dir: Path, index: IndexFormat) -> Path:
    return data_dir / index.value / "features.bed.gz"


def test_round_trips_with_every_thread_count(data_dir: Path, lines: list[str]) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    expected = gzip.decompress(path.read_bytes())
    assert expected.decode().endswith("".join(lines))
    for threads in (1, 2, 5):
        with BgzfReader(path, threads=threads) as reader:
            assert reader.readall() == expected
        with pybgzf.open_reader(path, threads=threads) as handle:
            assert handle.read() == expected.decode()


def test_matches_stdlib_gzip_line_by_line(data_dir: Path) -> None:
    path = bed_path(data_dir, IndexFormat.CSI)
    with gzip.open(path, "rt") as expected, pybgzf.open_reader(path, threads=4) as actual:
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


def test_readline_and_iteration(data_dir: Path, lines: list[str]) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    with BgzfReader(path) as reader:
        assert reader.readline() == b"#chrom\tstart\tend\tname\tscore\n"
        assert reader.readline(4) == b"chr1"
        assert reader.readline() == lines[0].encode()[4:]
        assert next(iter(reader)) == lines[1].encode()


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
    with pybgzf.open(path, index=IndexFormat.CSI, columns=pybgzf.INFER) as handle:
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
    with pybgzf.open_reader(f"{path}.gz", threads=2) as handle:
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
    with pytest.raises(OSError), pybgzf.open_reader(source(), threads=threads) as handle:
        handle.read()


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


@pytest.mark.parametrize("threads", [1, 4])
def test_seeking_after_the_end_reads_again(data_dir: Path, threads: int) -> None:
    path = bed_path(data_dir, IndexFormat.TBI)
    with BgzfReader(path, threads=threads) as reader:
        first = reader.readline()
        end = len(reader.readall()) + len(first)
        assert reader.read(10) == b""
        assert reader.tell() == reader.tell()
        reader.seek(0)
        assert reader.readline() == first
        assert len(first) + len(reader.readall()) == end


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
