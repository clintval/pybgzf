import csv
import gzip
import io
import os
import struct
import threading
from collections.abc import Callable
from pathlib import Path

import pybgzf
import pytest
from pybgzf import INFER
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexFormat

from tests.helpers import bed_text
from tests.indexes import ParsedIndex
from tests.indexes import read_index

BLOCK_SIZE = 65280


def block_starts(data: bytes) -> list[int]:
    """Return the compressed offset of every BGZF block."""
    starts: list[int] = []
    at = 0
    while at < len(data):
        starts.append(at)
        at += struct.unpack_from("<H", data, at + 16)[0] + 1
    return starts


def write_bytes(
    path: Path,
    data: bytes,
    chunk: int = 10_000,
    *,
    level: int = 6,
    threads: int = 1,
    index: IndexFormat | None = None,
    columns: Columns | None = None,
) -> None:
    with BgzfWriter(path, level=level, threads=threads, index=index, columns=columns) as writer:
        for start in range(0, len(data), chunk):
            writer.write(data[start : start + chunk])


def test_round_trips_through_gzip(tmp_path: Path) -> None:
    data = bed_text().encode()
    path = tmp_path / "out.bed.gz"
    write_bytes(path, data)
    assert gzip.decompress(path.read_bytes()) == data


def test_accepts_bytearray_and_memoryview(tmp_path: Path) -> None:
    path = tmp_path / "out.gz"
    with BgzfWriter(path) as writer:
        assert writer.write(bytearray(b"ab")) == 2
        assert writer.write(memoryview(b"cde")) == 3
        assert writer.write(memoryview(b"xxfgxx")[2:4]) == 2
    assert gzip.decompress(path.read_bytes()) == b"abcdefg"


@pytest.mark.parametrize("level", [0, 1, 6, 12])
def test_levels(tmp_path: Path, level: int) -> None:
    data = bed_text().encode()[:300_000]
    path = tmp_path / "out.gz"
    write_bytes(path, data, level=level)
    assert gzip.decompress(path.read_bytes()) == data


def test_output_is_identical_across_thread_counts(tmp_path: Path) -> None:
    data = bed_text().encode()
    outputs: list[tuple[bytes, ParsedIndex]] = []
    for threads in (1, 2, 4, 7):
        path = tmp_path / f"out{threads}.bed.gz"
        write_bytes(path, data, threads=threads, index=IndexFormat.TBI, columns=Columns.BED)
        outputs.append((path.read_bytes(), read_index(Path(f"{path}.tbi"))))
    assert all(output == outputs[0] for output in outputs)


def test_empty_file_is_only_the_eof_marker(tmp_path: Path) -> None:
    path = tmp_path / "empty.gz"
    BgzfWriter(path).close()
    assert len(path.read_bytes()) == 28
    assert gzip.decompress(path.read_bytes()) == b""


def test_tell_returns_virtual_offsets(tmp_path: Path) -> None:
    for threads in (1, 3):
        path = tmp_path / f"out{threads}.gz"
        told: list[int] = []
        with BgzfWriter(path, threads=threads) as writer:
            assert writer.tell() == 0
            for _ in range(5):
                writer.write(b"x" * 30_000)
                told.append(writer.tell())
        starts = block_starts(path.read_bytes())
        expected: list[int] = []
        for written in range(30_000, 150_001, 30_000):
            block, offset = divmod(written, BLOCK_SIZE)
            expected.append((starts[block] << 16) | offset)
        assert told == expected


def test_flush_ends_the_block(tmp_path: Path) -> None:
    path = tmp_path / "out.gz"
    with BgzfWriter(path) as writer:
        writer.write(b"abc")
        writer.flush()
        assert writer.tell() == len(path.read_bytes()) << 16
        writer.write(b"def")
    assert len(block_starts(path.read_bytes())) == 3


def test_closed_writer_rejects_io(tmp_path: Path) -> None:
    writer = BgzfWriter(tmp_path / "out.gz")
    assert writer.writable() and not writer.readable() and not writer.seekable()
    writer.close()
    writer.close()
    assert writer.closed
    for call in (lambda: writer.write(b"x"), writer.flush, writer.tell):
        with pytest.raises(ValueError, match="closed"):
            call()


CSI = IndexFormat.CSI
BED = Columns.BED


INVALID_OPTIONS: list[tuple[Callable[[Path], BgzfWriter], str]] = [
    (lambda path: BgzfWriter(path, level=13), "level"),
    (lambda path: BgzfWriter(path, threads=0), "threads"),
    (lambda path: BgzfWriter(path, threads=100_000), "threads must be between 1 and 1024"),
    (
        lambda path: BgzfWriter(path, index=IndexFormat.TBI),
        "columns is required when index is set",
    ),
    (lambda path: BgzfWriter(path, index_path="x.tbi"), "only used when index is set"),
    (lambda path: BgzfWriter(path, index=CSI, columns=BED, csi_min_shift=0), "csi_min_shift"),
    (lambda path: BgzfWriter(path, index=CSI, columns=BED, csi_depth=12), "csi_depth"),
]


@pytest.mark.parametrize(("make", "message"), INVALID_OPTIONS)
def test_invalid_options_create_nothing(
    tmp_path: Path, make: Callable[[Path], BgzfWriter], message: str
) -> None:
    path = tmp_path / "out.bed.gz"
    with pytest.raises(ValueError, match=message):
        make(path)
    assert not path.exists()


def test_missing_directory_is_an_os_error(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        BgzfWriter(tmp_path / "missing" / "out.gz")


def test_columns_without_an_index_are_unused(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "out.gz", columns=Columns.BED) as writer:
        writer.write(b"not\tsorted\tat all\n")
        assert writer.columns == Columns.BED
    assert not (tmp_path / "out.gz.tbi").exists()


def test_writes_to_a_file_like_object_with_an_index_path(tmp_path: Path) -> None:
    buffer = io.BytesIO()
    index_path = tmp_path / "stream.tbi"
    data = bed_text().encode()
    with BgzfWriter(
        buffer, index=IndexFormat.TBI, index_path=index_path, columns=Columns.BED
    ) as writer:
        writer.write(data)
    assert not buffer.closed
    assert gzip.decompress(buffer.getvalue()) == data
    path = tmp_path / "file.bed.gz"
    write_bytes(path, data, index=IndexFormat.TBI, columns=Columns.BED)
    assert buffer.getvalue() == path.read_bytes()
    assert read_index(index_path) == read_index(Path(f"{path}.tbi"))


def test_file_like_objects_need_an_index_path() -> None:
    with pytest.raises(ValueError, match="a file-like object, which is not a regular file"):
        BgzfWriter(io.BytesIO(), index=IndexFormat.TBI, columns=Columns.BED)


def test_writes_to_a_pipe(tmp_path: Path) -> None:
    read_fd, write_fd = os.pipe()
    received = bytearray()

    def drain() -> None:
        with os.fdopen(read_fd, "rb") as pipe:
            received.extend(pipe.read())

    reader = threading.Thread(target=drain)
    reader.start()
    data = bed_text().encode()
    with (
        os.fdopen(write_fd, "wb") as pipe,
        BgzfWriter(
            pipe, threads=2, index=IndexFormat.CSI, index_path=tmp_path / "p.csi", columns=INFER
        ) as writer,
    ):
        writer.write(data)
    reader.join()
    assert gzip.decompress(bytes(received)) == data
    assert writer.columns == Columns.BED
    assert (tmp_path / "p.csi").exists()


class _Broken:
    def write(self, data: bytes) -> int:
        raise BrokenPipeError(f"reader went away before {len(data)} bytes")


def test_sink_errors_keep_their_type() -> None:
    writer = BgzfWriter(_Broken())
    with pytest.raises(BrokenPipeError, match="reader went away"):
        writer.write(b"x" * (BLOCK_SIZE + 1))
    with pytest.raises(OSError, match="failed earlier"):
        writer.write(b"x")


def test_text_mode(tmp_path: Path) -> None:
    path = tmp_path / "out.bed.gz"
    with pybgzf.open(path, index=IndexFormat.TBI, columns=Columns.BED) as handle:
        assert isinstance(handle, io.TextIOWrapper)
        handle.write("chr1\t1\t10\tcafé\n")
        handle.write("chr1\t5\t10\tnaïve\n")
    assert gzip.decompress(path.read_bytes()).decode() == "chr1\t1\t10\tcafé\nchr1\t5\t10\tnaïve\n"
    assert Path(f"{path}.tbi").exists()


def test_csv_writer(tmp_path: Path) -> None:
    path = tmp_path / "out.bed.gz"
    with pybgzf.open(path, newline="", index=IndexFormat.TBI, columns=Columns.BED) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#chrom", "start", "end"])
        for start in range(0, 100_000, 10):
            writer.writerow(["chr1", start, start + 5])
    text = gzip.decompress(path.read_bytes()).decode()
    assert text.startswith("#chrom\tstart\tend\r\nchr1\t0\t5\r\n")
    assert len(text.splitlines()) == 10_001


def test_unsorted_records_raise_at_write_naming_the_line(tmp_path: Path) -> None:
    path = tmp_path / "out.bed.gz"
    with BgzfWriter(path, index=IndexFormat.TBI, columns=Columns.BED) as writer:
        writer.write(b"#header\nchr1\t100\t200\n")
        with pytest.raises(ValueError, match=r"^line 3: records are not sorted"):
            writer.write(b"chr1\t50\t60\n")
        with pytest.raises(ValueError, match="indexing failed earlier"):
            writer.write(b"chr1\t300\t400\n")
    assert gzip.decompress(path.read_bytes()) == b"#header\nchr1\t100\t200\n"
    assert not Path(f"{path}.tbi").exists()


def test_a_failed_write_removes_a_stale_index(tmp_path: Path) -> None:
    path = tmp_path / "out.bed.gz"
    Path(f"{path}.tbi").write_bytes(b"stale")
    with pytest.raises(ValueError, match="not contiguous"):
        with BgzfWriter(path, index=IndexFormat.TBI, columns=Columns.BED) as writer:
            writer.write(b"chr1\t1\t2\nchr2\t1\t2\nchr1\t3\t4\n")
    assert not Path(f"{path}.tbi").exists()


def test_bad_lines_are_value_errors(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "a.bed.gz", index=IndexFormat.TBI, columns=Columns.BED) as writer:
        with pytest.raises(ValueError, match=r"line 1: column 2 is not an integer"):
            writer.write(b"chr1\tone\t2\n")
    with BgzfWriter(tmp_path / "b.bed.gz", index=IndexFormat.TBI, columns=Columns.BED) as writer:
        with pytest.raises(ValueError, match=r"line 1: the end 5 is before the start 11"):
            writer.write(b"chr1\t10\t5\n")


def test_tabix_limit_suggests_csi(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "a.bed.gz", index=IndexFormat.TBI, columns=Columns.BED) as writer:
        with pytest.raises(ValueError, match=r"2\^29.*CSI"):
            writer.write(b"chr1\t1\t600000000\n")
    with BgzfWriter(tmp_path / "b.bed.gz", index=IndexFormat.CSI, columns=Columns.BED) as writer:
        writer.write(b"chr1\t1\t600000000\n")
    with BgzfWriter(
        tmp_path / "c.bed.gz", index=IndexFormat.CSI, columns=Columns.BED, csi_depth=5
    ) as writer:
        with pytest.raises(ValueError, match="increase csi_depth"):
            writer.write(b"chr1\t1\t600000000\n")


def test_an_unterminated_last_line_is_checked_at_close(tmp_path: Path) -> None:
    writer = BgzfWriter(tmp_path / "out.bed.gz", index=IndexFormat.TBI, columns=Columns.BED)
    writer.write(b"chr1\t10\t20\nchr1\t5\t9")
    with pytest.raises(ValueError, match=r"^line 2"):
        writer.close()
    assert writer.closed


def test_lines_split_across_writes_index_like_whole_writes(tmp_path: Path) -> None:
    data = bed_text().encode()[:200_000]
    whole = tmp_path / "whole.bed.gz"
    write_bytes(whole, data, chunk=len(data), index=IndexFormat.TBI, columns=Columns.BED)
    pieces = tmp_path / "pieces.bed.gz"
    with BgzfWriter(pieces, index=IndexFormat.TBI, columns=Columns.BED) as writer:
        for start in range(0, len(data), 3):
            writer.write(data[start : start + 3])
    assert pieces.read_bytes() == whole.read_bytes()
    assert read_index(Path(f"{pieces}.tbi")) == read_index(Path(f"{whole}.tbi"))


def test_crlf_lines_index_like_lf_lines(tmp_path: Path) -> None:
    lf = b"chr1\t1\t10\nchr1\t20\t30\nchr2\t5\t6\n"
    crlf = lf.replace(b"\n", b"\r\n")
    for name, data in (("lf", lf), ("crlf", crlf)):
        write_bytes(tmp_path / f"{name}.bed.gz", data, index=IndexFormat.TBI, columns=Columns.BED)
    ours, theirs = (read_index(tmp_path / f"{name}.bed.gz.tbi") for name in ("lf", "crlf"))
    assert ours.header == theirs.header
    assert [len(ref.bins) for ref in ours.references] == [
        len(ref.bins) for ref in theirs.references
    ]


def test_columns_are_inferred_from_the_path(tmp_path: Path) -> None:
    path = tmp_path / "calls.vcf.gz"
    with BgzfWriter(path, index=IndexFormat.TBI, columns=INFER) as writer:
        assert writer.columns == Columns.VCF


def test_columns_are_inferred_from_the_index_path(tmp_path: Path) -> None:
    with BgzfWriter(
        io.BytesIO(), index=IndexFormat.TBI, index_path=tmp_path / "x.gff3.gz.tbi", columns=INFER
    ) as w:
        assert w.columns == Columns.GFF


def test_suffix_beats_content(tmp_path: Path) -> None:
    path = tmp_path / "reads.sam.gz"
    with BgzfWriter(path, index=IndexFormat.TBI, columns=INFER) as writer:
        assert writer.columns == Columns.SAM
        with pytest.raises(ValueError, match="line 1: expected at least 4"):
            writer.write(b"chr1\t0\t5\n")


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t5\t.\tA\tT\n", Columns.VCF),
        ("@HD\tVN:1.6\n@SQ\tSN:1\tLN:100\nr\t0\t1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n", Columns.SAM),
        ("##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=a\n", Columns.GFF),
        ('chr1\tsrc\texon\t1\t9\t.\t-\t0\tgene_id "a";\n', Columns.GFF),
        (
            "browser position chr1\ntrack name=a\n#c\nchr1\t0\t5\n",
            Columns(1, 2, 3, True, "#", skip_lines=2),
        ),
        ("chr1\t0\t5\n", Columns.BED),
    ],
)
def test_columns_are_inferred_from_content(tmp_path: Path, text: str, expected: Columns) -> None:
    path = tmp_path / "stream"
    with pybgzf.open(path, index=IndexFormat.TBI, columns=INFER) as handle:
        handle.write(text)
        handle.flush()
        writer = handle.buffer
        assert isinstance(writer, BgzfWriter)
        assert writer.columns == expected
    explicit = tmp_path / "explicit"
    with pybgzf.open(explicit, index=IndexFormat.TBI, columns=expected) as handle:
        handle.write(text)
    assert read_index(Path(f"{path}.tbi")) == read_index(Path(f"{explicit}.tbi"))


def test_columns_are_undecided_until_the_first_data_line(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "stream", index=IndexFormat.TBI, columns=INFER) as writer:
        writer.write(b"# a comment\n")
        assert writer.columns is None
        writer.write(b"chr1\t1\t2\n")
        assert writer.columns == Columns.BED


def test_ambiguous_content_raises(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "stream", index=IndexFormat.TBI, columns=INFER) as writer:
        with pytest.raises(ValueError, match="line 2 does not look like"):
            writer.write(b"# comment\nhello world\n")


def test_inference_without_data_raises_at_close(tmp_path: Path) -> None:
    writer = BgzfWriter(tmp_path / "stream", index=IndexFormat.TBI, columns=INFER)
    writer.write(b"# nothing but comments\n")
    with pytest.raises(ValueError, match="no data lines"):
        writer.close()
    assert gzip.decompress((tmp_path / "stream").read_bytes()) == b"# nothing but comments\n"


def read_fifo_in_thread(path: Path) -> tuple[threading.Thread, bytearray]:
    received = bytearray()

    def drain() -> None:
        with path.open("rb") as fifo:
            received.extend(fifo.read())

    reader = threading.Thread(target=drain)
    reader.start()
    return reader, received


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="named pipes are not available")
def test_a_fifo_needs_an_index_path(tmp_path: Path) -> None:
    fifo = tmp_path / "features.vcf"
    os.mkfifo(fifo)
    with pytest.raises(ValueError, match="not a regular file; pass index_path"):
        BgzfWriter(fifo, index=IndexFormat.TBI, columns=Columns.BED)
    reader, received = read_fifo_in_thread(fifo)
    data = bed_text().encode()
    index_path = tmp_path / "fifo.tbi"
    with BgzfWriter(fifo, index=IndexFormat.TBI, index_path=index_path, columns=INFER) as writer:
        writer.write(data)
    reader.join()
    assert writer.columns == Columns.BED
    assert gzip.decompress(bytes(received)) == data
    regular = tmp_path / "regular.bed.gz"
    write_bytes(regular, data, index=IndexFormat.TBI, columns=Columns.BED)
    assert bytes(received) == regular.read_bytes()
    assert read_index(index_path) == read_index(Path(f"{regular}.tbi"))
    assert not Path(f"{fifo}.tbi").exists()


@pytest.mark.skipif(not os.path.exists("/dev/null"), reason="/dev/null is not available")
def test_a_character_device_needs_an_index_path(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="'/dev/null', which is not a regular file"):
        BgzfWriter("/dev/null", index=IndexFormat.CSI, columns=Columns.BED)
    index_path = tmp_path / "null.csi"
    with BgzfWriter(
        "/dev/null", index=IndexFormat.CSI, index_path=index_path, columns=BED
    ) as writer:
        writer.write(b"chr1\t1\t2\n")
    assert read_index(index_path).magic == b"CSI\x01"


def test_a_regular_file_gets_the_default_index_path(tmp_path: Path) -> None:
    path = tmp_path / "exists.bed.gz"
    path.write_bytes(b"old")
    with BgzfWriter(path, index=IndexFormat.CSI, columns=Columns.BED) as writer:
        writer.write(b"chr1\t1\t2\n")
    assert Path(f"{path}.csi").exists()


@pytest.mark.parametrize(
    ("name", "text", "expected"),
    [
        ("points.bed.gz", "#c\nchr1\t5\nchr1\t9\n", Columns.BED2),
        ("points.bed.gz", "track name=x\nchr1\t5\n", Columns(1, 2, None, True, "#", skip_lines=1)),
        ("features.bed.gz", "chr1\t5\t6\n", Columns.BED),
        ("features.bed.gz", "chr1\t5\t6\tname\n", Columns.BED),
        ("stream", "chr1\t5\n", Columns.BED2),
    ],
)
def test_two_column_bed_is_inferred_as_bed2(
    tmp_path: Path, name: str, text: str, expected: Columns
) -> None:
    with pybgzf.open(tmp_path / name, index=IndexFormat.TBI, columns=INFER) as handle:
        handle.write(text)
        writer = handle.buffer
        assert isinstance(writer, BgzfWriter)
    assert writer.columns == expected


def test_bed_accepts_lines_without_an_end(tmp_path: Path) -> None:
    for columns in (Columns.BED, Columns.BED2):
        path = tmp_path / f"{columns.end}.bed.gz"
        with pybgzf.open(path, index=IndexFormat.TBI, columns=columns) as handle:
            handle.write("chr1\t5\nchr1\t9\t20\n")
        assert Path(f"{path}.tbi").exists()


@pytest.mark.parametrize(
    "header",
    [
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4611686018427387904>\n",
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=9223372036854775807>\n",
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=99999999999999999999>\n",
    ],
)
def test_references_too_long_for_csi_raise(tmp_path: Path, header: str) -> None:
    path = tmp_path / "long.vcf.gz"
    with BgzfWriter(path, index=IndexFormat.CSI, columns=Columns.VCF) as writer:
        writer.write(header.encode())
        with pytest.raises(ValueError, match="too long for a CSI index"):
            writer.write(b"chr1\t100\t.\tA\tT\t.\t.\t.\n")
    assert not Path(f"{path}.csi").exists()


def test_sam_references_too_long_for_csi_raise(tmp_path: Path) -> None:
    path = tmp_path / "long.sam.gz"
    with BgzfWriter(path, index=IndexFormat.CSI, columns=Columns.SAM) as writer:
        writer.write(b"@SQ\tSN:chr1\tLN:4611686018427387904\n")
        with pytest.raises(ValueError, match="too long for a CSI index"):
            writer.write(b"r\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n")


def test_huge_cigar_lengths_saturate(tmp_path: Path) -> None:
    with BgzfWriter(tmp_path / "a.sam.gz", index=IndexFormat.CSI, columns=Columns.SAM) as writer:
        with pytest.raises(ValueError, match="beyond the limit"):
            writer.write(b"r\t0\tchr1\t5\t60\t9223372036854775807M2M\t*\t0\t0\tA\tI\n")
