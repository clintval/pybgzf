"""Compare output with htslib's `bgzip` and `tabix`, which index the same compressed bytes."""

import gzip
import random
import shutil
import subprocess
from pathlib import Path

import pybgzf
import pytest
from pybgzf import INFER
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexFormat

from tests.helpers import bed_lines
from tests.helpers import bed_text
from tests.helpers import requires_htslib
from tests.helpers import requires_htslib_1_23
from tests.helpers import tabix
from tests.indexes import read_index

pytestmark = requires_htslib

TBI = IndexFormat.TBI
CSI = IndexFormat.CSI


def write(
    path: Path,
    text: str,
    index: IndexFormat,
    columns: Columns,
    *,
    threads: int = 1,
    chunk: int = 9_973,
    csi_min_shift: int = 14,
) -> Path:
    data = text.encode()
    with BgzfWriter(
        path, threads=threads, index=index, columns=columns, csi_min_shift=csi_min_shift
    ) as writer:
        for start in range(0, len(data), chunk):
            writer.write(data[start : start + chunk])
    return Path(f"{path}.{index.value}")


def htslib_index(path: Path, index: IndexFormat, *args: str) -> Path:
    """Index a copy of `path` with tabix and return the index path."""
    copy = path.parent / "htslib" / path.name
    copy.parent.mkdir(exist_ok=True)
    shutil.copyfile(path, copy)
    tabix("-f", *(["-C"] if index is CSI else []), *args, copy)
    return Path(f"{copy}.{index.value}")


def assert_identical(ours: Path, theirs: Path) -> None:
    """Assert two indexes are the same, first as parsed values and then byte for byte."""
    assert read_index(ours) == read_index(theirs)
    assert gzip.decompress(ours.read_bytes()) == gzip.decompress(theirs.read_bytes())


def assert_same_index(
    path: Path, text: str, index: IndexFormat, columns: Columns, *args: str
) -> None:
    ours = write(path, text, index, columns)
    assert_identical(ours, htslib_index(path, index, *args))


@pytest.mark.parametrize("index", [TBI, CSI])
@pytest.mark.parametrize("threads", [1, 4])
def test_bed(tmp_path: Path, index: IndexFormat, threads: int) -> None:
    path = tmp_path / "a.bed.gz"
    ours = write(path, bed_text(), index, Columns.BED, threads=threads)
    assert_identical(ours, htslib_index(path, index, "-p", "bed"))


def test_decompressed_tabix_index_is_byte_identical(tmp_path: Path) -> None:
    path = tmp_path / "a.bed.gz"
    ours = write(path, bed_text(), TBI, Columns.BED)
    theirs = htslib_index(path, TBI, "-p", "bed")
    assert gzip.decompress(ours.read_bytes()) == gzip.decompress(theirs.read_bytes())


def test_csi_min_shift(tmp_path: Path) -> None:
    path = tmp_path / "a.bed.gz"
    ours = write(path, bed_text(), CSI, Columns.BED, csi_min_shift=12)
    assert_identical(ours, htslib_index(path, CSI, "-p", "bed", "-m", "12"))


def test_decompresses_with_bgzip(tmp_path: Path) -> None:
    text = bed_text()
    path = tmp_path / "a.bed.gz"
    write(path, text, TBI, Columns.BED, threads=3)
    subprocess.run(["bgzip", "-t", str(path)], check=True)
    result = subprocess.run(["bgzip", "-dc", str(path)], check=True, capture_output=True)
    assert result.stdout == text.encode()


def test_region_queries_match(tmp_path: Path) -> None:
    lines = bed_lines(seed=11)
    path = tmp_path / "q.bed.gz"
    write(path, "".join(lines), TBI, Columns.BED, threads=2)
    copy_index = htslib_index(path, TBI, "-p", "bed")
    rng = random.Random(3)
    for _ in range(40):
        name = rng.choice(["chr1", "chr10", "chr2", "chr3"])
        start = rng.randint(1, 20_000_000)
        end = start + rng.choice([1, 100, 10_000, 1_000_000])
        region = f"{name}:{start}-{end}"
        expected = [
            line
            for line in lines
            if line.split("\t")[0] == name
            and int(line.split("\t")[1]) < end
            and int(line.split("\t")[2]) > start - 1
        ]
        assert tabix(path, region) == "".join(expected)
        assert tabix(copy_index.parent / path.name, region) == "".join(expected)


GFF = """##gff-version 3
chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1
chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=g1
chr1\tsrc\texon\t700\t900\t.\t+\t.\tParent=g1
chr1\tsrc\tgene\t70000\t2000000\t.\t-\t.\tID=g2
chr2\tsrc\tgene\t1\t1\t.\t.\t.\tID=g3
"""


@pytest.mark.parametrize("index", [TBI, CSI])
def test_gff(tmp_path: Path, index: IndexFormat) -> None:
    assert_same_index(tmp_path / "a.gff.gz", GFF, index, Columns.GFF, "-p", "gff")


VCF = """##fileformat=VCFv4.3
##contig=<ID=1,length=249250621>
##contig=<ID=X,length=156040895>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t0\t.\tN\tA\t.\t.\t.\tGT\t0/1\t0/0
1\t100\t.\tACGT\tA\t.\t.\t.\tGT\t0/1\t0/0
1\t150\tdel\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=5000\tGT\t0/1\t0/0
1\t160\tsv\tA\t<DUP:TANDEM>,T\t.\t.\tSVLEN=40000,.\tGT\t0/1\t0/0
1\t170\tbad\tA\t<DEL>\t.\t.\tEND=10\tGT\t0/1\t0/0
1\t200\tblock\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:300\t0/0:7000
1\t9000\tnolen\tA\t<NON_REF>\t.\t.\t.\tGT\t0/0\t0/0
X\t1\t.\tA\tG\t.\t.\t.\tGT\t1/1\t0/1
X\t156040000\t.\tA\tG\t.\t.\t.\tGT\t1/1\t0/1
"""


@requires_htslib_1_23
@pytest.mark.parametrize("index", [TBI, CSI])
def test_vcf(tmp_path: Path, index: IndexFormat) -> None:
    assert_same_index(tmp_path / "a.vcf.gz", VCF, index, Columns.VCF, "-p", "vcf")


@requires_htslib_1_23
def test_vcf_csi_depth_follows_contig_lengths(tmp_path: Path) -> None:
    long = VCF.replace("length=249250621", "length=90000000000").replace(
        "X\t156040000", "X\t100000000"
    )
    assert_same_index(tmp_path / "a.vcf.gz", long, CSI, Columns.VCF, "-p", "vcf")
    assert read_index(tmp_path / "a.vcf.gz.csi").depth == 8


@requires_htslib_1_23
@pytest.mark.parametrize(("length", "min_shift"), [(2**60, 34), (2**62 - 256, 35)])
def test_csi_min_shift_grows_for_long_contigs(tmp_path: Path, length: int, min_shift: int) -> None:
    text = VCF.replace("length=249250621", f"length={length}")
    ours = write(tmp_path / "a.vcf.gz", text, CSI, Columns.VCF)
    assert (read_index(ours).min_shift, read_index(ours).depth) == (min_shift, 9)
    assert_identical(ours, htslib_index(tmp_path / "a.vcf.gz", CSI, "-p", "vcf"))


SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:248956422
@SQ\tSN:chr2\tLN:242193529
r1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\tIIII
r2\t0\tchr1\t120\t60\t10M20000N10M\t*\t0\t0\tACGT\tIIII
r3\t16\tchr1\t130\t60\t5S10M3D5I10M\t*\t0\t0\tACGT\tIIII
r4\t0\tchr2\t5\t60\t*\t*\t0\t0\tACGT\tIIII
r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
"""


@pytest.mark.parametrize("index", [TBI, CSI])
def test_sam(tmp_path: Path, index: IndexFormat) -> None:
    assert_same_index(tmp_path / "a.sam.gz", SAM, index, Columns.SAM, "-p", "sam")


POINTS = "".join(
    f"c{name}\t{position}\tx\n" for name in (1, 2) for position in range(1, 200_000, 37)
)


def test_one_based_point_columns(tmp_path: Path) -> None:
    columns = Columns(refname=1, start=2, end=None, zero_based=False, meta_char="#")
    assert_same_index(tmp_path / "p.gz", POINTS, TBI, columns, "-s", "1", "-b", "2", "-e", "2")


def test_zero_based_point_columns(tmp_path: Path) -> None:
    columns = Columns(refname=1, start=2, end=None, zero_based=True, meta_char="#")
    assert_same_index(
        tmp_path / "p.gz", POINTS, TBI, columns, "-0", "-s", "1", "-b", "2", "-e", "2"
    )


def test_reordered_columns_meta_char_and_skip_lines(tmp_path: Path) -> None:
    text = "first line\nsecond\tline\n%comment\n" + "".join(
        f"{position}\t{position + 10}\tchr{name}\n"
        for name in (1, 2)
        for position in range(1, 100_000, 97)
    )
    columns = Columns(refname=3, start=1, end=2, zero_based=False, meta_char="%", skip_lines=2)
    args = ("-s", "3", "-b", "1", "-e", "2", "-c", "%", "-S", "2")
    assert_same_index(tmp_path / "r.gz", text, TBI, columns, *args)


def test_crlf(tmp_path: Path) -> None:
    text = bed_text().replace("\n", "\r\n")
    assert_same_index(tmp_path / "crlf.bed.gz", text, TBI, Columns.BED, "-p", "bed")


def test_unterminated_last_line(tmp_path: Path) -> None:
    text = bed_text().removesuffix("\n")
    assert_same_index(tmp_path / "a.bed.gz", text, TBI, Columns.BED, "-p", "bed")


@pytest.mark.parametrize("index", [TBI, CSI])
@pytest.mark.parametrize("text", ["", "#only a header\n"])
def test_files_without_records(tmp_path: Path, index: IndexFormat, text: str) -> None:
    assert_same_index(tmp_path / "e.bed.gz", text, index, Columns.BED, "-p", "bed")


def test_lines_ending_on_block_boundaries(tmp_path: Path) -> None:
    lines = [f"chr1\t{position:>12}\t{position + 5:>12}\t{'x' * 32}\n" for position in range(5000)]
    assert {len(line) for line in lines} == {64}
    for threads in (1, 3):
        path = tmp_path / f"b{threads}.bed.gz"
        ours = write(path, "".join(lines), TBI, Columns.BED, threads=threads, chunk=64)
        assert_identical(ours, htslib_index(path, TBI, "-p", "bed"))


def test_flushes_between_lines(tmp_path: Path) -> None:
    path = tmp_path / "f.bed.gz"
    with BgzfWriter(path, threads=2, index=TBI, columns=Columns.BED) as writer:
        for number, line in enumerate(bed_lines()[:20_000]):
            writer.write(line.encode())
            if number % 997 == 0:
                writer.flush()
    ours = Path(f"{path}.tbi")
    assert_identical(ours, htslib_index(path, TBI, "-p", "bed"))


@requires_htslib_1_23
def test_inferred_columns_index_like_tabix(tmp_path: Path) -> None:
    path = tmp_path / "stream"
    with pybgzf.open_writer(
        path, index=TBI, index_path=tmp_path / "stream.tbi", columns=INFER
    ) as handle:
        handle.write(VCF)
    assert read_index(tmp_path / "stream.tbi") == read_index(htslib_index(path, TBI, "-p", "vcf"))


QUERY_REGIONS = [
    ("chr1", 0, 1),
    ("chr1", 99, 100),
    ("chr1", 120, 20_200),
    ("chr1", 800, 70_000),
    ("chr2", 0, 10),
    ("1", 0, 200),
    ("1", 160, 170),
    ("1", 4_000, 10_000),
    ("X", 156_039_000, 156_041_000),
]


def assert_same_queries(path: Path, text: str, columns: Columns, *args: str) -> None:
    ours = write(path, text, TBI, columns)
    theirs = htslib_index(path, TBI, *args)
    for index_path in (ours, theirs):
        with pybgzf.IndexedReader(path, index_path=index_path) as reader:
            for refname, start, end in QUERY_REGIONS:
                expected = tabix(path, f"{refname}:{start + 1}-{end}").splitlines()
                assert list(reader.query(refname, start, end)) == expected, (refname, start, end)


def test_gff_queries_match_tabix(tmp_path: Path) -> None:
    assert_same_queries(tmp_path / "a.gff.gz", GFF, Columns.GFF, "-p", "gff")


def test_sam_queries_match_tabix(tmp_path: Path) -> None:
    assert_same_queries(tmp_path / "a.sam.gz", SAM, Columns.SAM, "-p", "sam")


@requires_htslib_1_23
def test_vcf_queries_match_tabix(tmp_path: Path) -> None:
    assert_same_queries(tmp_path / "a.vcf.gz", VCF, Columns.VCF, "-p", "vcf")


BED2 = "#chrom\tposition\n" + "".join(
    f"chr{name}\t{position}\n" for name in (1, 2) for position in range(0, 3_000_000, 173)
)
MIXED_BED = "".join(
    f"chr1\t{position}\n" if position % 3 else f"chr1\t{position}\t{position + 500}\n"
    for position in range(0, 3_000_000, 211)
)


@pytest.mark.parametrize("text", [BED2, MIXED_BED], ids=["bed2", "mixed"])
def test_bed_without_an_end_column_is_byte_identical(tmp_path: Path, text: str) -> None:
    path = tmp_path / "points.bed.gz"
    ours = write(path, text, TBI, Columns.BED)
    theirs = htslib_index(path, TBI, "-p", "bed")
    assert gzip.decompress(ours.read_bytes()) == gzip.decompress(theirs.read_bytes())


@pytest.mark.parametrize("index", [TBI, CSI])
def test_bed2_matches_tabix_point_columns(tmp_path: Path, index: IndexFormat) -> None:
    path = tmp_path / "points.bed.gz"
    ours = write(path, BED2, index, Columns.BED2)
    theirs = htslib_index(path, index, "-0", "-s", "1", "-b", "2", "-e", "2")
    assert_identical(ours, theirs)
    with pybgzf.IndexedReader(path, index_path=ours) as reader:
        assert reader.columns == Columns.BED2
        assert list(reader.query("chr1", 173, 174)) == ["chr1\t173"]
        assert list(reader.query("chr1", 174, 346)) == []
        assert list(reader.query("chr2", 0, 347)) == ["chr2\t0", "chr2\t173", "chr2\t346"]
        assert list(reader.query("chr1", 173, 174)) == tabix(path, "chr1:174-174").splitlines()
