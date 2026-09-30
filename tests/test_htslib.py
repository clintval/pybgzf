"""Compare output with htslib's `bgzip` and `tabix`, which index the same compressed bytes."""

import bisect
import gzip
import itertools
import random
import shutil
import struct
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path

import pybgzf
import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from pybgzf import INFER
from pybgzf import BgzfWriter
from pybgzf import Columns
from pybgzf import IndexFormat

from tests.helpers import bed_lines
from tests.helpers import bed_text
from tests.helpers import htslib_version
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
    return Path(f"{path}.{index.name.lower()}")


def htslib_index(path: Path, index: IndexFormat, *args: str) -> Path:
    """Index a copy of `path` with tabix and return the index path."""
    copy = path.parent / "htslib" / path.name
    copy.parent.mkdir(exist_ok=True)
    shutil.copyfile(path, copy)
    tabix("-f", *(["-C"] if index is CSI else []), *args, copy)
    return Path(f"{copy}.{index.name.lower()}")


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


@requires_htslib_1_23
def test_vcf_records_wider_than_a_block(tmp_path: Path) -> None:
    samples = 40_000
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
    lines = [
        "##fileformat=VCFv4.3\n",
        "##contig=<ID=1,length=1000000>\n",
        header + "".join(f"\tS{sample}" for sample in range(samples)) + "\n",
        *(
            f"1\t{position}\t.\tA\tG\t.\t.\t.\tGT" + "\t0/0" * samples + "\n"
            for position in range(100, 105)
        ),
    ]
    path = tmp_path / "wide.vcf.gz"
    with BgzfWriter(path, index=TBI, columns=Columns.VCF) as writer:
        for line in lines:
            writer.write(line.encode())
    assert_identical(Path(f"{path}.tbi"), htslib_index(path, TBI, "-p", "vcf"))


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
    with pybgzf.writer(
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


def random_bed_stream(rng: random.Random) -> bytes:
    """Sorted BED text with zero-length and bin-boundary features, long lines, and comments."""
    lines = ["#chrom\tstart\tend\n"]
    for reference in range(rng.choice([1, 3, 60])):
        position = 0
        for _ in range(rng.randint(1, 300)):
            position += rng.choice([0, 1, 16384 - position % 16384, rng.randint(0, 50_000)])
            width = rng.choice([0, 1, 16384, rng.randint(1, 2_000_000)])
            padding = "x" * rng.choice([0, 0, 0, 70_000])
            lines.append(f"c{reference}\t{position}\t{position + width}\t{padding}\n")
            if rng.random() < 0.02:
                lines.append("#between records\n")
    text = "".join(lines)
    if rng.random() < 0.3:
        text = text.replace("\n", "\r\n")
    if rng.random() < 0.3:
        text = text.rstrip("\r\n")
    return text.encode()


@pytest.mark.parametrize("seed", range(8))
def test_random_streams_index_like_tabix(tmp_path: Path, seed: int) -> None:
    rng = random.Random(seed)
    data = random_bed_stream(rng)
    index = rng.choice([TBI, CSI])
    for threads in (1, 3):
        path = tmp_path / f"s{threads}.bed.gz"
        with BgzfWriter(path, threads=threads, index=index, columns=Columns.BED) as writer:
            at = 0
            while at < len(data):
                size = rng.choice([1, 7, 1000, 65280, 300_000])
                writer.write(data[at : at + size])
                at += size
                if rng.random() < 0.05:
                    writer.flush()
        ours = Path(f"{path}.{index.name.lower()}")
        assert_identical(ours, htslib_index(path, index, "-p", "bed"))


def test_the_tabix_limit_matches_tabix(tmp_path: Path) -> None:
    path = tmp_path / "edge.bed.gz"
    assert_same_index(path, "chr1\t1\t536870912\n", TBI, Columns.BED, "-p", "bed")
    with BgzfWriter(tmp_path / "over.bed.gz", index=TBI, columns=Columns.BED) as writer:
        with pytest.raises(ValueError, match="beyond the tabix limit"):
            writer.write(b"chr1\t1\t536870913\n")


TABIX_SKIPS_COMMENTS = htslib_version() >= (1, 23)
SYMBOLIC = ["<DEL>", "<DUP>", "<DUP:TANDEM>", "<INV>", "<CNV>", "<INS>", "<*>", "<NON_REF>"]

Place = Callable[[str, int], str]


@st.composite
def reference_positions(draw: st.DrawFn, limit: int) -> list[int]:
    """Sorted positions from a reference start, a bin edge, or halfway to `limit`."""
    first = draw(st.sampled_from([0, 1, 5_000, 16_384, limit // 2]))
    step = st.one_of(st.integers(0, 3_000), st.sampled_from([0, 16_384, 1 << 16]))
    steps = draw(st.lists(step, min_size=1, max_size=20))
    count = draw(st.integers(1_000, 3_000))
    return list(
        itertools.accumulate(itertools.islice(itertools.cycle(steps), count), initial=first)
    )


@st.composite
def placed_lines(draw: st.DrawFn, lines: st.SearchStrategy[Place], limit: int) -> list[str]:
    """Lines drawn from `lines` and placed, in turn, on three references."""
    places = draw(st.lists(lines, min_size=4, max_size=20))
    return [
        place(name, position)
        for name in ["chr1", "chr2", "chr10"]
        for position, place in zip(
            draw(reference_positions(limit)), itertools.cycle(places), strict=False
        )
    ]


@st.composite
def vcf_line(draw: st.DrawFn) -> Place:
    """A VCF record ending at its REF, END, SVLEN, or FORMAT LEN, to be placed."""
    ref = draw(st.text("ACGT", min_size=1, max_size=60))
    alleles = draw(st.lists(st.sampled_from(["A", "T", ".", *SYMBOLIC]), min_size=1, max_size=3))
    info = draw(st.lists(st.sampled_from(["SVTYPE=DEL", "XEND=9", "CIEND=-5,5"]), max_size=2))
    if draw(st.booleans()):
        lengths = st.one_of(st.just("."), st.integers(-200_000, 200_000).map(str))
        info.append(f"SVLEN={','.join(draw(st.lists(lengths, min_size=1, max_size=3)))}")
    end = draw(st.one_of(st.none(), st.just("."), st.integers(-5, 200_000).map(str)))
    at = draw(st.integers(0, len(info)))
    keys = draw(st.sampled_from(["GT", "GT:LEN", "LEN:GT"]))
    length = st.one_of(st.just("."), st.integers(-5, 20_000).map(str))
    samples = [
        ":".join(draw(length) if key == "LEN" else "0/1" for key in keys.split(":"))
        for _ in range(2)
    ]

    def place(name: str, position: int) -> str:
        entries = list(info)
        if end is not None:
            entries.insert(at, f"END={end if end == '.' else position + int(end)}")
        fields = [name, str(position), ".", ref, ",".join(alleles), ".", "PASS"]
        return "\t".join([*fields, ";".join(entries) or ".", keys, *samples]) + "\n"

    return place


@st.composite
def vcf_file(draw: st.DrawFn, limit: int) -> str:
    """A VCF with SNVs, indels, symbolic alleles, and gVCF blocks on three contigs."""
    lines = draw(placed_lines(vcf_line(), limit))
    longest = max(int(line.split("\t")[1]) for line in lines) + 1_000_000
    header = ["##fileformat=VCFv4.3\n"]
    header += [f"##contig=<ID={name},length={longest}>\n" for name in ["chr1", "chr2", "chr10"]]
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
    return "".join(header + lines)


CIGAR_OPERATIONS = st.one_of(
    st.builds("{}{}".format, st.integers(1, 150), st.sampled_from("MIDSH=X")),
    st.builds("{}N".format, st.integers(50, 200_000)),
)
CIGARS = st.one_of(st.just("*"), st.lists(CIGAR_OPERATIONS, min_size=1, max_size=6).map("".join))


@st.composite
def sam_line(draw: st.DrawFn) -> Place:
    """A SAM record, spliced, clipped, or unmapped, to be placed."""
    flag = draw(st.sampled_from([0, 16, 99, 147, 4, 69, 133]))
    cigar = "*" if flag & 4 else draw(CIGARS)

    def place(name: str, position: int) -> str:
        return f"r\t{flag}\t{name}\t{position}\t60\t{cigar}\t*\t0\t0\tACGT\tIIII\n"

    return place


@st.composite
def sam_file(draw: st.DrawFn, limit: int) -> str:
    """A SAM with spliced, clipped, and unmapped reads, and unplaced reads at the end."""
    lines = ["@HD\tVN:1.6\tSO:coordinate\n"]
    lines += [f"@SQ\tSN:{name}\tLN:{limit}\n" for name in ["chr1", "chr2", "chr10"]]
    lines += draw(placed_lines(sam_line(), limit))
    lines += ["u\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"] * draw(st.integers(0, 30))
    return "".join(lines)


@st.composite
def gff_line(draw: st.DrawFn) -> Place:
    """A GFF3 feature up to 400 kb wide, to be placed, and maybe a comment after it."""
    kind = draw(st.sampled_from(["gene", "exon", "CDS"]))
    strand = draw(st.sampled_from("+-."))
    width = draw(
        st.one_of(st.sampled_from([0, 1, 10, 300, 5_000, 400_000]), st.integers(0, 400_000))
    )
    comment = draw(st.sampled_from(["", "###\n", "# a comment\n"])) if TABIX_SKIPS_COMMENTS else ""

    def place(name: str, position: int) -> str:
        fields = [name, "src", kind, str(position + 1), str(position + width), ".", strand, "."]
        return "\t".join([*fields, "ID=f"]) + "\n" + comment

    return place


@st.composite
def gff_file(draw: st.DrawFn, limit: int) -> str:
    """A GFF3 with features up to 400 kb wide, and comments between them if tabix skips them."""
    return "##gff-version 3\n" + "".join(draw(placed_lines(gff_line(), limit)))


def uncompressed_block_starts(path: Path) -> list[int]:
    """Return the uncompressed offset at which every BGZF block of `path` starts."""
    data = path.read_bytes()
    starts: list[int] = []
    at = total = 0
    while at < len(data):
        starts.append(total)
        at += struct.unpack_from("<H", data, at + 16)[0] + 1
        total += struct.unpack_from("<I", data, at - 4)[0]
    return starts


def data_lines(text: str, columns: Columns) -> list[tuple[int, str, int]]:
    """Return each data line's offset, reference name, and 0-based start."""
    found: list[tuple[int, str, int]] = []
    offset = 0
    for line in text.splitlines(keepends=True):
        if not line.startswith(columns.meta_char):
            fields = line.split("\t")
            start = int(fields[columns.start - 1]) - (0 if columns.zero_based else 1)
            found.append((offset, fields[columns.refname - 1], max(start, 0)))
        offset += len(line)
    return found


@st.composite
def query_regions(
    draw: st.DrawFn, rows: list[tuple[int, str, int]], blocks: list[int]
) -> list[tuple[str, int, int]]:
    """Regions at reference starts and ends, block boundaries, and random lines."""
    last = {name: start for _, name, start in rows}
    offsets = [offset for offset, _, _ in rows]
    anchors = [(name, 0) for name in last] + [(name, start + 1) for name, start in last.items()]
    for block in blocks[1:]:
        at = max(bisect.bisect_right(offsets, block) - 1, 0)
        anchors += [(name, start) for _, name, start in rows[at : at + 2]]
    nudges = st.one_of(st.sampled_from([-1, 0, 1]), st.integers(2, 1_000))
    for _, name, start in draw(st.lists(st.sampled_from(rows), min_size=100, max_size=300)):
        anchors.append((name, max(start + draw(nudges), 0)))
    found = [(name, start, 1 << 40) for name, start in last.items()]
    found += [("chrUn", 0, 1 << 29), ("CHR1", 0, 1 << 29), ("chr", 0, 100)]
    widths = st.sampled_from([0, 1, 2, 100, 16_384, 1_000_000])
    for name, position in anchors:
        width = draw(widths)
        found += [(name, position, position + width), (name, max(position - width, 0), position)]
    return draw(st.permutations(found))


def assert_queries_match_tabix(
    path: Path, readers: list[pybgzf.IndexedReader], regions: list[tuple[str, int, int]]
) -> None:
    """Query all regions with one tabix call, and region by region only to report a difference.

    Tabix reads `ref:1-0` as the whole reference, so regions ending at 0 must be empty instead.
    """
    empty = [region for region in regions if region[2] == 0]
    queried = [region for region in regions if region[2] > 0]
    arguments = [f"{name}:{start + 1}-{end}" for name, start, end in queried]
    expected = tabix(path, *arguments).splitlines()
    for reader in readers:
        assert not any(list(reader.query(*region)) for region in empty)
        found = [line for name, start, end in queried for line in reader.query(name, start, end)]
        if found != expected:
            for (name, start, end), argument in zip(queried, arguments, strict=True):
                lines = list(reader.query(name, start, end))
                assert lines == tabix(path, argument).splitlines(), argument
        assert found == expected


RANDOM_FILES = {
    "vcf": (vcf_file, Columns.VCF),
    "sam": (sam_file, Columns.SAM),
    "gff": (gff_file, Columns.GFF),
}


@pytest.mark.parametrize("index", [TBI, CSI])
@pytest.mark.parametrize("preset", [pytest.param("vcf", marks=requires_htslib_1_23), "sam", "gff"])
@settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data=st.data())
def test_random_files_index_and_query_like_tabix(
    preset: str, index: IndexFormat, data: st.DataObject
) -> None:
    generate, columns = RANDOM_FILES[preset]
    text = data.draw(generate(1 << 29 if index is TBI else 1 << 31), label="text")
    sizes = data.draw(
        st.lists(st.sampled_from([1, 7, 1000, 65280, 300_000]), min_size=1), label="sizes"
    )
    threads = data.draw(st.sampled_from([1, 3]), label="threads")
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / f"r.{preset}.gz"
        with pybgzf.writer(path, threads=threads, index=index, columns=columns) as handle:
            at = 0
            for size in itertools.cycle(sizes):
                if at >= len(text):
                    break
                handle.write(text[at : at + size])
                at += size
        ours = Path(f"{path}.{index.name.lower()}")
        theirs = htslib_index(path, index, "-p", preset)
        assert_identical(ours, theirs)
        blocks = uncompressed_block_starts(path)
        assert len(blocks) > 1
        found = data.draw(query_regions(data_lines(text, columns), blocks), label="regions")
        with (
            pybgzf.IndexedReader(path, index_path=ours) as first,
            pybgzf.IndexedReader(path, index_path=theirs, threads=3) as second,
        ):
            assert_queries_match_tabix(theirs.parent / path.name, [first, second], found)
