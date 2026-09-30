from collections.abc import Callable
from collections.abc import Iterator
from pathlib import Path

import pytest
from pybgzf import Columns
from pybgzf import LineFormat


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("a.bed", Columns.BED),
        ("a.bed.gz", Columns.BED),
        ("dir.vcf/A.BED.BGZ", Columns.BED),
        ("a.gff.bgzf", Columns.GFF),
        ("a.gff3.gz", Columns.GFF),
        ("a.gtf", Columns.GFF),
        ("a.vcf.gz", Columns.VCF),
        ("a.VCF", Columns.VCF),
        ("a.sam.gz", Columns.SAM),
    ],
)
def test_from_path(name: str, expected: Columns) -> None:
    assert Columns.from_path(name) == expected
    assert Columns.from_path(Path(name)) == expected


@pytest.mark.parametrize("name", ["a.txt", "a.gz", "a.bed.gz.gz", "bed", "a.bam", "a.tsv.gz"])
def test_from_path_rejects_unknown_suffixes(name: str) -> None:
    with pytest.raises(ValueError, match="pass columns explicitly"):
        Columns.from_path(name)


def test_sniff_accepts_str_and_bytes_with_terminators() -> None:
    assert Columns.sniff(["##fileformat=VCFv4.3\n"]) == Columns.VCF
    assert Columns.sniff([b"# x\r\n", b"chr1\t0\t1\r\n"]) == Columns.BED


def test_sniff_reads_no_further_than_needed() -> None:
    def lines() -> Iterator[str]:
        yield "@HD\tVN:1.6"
        raise AssertionError("read too far")

    assert Columns.sniff(lines()) == Columns.SAM


def test_sniff_without_data_raises() -> None:
    with pytest.raises(ValueError, match="no data line"):
        Columns.sniff(["# only", "track name=x"])


def test_sniff_with_ambiguous_data_raises() -> None:
    with pytest.raises(ValueError, match="line 1 does not look like"):
        Columns.sniff(["chr1 0 5"])


def test_presets() -> None:
    assert Columns.BED == Columns(1, 2, 3, True, "#")
    assert Columns.GFF == Columns(1, 4, 5, False, "#")
    assert Columns.VCF.format is LineFormat.VCF
    assert Columns.SAM.meta_char == "@"


@pytest.mark.parametrize(
    "make",
    [
        lambda: Columns(0, 2, 3, True, "#"),
        lambda: Columns(1, 2, 0, True, "#"),
        lambda: Columns(1, 2, 3, True, "##"),
        lambda: Columns(1, 2, 3, True, ""),
        lambda: Columns(1, 2, 3, False, "#", format=LineFormat.VCF),
        lambda: Columns(1, 2, None, True, "#", format=LineFormat.SAM),
    ],
)
def test_invalid_columns(make: Callable[[], Columns]) -> None:
    with pytest.raises(ValueError):
        make()
