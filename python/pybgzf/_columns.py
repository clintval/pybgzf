from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from enum import auto
from pathlib import Path
from typing import Any
from typing import ClassVar
from typing import Final

from typing_extensions import override

from pybgzf import _pybgzf

ColumnsTuple = tuple[int, int, int | None, bool, str, int, _pybgzf.LineKind]

BGZF_SUFFIXES: Final[tuple[str, ...]] = (".gz", ".bgz", ".bgzf")
"""The file name suffixes of BGZF files, which any gzip reader can also read."""


class LineFormat(Enum):
    """How a line's end is found, which tabix also records as the file format."""

    @override
    def __reduce_ex__(self, proto: object) -> tuple[Any, ...]:
        """Pickle by name, since the value is an extension object."""
        return getattr, (type(self), self.name)

    GENERIC = _pybgzf.LineKind.GENERIC
    """The end is read from the end column, or the line covers one base if there is none."""

    SAM = _pybgzf.LineKind.SAM
    """The end is computed from the CIGAR string."""

    VCF = _pybgzf.LineKind.VCF
    """The end is computed from REF, `INFO/END`, `INFO/SVLEN`, and `FORMAT/LEN`."""


class Infer(Enum):
    """A sentinel that asks a writer to infer its columns while it streams."""

    INFER = auto()


INFER: Final = Infer.INFER
"""Infer columns from the destination's file name or, failing that, from the content."""


@dataclass(frozen=True, slots=True)
class Columns:
    """Where each line keeps its reference name, start, and end, like `tabix -s -b -e -0 -c -S`.

    Positions are read as base-10 integers, so, unlike in tabix, `010` is 10 and `0x10` is 0.

    Attributes:
        refname: The 1-based column holding the reference name.
        start: The 1-based column holding the start position.
        end: The 1-based column holding the end position, or None when each line covers one base.
        zero_based: True for 0-based, half-open positions (BED), False for 1-based, closed ones.
        meta_char: Lines starting with this character are headers and are not indexed.
        skip_lines: The number of leading lines that are never indexed.
        format: How the end is found; SAM and VCF compute it from other columns.
    """

    refname: int
    start: int
    end: int | None
    zero_based: bool
    meta_char: str
    skip_lines: int = 0
    format: LineFormat = LineFormat.GENERIC

    BED: ClassVar[Columns]
    BED2: ClassVar[Columns]
    GFF: ClassVar[Columns]
    VCF: ClassVar[Columns]
    SAM: ClassVar[Columns]

    def __post_init__(self) -> None:
        """Check that tabix could record this layout."""
        _pybgzf.validate_columns(columns_to_tuple(self))

    @classmethod
    def from_path(cls, path: str | os.PathLike[str]) -> Columns:
        """Choose columns from a file name, as `tabix` chooses a preset.

        A trailing `.gz`, `.bgz`, or `.bgzf` is ignored, then `.bed` is BED, `.gff`, `.gff3`, and
        `.gtf` are GFF, `.vcf` is VCF, and `.sam` is SAM, ignoring case.

        Raises:
            ValueError: If the file name has none of these suffixes.
        """
        name = Path(path).name.lower()
        name = next((name.removesuffix(s) for s in BGZF_SUFFIXES if name.endswith(s)), name)
        presets = {".bed": cls.BED, ".gff": cls.GFF, ".gff3": cls.GFF, ".gtf": cls.GFF}
        presets |= {".sam": cls.SAM, ".vcf": cls.VCF}
        if (preset := presets.get(Path(name).suffix)) is not None:
            return preset
        raise ValueError(
            f"cannot infer columns from the file name {os.fspath(path)!r}; pass columns explicitly"
        )

    @classmethod
    def sniff(cls, lines: Iterable[str | bytes]) -> Columns:
        """Choose columns from the first lines of a file, reading no further than needed.

        Header lines decide straight away: `##fileformat=VCF` is VCF, `@HD`, `@SQ`, `@RG`, `@PG`,
        and `@CO` are SAM, and `##gff-version` is GFF.
        Lines starting with `track `, `browser `, or `#` are skipped.
        Otherwise the first data line decides: GFF if it has nine fields that look like GFF or GTF,
        else BED if its second and third fields are ordered, non-negative integers, or BED2 if it
        has only two fields and the second is a non-negative integer.
        Leading `track ` and `browser ` lines of a BED or GFF file become `skip_lines`.

        Raises:
            ValueError: If the first data line looks like none of these, or no data line is found.
        """
        sniffer = _pybgzf.Sniffer()
        for line in lines:
            data = line.encode() if isinstance(line, str) else line
            decided = sniffer.push(data.removesuffix(b"\n").removesuffix(b"\r"))
            if decided is not None:
                return columns_from_tuple(decided)
        raise ValueError(
            "could not infer columns because no data line was found; pass columns explicitly"
        )


def columns_to_tuple(columns: Columns) -> ColumnsTuple:
    """Flatten columns for the extension module."""
    c = columns
    return (c.refname, c.start, c.end, c.zero_based, c.meta_char, c.skip_lines, c.format.value)


def columns_from_tuple(values: ColumnsTuple) -> Columns:
    """Rebuild columns from the extension module."""
    return Columns(*values[:-1], LineFormat(values[-1]))


Columns.BED = Columns(refname=1, start=2, end=3, zero_based=True, meta_char="#")
Columns.BED2 = Columns(refname=1, start=2, end=None, zero_based=True, meta_char="#")
Columns.GFF = Columns(refname=1, start=4, end=5, zero_based=False, meta_char="#")
Columns.VCF = Columns(
    refname=1, start=2, end=None, zero_based=False, meta_char="#", format=LineFormat.VCF
)
Columns.SAM = Columns(
    refname=3, start=4, end=None, zero_based=False, meta_char="@", format=LineFormat.SAM
)
