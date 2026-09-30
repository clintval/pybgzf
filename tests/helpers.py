"""Helpers shared by the tests."""

import random
import shutil
import subprocess
from pathlib import Path

import pytest

HAS_HTSLIB = shutil.which("tabix") is not None and shutil.which("bgzip") is not None

requires_htslib = pytest.mark.skipif(not HAS_HTSLIB, reason="tabix and bgzip are not installed")


def htslib_version() -> tuple[int, ...]:
    """Return the version of the installed tabix, or an empty tuple if there is none."""
    if not HAS_HTSLIB:
        return ()
    result = subprocess.run(["tabix", "--version"], check=True, capture_output=True, text=True)
    version = result.stdout.splitlines()[0].split()[-1]
    return tuple(int(part) for part in version.split("+")[0].split(".") if part.isdigit())


requires_htslib_1_23 = pytest.mark.skipif(
    htslib_version() < (1, 23), reason="VCF end positions follow htslib 1.23 and later"
)


def tabix(*args: str | Path) -> str:
    """Run tabix and return its standard output."""
    result = subprocess.run(["tabix", *map(str, args)], check=True, capture_output=True, text=True)
    return result.stdout


def bed_lines(seed: int = 7, references: tuple[str, ...] = ("chr1", "chr10", "chr2")) -> list[str]:
    """Sorted BED lines spanning many BGZF blocks, with widths from 1 base to 300 kb."""
    rng = random.Random(seed)
    lines: list[str] = []
    for name in references:
        position = 0
        for index in range(rng.randint(3000, 9000)):
            position += rng.randint(0, 2500)
            width = rng.choice([1, 40, 700, 20_000, 300_000])
            lines.append(
                f"{name}\t{position}\t{position + width}\tfeature{index}\t{rng.random()}\n"
            )
    return lines


def bed_text() -> str:
    """A BED file with a header line."""
    return "#chrom\tstart\tend\tname\tscore\n" + "".join(bed_lines())
