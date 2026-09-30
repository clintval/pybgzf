"""Compare pybgzf with the standard library's gzip and, if installed, htslib's bgzip and tabix.

Writing compresses the data in 1 MiB chunks; reading iterates over the lines of the result.

Run with `uv run python benchmarks/benchmark.py [megabytes]`.
"""

import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

import pybgzf


def bed_data(megabytes: int) -> bytes:
    rng = random.Random(42)
    lines: list[str] = []
    size = 0
    references = iter(f"chr{number}" for number in range(1, 10_000))
    name = next(references)
    position = 0
    while size < megabytes * 1_000_000:
        position += rng.randint(0, 400)
        if position > 50_000_000:
            name, position = next(references), 0
        width = rng.choice([1, 50, 300, 5000])
        score = rng.randint(0, 1000)
        line = f"{name}\t{position}\t{position + width}\tfeature{len(lines)}\t{score}\t+\n"
        lines.append(line)
        size += len(line)
    return "".join(lines).encode()


def timed(label: str, action: Callable[[], None], size: int) -> None:
    best = float("inf")
    for _ in range(3):
        start = time.perf_counter()
        action()
        best = min(best, time.perf_counter() - start)
    print(f"{label:<42} {best:8.3f} s {size / best / 1e6:9.1f} MB/s")


def main() -> None:
    megabytes = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    data = bed_data(megabytes)
    cores = os.cpu_count() or 1
    print(f"{len(data) / 1e6:.0f} MB of BED on {cores} cores, level 6, best of 3\n")
    with tempfile.TemporaryDirectory() as directory:
        write_benchmarks(data, Path(directory), cores)
        print()
        read_benchmarks(data, Path(directory), cores)


def write_benchmarks(data: bytes, out: Path, cores: int) -> None:

    def stdlib_gzip() -> None:
        with gzip.open(out / "stdlib.gz", "wb", compresslevel=6) as handle:
            handle.write(data)

    timed("gzip (stdlib)", stdlib_gzip, len(data))

    for threads in sorted({1, 2, 4, 8, cores}):
        for index in (None, pybgzf.IndexFormat.TBI):

            def write(threads: int = threads, index: pybgzf.IndexFormat | None = index) -> None:
                columns = pybgzf.Columns.BED if index else None
                with pybgzf.BgzfWriter(
                    out / "pybgzf.bed.gz", threads=threads, index=index, columns=columns
                ) as writer:
                    for start in range(0, len(data), 1 << 20):
                        writer.write(data[start : start + (1 << 20)])

            label = f"pybgzf threads={threads}" + (" + tabix index" if index else "")
            timed(label, write, len(data))

    (out / "input.bed").write_bytes(data)
    if shutil.which("bgzip"):
        for threads in sorted({1, 2, 4, 8, cores}):

            def bgzip(threads: int = threads) -> None:
                with (out / "bgzip.bed.gz").open("wb") as handle:
                    subprocess.run(
                        ["bgzip", "-l", "6", "-@", str(threads), "-c", str(out / "input.bed")],
                        stdout=handle,
                        check=True,
                    )

            timed(f"bgzip -@{threads}", bgzip, len(data))

            def bgzip_and_tabix(threads: int = threads) -> None:
                bgzip(threads)
                subprocess.run(["tabix", "-f", "-p", "bed", str(out / "bgzip.bed.gz")], check=True)

            if shutil.which("tabix"):
                timed(f"bgzip -@{threads} then tabix", bgzip_and_tabix, len(data))


def read_benchmarks(data: bytes, out: Path, cores: int) -> None:
    with pybgzf.BgzfWriter(out / "read.bed.gz", threads=cores) as writer:
        writer.write(data)
    path = out / "read.bed.gz"

    def gzip_lines() -> None:
        with gzip.open(path, "rt") as handle:
            for _ in handle:
                pass

    timed("read lines, gzip (stdlib)", gzip_lines, len(data))

    for threads in sorted({1, 2, 4, 8, cores}):

        def pybgzf_lines(threads: int = threads) -> None:
            with pybgzf.open_reader(path, threads=threads) as handle:
                for _ in handle:
                    pass

        timed(f"read lines, pybgzf threads={threads}", pybgzf_lines, len(data))

        def pybgzf_bytes(threads: int = threads) -> None:
            with pybgzf.BgzfReader(path, threads=threads) as reader:
                reader.readall()

        timed(f"read bytes, pybgzf threads={threads}", pybgzf_bytes, len(data))

    if shutil.which("bgzip"):
        for threads in sorted({1, 4, cores}):

            def bgzip_decompress(threads: int = threads) -> None:
                subprocess.run(
                    ["bgzip", "-d", "-c", "-@", str(threads), str(path)],
                    stdout=subprocess.DEVNULL,
                    check=True,
                )

            timed(f"read bytes, bgzip -d -@{threads}", bgzip_decompress, len(data))


if __name__ == "__main__":
    main()
