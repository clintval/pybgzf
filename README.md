# pybgzf

[![Build Status](https://github.com/clintval/pybgzf/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/pybgzf/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/pybgzf)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Streaming BGZF compression with on-the-fly tabix and CSI indexing.

Install with pip or uv:

```bash
❯ pip install pybgzf
```

Until the first release is on PyPI, install from GitHub, which builds from source with a [Rust toolchain](https://rustup.rs):

```bash
❯ pip install git+https://github.com/clintval/pybgzf
```

## Introduction

`pybgzf` writes blocked gzip (BGZF) files from Python and builds their tabix (`.tbi`) or CSI (`.csi`) index while it writes, so there is no second pass with `tabix`.
Its core is Rust: blocks are compressed with the multithreaded, libdeflate-backed [`bgzf`](https://github.com/fulcrumgenomics/bgzf) crate and indexes are written with [noodles](https://github.com/zaeleus/noodles).
It needs neither pysam nor htslib.

## Quick Start

Write sorted BED lines as text and the index appears next to the file on close:

```python
>>> from pathlib import Path
>>> from tempfile import mkdtemp
>>>
>>> import pybgzf
>>> from pybgzf import Columns, IndexFormat
>>>
>>> directory = Path(mkdtemp())
>>> path = directory / "features.bed.gz"
>>>
>>> with pybgzf.open(path, index=IndexFormat.TBI, columns=Columns.BED) as handle:
...     _ = handle.write("chr1\t100\t200\tgene-a\n")
...     _ = handle.write("chr1\t150\t300\tgene-b\n")
...     _ = handle.write("chr2\t10\t20\tgene-c\n")
>>>
>>> sorted(file.name for file in directory.iterdir())
['features.bed.gz', 'features.bed.gz.tbi']

```

The output is ordinary BGZF, so `gzip -d`, `bgzip -d`, and `tabix features.bed.gz chr1:120-160` all read it.

Stream bytes to a pipe or any binary file object, with columns inferred from the content:

```python
>>> import io
>>>
>>> stream = io.BytesIO()
>>> with pybgzf.BgzfWriter(stream, threads=4, index=IndexFormat.CSI, index_path=directory / "calls.csi", columns=pybgzf.INFER) as writer:
...     _ = writer.write(b"##fileformat=VCFv4.3\n")
...     writer.columns == Columns.VCF
True

```

## Features

- Indexes are identical to what `tabix` builds from the same file, for BED, GFF, VCF, SAM, and custom columns (`tabix -s -b -e -0 -c -S`).
- Compression runs on any number of threads and the output is byte-for-byte the same for every thread count.
- Lines may be split across writes, so `csv` and other text writers work through `pybgzf.open`.
- Unsorted input raises `ValueError` from the write that completes the offending line, naming its line number.
- `columns=pybgzf.INFER` picks BED, GFF, VCF, or SAM from the file name or, when streaming, from the first lines.
- The GIL is released while compressing and while waiting on compression threads.

## Benchmarks

Writing 200 MB of BED lines at level 6 on an Apple M3 Max with `benchmarks/benchmark.py`:

| Command                    | Seconds | MB/s  |
|----------------------------|--------:|------:|
| `gzip` (stdlib)            |   4.80  |    42 |
| pybgzf, 1 thread, `.tbi`   |   1.84  |   109 |
| pybgzf, 4 threads, `.tbi`  |   0.43  |   466 |
| pybgzf, 16 threads, `.tbi` |   0.31  |   641 |
| `bgzip -@1` then `tabix`   |   3.01  |    66 |
| `bgzip -@4` then `tabix`   |   1.31  |   153 |
| `bgzip -@16` then `tabix`  |   0.90  |   222 |

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.

The multithreaded, position-tracking writer is adapted from [fgumi](https://github.com/fulcrumgenomics/fgumi), and indexing follows [htslib](https://github.com/samtools/htslib); see [NOTICE](NOTICE).
