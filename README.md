# pybgzf

[![Build Status](https://github.com/clintval/pybgzf/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/pybgzf/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/pybgzf)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/clintval/pybgzf/blob/main/LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Streaming BGZF compression with on-the-fly tabix and CSI indexing.

Install with pip or uv:

```console
pip install pybgzf
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
>>> with pybgzf.writer(path, index=IndexFormat.TBI, columns=Columns.BED) as handle:
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

Read the BED file back on several threads, or query a region with 0-based, half-open coordinates:

```python
>>> with pybgzf.reader(path, threads=4) as handle:
...     handle.readline()
'chr1\t100\t200\tgene-a\n'

>>> with pybgzf.IndexedReader(path) as reader:
...     list(reader.query("chr1", 180, 190))
['chr1\t100\t200\tgene-a', 'chr1\t150\t300\tgene-b']

```

## Development and Testing

See the [contributing guide](https://github.com/clintval/pybgzf/blob/main/CONTRIBUTING.md) for more information.

The multithreaded, position-tracking writer is adapted from [fgumi](https://github.com/fulcrumgenomics/fgumi), and indexing follows [htslib](https://github.com/samtools/htslib); see [NOTICE](https://github.com/clintval/pybgzf/blob/main/NOTICE).
