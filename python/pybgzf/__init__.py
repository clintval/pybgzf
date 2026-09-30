"""Streaming BGZF compression with on-the-fly tabix and CSI indexing."""

import importlib.metadata

from pybgzf._columns import INFER
from pybgzf._columns import Columns
from pybgzf._columns import Infer
from pybgzf._columns import LineFormat
from pybgzf._reader import BgzfReader
from pybgzf._reader import IndexedReader
from pybgzf._reader import ReadableBinary
from pybgzf._reader import open_reader
from pybgzf._writer import BgzfWriter
from pybgzf._writer import IndexFormat
from pybgzf._writer import WritableBinary
from pybgzf._writer import open_writer

__version__ = importlib.metadata.version("pybgzf")

__all__ = [
    "INFER",
    "BgzfReader",
    "BgzfWriter",
    "Columns",
    "IndexFormat",
    "IndexedReader",
    "Infer",
    "LineFormat",
    "ReadableBinary",
    "WritableBinary",
    "open_reader",
    "open_writer",
]
