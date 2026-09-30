"""Streaming BGZF compression with on-the-fly tabix and CSI indexing."""

from pybgzf._columns import INFER
from pybgzf._columns import Columns
from pybgzf._columns import Infer
from pybgzf._columns import LineFormat
from pybgzf._writer import BgzfWriter
from pybgzf._writer import IndexFormat
from pybgzf._writer import WritableBinary
from pybgzf._writer import open  # noqa: A004

__all__ = [
    "INFER",
    "BgzfWriter",
    "Columns",
    "IndexFormat",
    "Infer",
    "LineFormat",
    "WritableBinary",
    "open",
]
