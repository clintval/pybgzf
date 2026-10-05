import os
from collections.abc import Iterator
from typing import ClassVar
from typing import final

from typing_extensions import Buffer
from typing_extensions import override

from pybgzf._columns import ColumnsTuple
from pybgzf._reader import ReadableBinary
from pybgzf._writer import WritableBinary

__all__ = [
    "BLOCK_SIZE",
    "IndexKind",
    "IndexedReader",
    "LineKind",
    "QueryIterator",
    "Reader",
    "Sniffer",
    "TruncatedWarning",
    "Writer",
    "validate_columns",
]

BLOCK_SIZE: int

class TruncatedWarning(UserWarning):
    """A BGZF file ends without its end-of-file marker, so it may be truncated."""

@final
class IndexKind:
    TBI: ClassVar[IndexKind]
    CSI: ClassVar[IndexKind]

@final
class LineKind:
    GENERIC: ClassVar[LineKind]
    SAM: ClassVar[LineKind]
    VCF: ClassVar[LineKind]

@final
class Writer:
    def __new__(
        cls,
        dest: str | WritableBinary,
        *,
        level: int,
        threads: int,
        index: IndexKind | None,
        index_path: str | None,
        columns: ColumnsTuple | None,
        infer: bool,
        infer_bed: bool,
        csi_min_shift: int,
        csi_depth: int | None,
    ) -> Writer: ...
    def write(self, data: Buffer) -> int: ...
    def flush(self) -> None: ...
    def tell(self) -> int: ...
    def close(self) -> None: ...
    def abandon(self) -> None: ...
    @property
    def closed(self) -> bool: ...
    @property
    def columns(self) -> ColumnsTuple | None: ...

@final
class Reader:
    def __new__(cls, src: str | ReadableBinary, *, threads: int) -> Reader: ...
    def read(self, size: int) -> bytes: ...
    def readall(self) -> bytes: ...
    def readline(self, size: int) -> bytes: ...
    def tell(self) -> int: ...
    def seek(self, position: int) -> int: ...
    def close(self) -> None: ...
    @property
    def closed(self) -> bool: ...

@final
class IndexedReader:
    def __new__(
        cls, path: str, index_path: str | os.PathLike[str], *, threads: int
    ) -> IndexedReader: ...
    def query(self, refname: str, start: int, end: int) -> QueryIterator: ...
    def close(self) -> None: ...
    @property
    def refnames(self) -> list[str]: ...
    @property
    def columns(self) -> ColumnsTuple: ...
    @property
    def closed(self) -> bool: ...

@final
class QueryIterator(Iterator[str]):
    @override
    def __iter__(self) -> QueryIterator: ...
    @override
    def __next__(self) -> str: ...

@final
class Sniffer:
    def __new__(cls) -> Sniffer: ...
    def push(self, line: bytes) -> ColumnsTuple | None: ...

def validate_columns(columns: ColumnsTuple) -> None: ...
