from typing import final

from typing_extensions import Buffer

from pybgzf._columns import ColumnsTuple
from pybgzf._writer import WritableBinary

BLOCK_SIZE: int

@final
class Writer:
    def __new__(
        cls,
        dest: str | WritableBinary,
        *,
        level: int,
        threads: int,
        index: str | None,
        index_path: str | None,
        columns: ColumnsTuple | None,
        infer: bool,
        csi_min_shift: int,
        csi_depth: int | None,
    ) -> Writer: ...
    def write(self, data: Buffer) -> int: ...
    def flush(self) -> None: ...
    def tell(self) -> int: ...
    def close(self) -> None: ...
    @property
    def closed(self) -> bool: ...
    @property
    def columns(self) -> ColumnsTuple | None: ...

@final
class Sniffer:
    def __new__(cls) -> Sniffer: ...
    def push(self, line: bytes) -> ColumnsTuple | None: ...

def validate_columns(columns: ColumnsTuple) -> None: ...
