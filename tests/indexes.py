"""Read tabix and CSI indexes into comparable values, ignoring the order bins are stored in."""

import gzip
import struct
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Reference:
    bins: tuple[tuple[int, int, tuple[tuple[int, int], ...]], ...]
    linear: tuple[int, ...]


@dataclass(frozen=True)
class ParsedIndex:
    magic: bytes
    min_shift: int
    depth: int
    header: bytes
    references: tuple[Reference, ...]
    unplaced: int | None


class _Cursor:
    def __init__(self, data: bytes) -> None:
        self.data: bytes = data
        self.at: int = 0

    def take(self, fmt: str) -> tuple[int, ...]:
        values = struct.unpack_from(fmt, self.data, self.at)
        self.at += struct.calcsize(fmt)
        return values

    def one(self, fmt: str) -> int:
        return self.take(fmt)[0]

    def raw(self, n: int) -> bytes:
        chunk = self.data[self.at : self.at + n]
        self.at += n
        return chunk

    def done(self) -> bool:
        return self.at >= len(self.data)


def _bins(cursor: _Cursor, csi: bool) -> tuple[tuple[int, int, tuple[tuple[int, int], ...]], ...]:
    bins: list[tuple[int, int, tuple[tuple[int, int], ...]]] = []
    for _ in range(cursor.one("<i")):
        bin_id = cursor.one("<I")
        loff = cursor.one("<Q") if csi else 0
        count = cursor.one("<i")
        chunks = tuple(cursor.take("<QQ") for _ in range(count))
        bins.append((bin_id, loff, tuple((int(a), int(b)) for a, b in chunks)))
    return tuple(sorted(bins))


def read_index(path: Path) -> ParsedIndex:
    """Parse a BGZF-compressed tabix or CSI index."""
    cursor = _Cursor(gzip.decompress(path.read_bytes()))
    magic = cursor.raw(4)
    csi = magic == b"CSI\x01"
    if csi:
        min_shift, depth = cursor.take("<ii")
        header = cursor.raw(cursor.one("<i"))
        n_ref = cursor.one("<i")
    else:
        assert magic == b"TBI\x01", magic
        min_shift, depth = 14, 5
        n_ref = cursor.one("<i")
        fields = cursor.raw(28)
        names = cursor.raw(struct.unpack_from("<i", fields, 24)[0])
        header = fields + names
    references: list[Reference] = []
    for _ in range(n_ref):
        bins = _bins(cursor, csi)
        linear: tuple[int, ...] = ()
        if not csi:
            linear = tuple(cursor.take(f"<{cursor.one('<i')}Q"))
        references.append(Reference(bins, linear))
    unplaced = None if cursor.done() else cursor.one("<Q")
    assert cursor.done()
    return ParsedIndex(magic, min_shift, depth, header, tuple(references), unplaced)
