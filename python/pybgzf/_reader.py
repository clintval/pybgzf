from __future__ import annotations

import codecs
import errno
import io
import os
from collections.abc import Iterator
from types import TracebackType
from typing import Protocol

from typing_extensions import Buffer
from typing_extensions import Self
from typing_extensions import override

from pybgzf import _pybgzf
from pybgzf._columns import Columns
from pybgzf._columns import columns_from_tuple

_TEXT_BUFFER_SIZE = 1 << 20
_MAX_POSITION = (1 << 63) - 1


class ReadableBinary(Protocol):
    """Anything with a `read` method that returns bytes, such as `sys.stdin.buffer`."""

    def read(self, size: int, /) -> bytes:
        """Read up to `size` bytes."""
        ...


class BgzfReader(io.RawIOBase):
    """A binary file that decompresses BGZF as it is read.

    Positions are BGZF virtual offsets, the compressed offset of a block shifted left 16 bits
    plus the offset into that block's decompressed data, as stored in tabix and CSI indexes.
    `tell()` returns one and `seek()` accepts one, but `seekable()` is False because virtual
    offsets cannot be added to or subtracted from like byte offsets.
    Data that ends without the BGZF end-of-file marker is read, with a `TruncatedWarning`.

    Args:
        src: A path to open, or a readable binary file-like object such as a pipe.
        threads: The number of threads decompressing blocks, from 1 to 1024; 1 decompresses in
            the calling thread.

    Raises:
        ValueError: If `threads` is not between 1 and 1024.
        OSError: If `src` cannot be opened.
    """

    def __init__(self, src: str | os.PathLike[str] | ReadableBinary, *, threads: int = 1) -> None:
        super().__init__()
        if isinstance(src, (str, os.PathLike)):
            source: str | ReadableBinary = os.fspath(src)
            self._name: object = source
        elif callable(getattr(src, "read", None)):
            source = src
            self._name = getattr(src, "name", None)
        else:
            raise TypeError(f"expected a path or an object with a read method, not {src!r}")
        self._inner: _pybgzf.Reader = _pybgzf.Reader(source, threads=threads)

    @property
    def name(self) -> object:
        """The path read from, or the name of the file object read from, if it has one."""
        return self._name

    @override
    def readable(self) -> bool:
        """Return True: this file is readable."""
        return True

    @override
    def readinto(self, buffer: Buffer, /) -> int:
        """Read into a writable buffer and return the number of bytes read."""
        view = memoryview(buffer).cast("B")
        data = self._inner.read(len(view))
        view[: len(data)] = data
        return len(data)

    @override
    def readall(self) -> bytes:
        """Read and return everything left."""
        return self._inner.readall()

    @override
    def readline(self, size: int | None = -1, /) -> bytes:
        """Read and return the next line, including its newline."""
        return self._inner.readline(-1 if size is None else size)

    @override
    def tell(self) -> int:
        """Return the virtual offset of the next byte to be read."""
        return self._inner.tell()

    @override
    def seek(self, offset: int, whence: int = io.SEEK_SET, /) -> int:
        """Move to a virtual offset from `tell()` or an index, and return it.

        Only `io.SEEK_SET` is supported, and only when the source can seek.

        Raises:
            ValueError: If the offset is in no block of the file and is not its end.
        """
        if whence != io.SEEK_SET:
            raise io.UnsupportedOperation("only seeking to a virtual offset is supported")
        if offset < 0:
            raise ValueError(f"negative seek position {offset}")
        if offset >= 1 << 64:
            raise ValueError(f"virtual offset {offset} is too large")
        return self._inner.seek(offset)

    @override
    def close(self) -> None:
        """Stop any decompression threads and close the source if it was opened from a path.

        With more than one thread, a pipe or file-like source may be read from once more in the
        background after this returns, if a read was already waiting for data.
        """
        if self.closed:
            return
        inner: _pybgzf.Reader | None = getattr(self, "_inner", None)
        try:
            if inner is not None:
                inner.close()
        finally:
            super().close()


def reader(
    src: str | os.PathLike[str] | ReadableBinary,
    *,
    threads: int = 1,
    encoding: str | None = None,
    errors: str | None = None,
    newline: str | None = None,
) -> io.TextIOWrapper:
    """Open a BGZF file for reading text.

    The encoding defaults to UTF-8, and `errors` and `newline` work as in `io.TextIOWrapper`.
    """
    encoding = codecs.lookup("utf-8" if encoding is None else encoding).name
    raw = BgzfReader(src, threads=threads)
    try:
        buffered = io.BufferedReader(raw, buffer_size=_TEXT_BUFFER_SIZE)
        return io.TextIOWrapper(buffered, encoding=encoding, errors=errors, newline=newline)
    except BaseException:
        raw.close()
        raise


class IndexedReader:
    """A BGZF file and its tabix or CSI index, for reading the lines in a region.

    Lines are parsed with the columns and header character recorded in the index, and a query
    returns what `tabix path ref:start+1-end` prints, in the same order.
    A file that ends without the BGZF end-of-file marker warns with `TruncatedWarning` when opened.

    Args:
        path: The BGZF file.
        index_path: Its index; defaults to `path` plus `.csi` or, if there is none, `.tbi`, the
            order htslib looks in.
        threads: The number of threads decompressing blocks, from 1 to 1024; 1 decompresses in
            the calling thread and is usually fastest for many small queries, since more threads
            restart their read-ahead at every seek.

    Raises:
        FileNotFoundError: If the file, or an index for it, cannot be found.
        ValueError: If the index is not a tabix or CSI index with a tabix header.
    """

    def __init__(
        self,
        path: str | os.PathLike[str],
        *,
        index_path: str | os.PathLike[str] | None = None,
        threads: int = 1,
    ) -> None:
        file = os.fspath(path)
        if index_path is None:
            candidates = [f"{file}.csi", f"{file}.tbi"]
            found = [candidate for candidate in candidates if os.path.exists(candidate)]
            if not found:
                raise FileNotFoundError(
                    errno.ENOENT,
                    f"no index next to {file!r}: expected {candidates[0]!r} or {candidates[1]!r}; "
                    + "pass index_path",
                )
            index_path = found[0]
        self._inner: _pybgzf.IndexedReader = _pybgzf.IndexedReader(
            file, index_path, threads=threads
        )

    @property
    def refnames(self) -> list[str]:
        """The reference names in the index, in the order they appear in the file.

        Raises:
            ValueError: If a name is not UTF-8.
        """
        return self._inner.refnames

    @property
    def columns(self) -> Columns:
        """The columns, header character, and skipped lines recorded in the index."""
        return columns_from_tuple(self._inner.columns)

    @property
    def closed(self) -> bool:
        """True once the reader is closed."""
        return self._inner.closed

    def query(self, refname: str, start: int, end: int) -> Iterator[str]:
        """Return the lines, without their newlines, overlapping `[start, end)` on `refname`.

        Coordinates are 0-based and half-open, as in BED.
        An empty region, or a reference name not in the index, returns no lines.

        Raises:
            ValueError: If `start` is negative or `end` is less than `start`, or, once the lines
                before it have been returned, a line in the region cannot be parsed or is not
                UTF-8, naming its virtual offset.
        """
        if start < 0 or end < start:
            raise ValueError(
                f"start must be at least 0 and end at least start, not {start} and {end}"
            )
        return self._inner.query(refname, min(start, _MAX_POSITION), min(end, _MAX_POSITION))

    def close(self) -> None:
        """Stop any decompression threads and close the file."""
        self._inner.close()

    def __enter__(self) -> Self:
        """Return this reader."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Close this reader."""
        self.close()
