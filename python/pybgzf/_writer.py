from __future__ import annotations

import io
import os
import stat
from enum import Enum
from typing import Protocol

from typing_extensions import Buffer
from typing_extensions import override

from pybgzf import _pybgzf
from pybgzf._columns import Columns
from pybgzf._columns import Infer
from pybgzf._columns import columns_from_tuple
from pybgzf._columns import columns_to_tuple


class IndexFormat(Enum):
    """The kind of index to build while writing."""

    TBI = _pybgzf.IndexKind.TBI
    """A tabix index, for references up to 2^29 bases."""

    CSI = _pybgzf.IndexKind.CSI
    """A coordinate-sorted index, for references of any length."""


class WritableBinary(Protocol):
    """Anything with a `write` method that accepts bytes, such as `sys.stdout.buffer`."""

    def write(self, data: bytes, /) -> object:
        """Write bytes."""
        ...


_INDEX_SUFFIXES = {IndexFormat.TBI: ".tbi", IndexFormat.CSI: ".csi"}


def _is_regular_file(path: str) -> bool:
    """Return True if `path` is, or will be created as, a regular file."""
    try:
        return stat.S_ISREG(os.stat(path).st_mode)
    except FileNotFoundError:
        return True


def _suffix_columns(*names: str) -> Columns | None:
    for name in names:
        try:
            return Columns.from_path(name)
        except ValueError:
            continue
    return None


class BgzfWriter(io.RawIOBase):
    """A binary file that compresses what is written to it as BGZF.

    When `index` is set, every complete line is parsed with `columns` and added to a tabix or CSI
    index, which is written to `index_path` on close.
    Lines must be sorted by start within each reference, and each reference must be contiguous.
    If a line cannot be indexed, or columns cannot be inferred, no index is written and the file is
    left without the BGZF end-of-file marker, so that readers see it as truncated.
    Otherwise closing, including leaving a `with` block because of an exception or the writer being
    garbage collected, finishes the file and writes the index for what was written.

    Args:
        dest: A path to create, or a writable binary file-like object such as a pipe.
        level: The compression level, from 0 (stored) to 12 (smallest); 1 to 12 use libdeflate.
        threads: The number of threads compressing blocks, from 1 to 1024; 1 compresses in the
            calling thread.
        index: The kind of index to build, or None to build none.
        index_path: Where to write the index; defaults to `dest` plus `.tbi` or `.csi`, and is
            required when `index` is set and `dest` is a file-like object.
        columns: Where each line keeps its reference, start, and end; required when `index` is set.
            `INFER` infers them from the file name, or else from the content as it streams; a
            BED file whose first data line has two fields is BED2.
            Without `index`, columns are not used.
        csi_min_shift: The width, as a power of two, of the smallest CSI bin.
        csi_depth: The number of CSI bin levels; None chooses as `tabix -C` does.

    Raises:
        ValueError: If the options are invalid, checked before anything is created.
        OSError: If `dest` cannot be created.
    """

    def __init__(
        self,
        dest: str | os.PathLike[str] | WritableBinary,
        *,
        level: int = 6,
        threads: int = 1,
        index: IndexFormat | None = None,
        index_path: str | os.PathLike[str] | None = None,
        columns: Columns | Infer | None = None,
        csi_min_shift: int = 14,
        csi_depth: int | None = None,
    ) -> None:
        super().__init__()
        if isinstance(dest, (str, os.PathLike)):
            path: str | None = os.fspath(dest)
            sink: str | WritableBinary = os.fspath(dest)
        else:
            path = None
            sink = dest
        regular = path is not None and _is_regular_file(path)
        if index is None and index_path is not None:
            raise ValueError("index_path is only used when index is set")
        if index is not None and columns is None:
            raise ValueError("columns is required when index is set")
        if index is not None and index_path is None:
            if path is None or not regular:
                described = repr(path) if path is not None else "a file-like object"
                raise ValueError(
                    f"cannot place an index next to {described}, which is not a regular file; "
                    + "pass index_path"
                )
            index_path = f"{path}{_INDEX_SUFFIXES[index]}"
        index_file = None if index_path is None else os.fspath(index_path)

        infer = isinstance(columns, Infer)
        infer_bed = False
        explicit: Columns | None = columns if isinstance(columns, Columns) else None
        if infer and index_file is not None:
            index_name = index_file
            for suffix in _INDEX_SUFFIXES.values():
                index_name = index_name.removesuffix(suffix)
            names = [path, index_name] if path is not None and regular else [index_name]
            explicit = _suffix_columns(*names)
            infer = explicit is None
            if explicit == Columns.BED:
                explicit, infer_bed = None, True

        self._columns: Columns | None = explicit
        self._name: object = path if path is not None else getattr(dest, "name", None)
        self._inner: _pybgzf.Writer = _pybgzf.Writer(
            sink,
            level=level,
            threads=threads,
            index=None if index is None else index.value,
            index_path=index_file,
            columns=None if explicit is None or index is None else columns_to_tuple(explicit),
            infer=infer and index is not None,
            infer_bed=infer_bed and index is not None,
            csi_min_shift=csi_min_shift,
            csi_depth=csi_depth,
        )

    @property
    def name(self) -> object:
        """The path written to, or the name of the file object written to, if it has one."""
        return self._name

    @property
    def columns(self) -> Columns | None:
        """The columns in use, or None while they are still being inferred."""
        if self._columns is None:
            decided = self._inner.columns
            if decided is not None:
                self._columns = columns_from_tuple(decided)
        return self._columns

    @override
    def writable(self) -> bool:
        """Return True: this file is writable."""
        return True

    @override
    def write(self, data: Buffer, /) -> int:
        """Write bytes and return how many were written, which is always all of them.

        Raises:
            ValueError: If the file is closed, or a line cannot be indexed, naming the line.
                After an indexing error, further writes raise.
        """
        if self.closed:
            raise ValueError("I/O operation on closed file.")
        return self._inner.write(data)

    @override
    def flush(self) -> None:
        """End the current BGZF block and flush everything written so far."""
        inner: _pybgzf.Writer | None = getattr(self, "_inner", None)
        if inner is None or inner.closed:
            if self.closed:
                raise ValueError("I/O operation on closed file.")
            return
        inner.flush()

    @override
    def tell(self) -> int:
        """Return the virtual offset of the next byte, waiting for pending blocks to be written."""
        if self.closed:
            raise ValueError("I/O operation on closed file.")
        return self._inner.tell()

    @override
    def close(self) -> None:
        """Write the remaining data, the BGZF end-of-file marker, and then the index.

        After an indexing error, the end-of-file marker and the index are not written.
        A file-like `dest` is flushed but not closed.

        Raises:
            ValueError: If the last line cannot be indexed, or columns could not be inferred.
        """
        if self.closed:
            return
        inner: _pybgzf.Writer | None = getattr(self, "_inner", None)
        try:
            if inner is not None:
                inner.close()
        finally:
            super().close()


def open_writer(
    dest: str | os.PathLike[str] | WritableBinary,
    *,
    encoding: str | None = None,
    errors: str | None = None,
    newline: str | None = None,
    level: int = 6,
    threads: int = 1,
    index: IndexFormat | None = None,
    index_path: str | os.PathLike[str] | None = None,
    columns: Columns | Infer | None = None,
    csi_min_shift: int = 14,
    csi_depth: int | None = None,
) -> io.TextIOWrapper:
    """Open a BGZF file for writing text, so that `csv` and other text writers can write to it.

    The encoding defaults to UTF-8 and `errors` works as in `io.TextIOWrapper`.
    Lines end in a line feed on every platform unless `newline` is given, which then works as in
    `io.TextIOWrapper`.
    Every other argument is passed to `BgzfWriter`.
    """
    writer = BgzfWriter(
        dest,
        level=level,
        threads=threads,
        index=index,
        index_path=index_path,
        columns=columns,
        csi_min_shift=csi_min_shift,
        csi_depth=csi_depth,
    )
    return io.TextIOWrapper(
        writer,
        encoding="utf-8" if encoding is None else encoding,
        errors=errors,
        newline="\n" if newline is None else newline,
    )
