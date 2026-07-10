#!/usr/bin/env python

from pathlib import Path
from typing import Optional, Union

from anndata import AnnData

from ..._warnings import _warn_deprecated
from ..input_output import _write as _io_write

PathInput = Union[str, Path]


def to_csv(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData matrix to a CSV file.

    Deprecated. Use `bt.omics.io.to_csv` instead.
    """

    _warn_deprecated(
        "`bt.omics.tl.to_csv()`",
        replacement="`bt.omics.io.to_csv()`",
        stacklevel=2,
    )
    _io_write.to_csv(adata, filename, layer=layer)


def to_mtx(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData matrix to a Matrix Market file.

    Deprecated. Use `bt.omics.io.to_mtx` instead.
    """

    _warn_deprecated(
        "`bt.omics.tl.to_mtx()`",
        replacement="`bt.omics.io.to_mtx()`",
        stacklevel=2,
    )
    _io_write.to_mtx(adata, filename, layer=layer)


def to_npz(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData sparse matrix to a NumPy `.npz` file.

    Deprecated. Use `bt.omics.io.to_npz` instead.
    """

    _warn_deprecated(
        "`bt.omics.tl.to_npz()`",
        replacement="`bt.omics.io.to_npz()`",
        stacklevel=2,
    )
    _io_write.to_npz(adata, filename, layer=layer)
