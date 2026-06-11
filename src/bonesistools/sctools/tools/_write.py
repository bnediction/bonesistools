#!/usr/bin/env python

from pathlib import Path
from typing import Any, Optional, Union, cast

import pandas as pd
from anndata import AnnData
from scipy import io, sparse

PathInput = Union[str, Path]


def _get_matrix(adata: AnnData, layer: Optional[str]) -> Any:
    return adata.X if layer is None else adata.layers[layer]


def _get_dense_matrix(adata: AnnData, layer: Optional[str]) -> Any:
    X = _get_matrix(adata, layer)

    if sparse.issparse(X):
        return cast(Any, X).toarray()

    return cast(Any, X)


def to_csv(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData matrix to a CSV file.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix to export.
    filename: str or Path
        Output path. The `.csv` suffix is added if missing.
    layer: str, optional
        Layer to export instead of `adata.X`.
    """

    filename = str(filename)
    if not filename.endswith(".csv"):
        filename = Path(f"{filename}.csv")

    X = _get_dense_matrix(adata, layer)
    pd.DataFrame(X).to_csv(path_or_buf=filename, sep=",")


def to_mtx(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData matrix to a Matrix Market file.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix to export.
    filename: str or Path
        Output path. The `.mtx` suffix is added if missing.
    layer: str, optional
        Layer to export instead of `adata.X`.
    """

    filename = str(filename)
    if not filename.endswith(".mtx"):
        filename = Path(f"{filename}.mtx")

    X = _get_matrix(adata, layer)
    io.mmwrite(filename, X)


def to_npz(adata: AnnData, filename: PathInput, layer: Optional[str] = None) -> None:
    """
    Write an AnnData sparse matrix to a NumPy `.npz` file.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix to export.
    filename: str or Path
        Output path. The `.npz` suffix is added if missing.
    layer: str, optional
        Layer to export instead of `adata.X`.
    """

    filename = str(filename)
    if not filename.endswith(".npz"):
        filename = Path(f"{filename}.npz")

    X = _get_matrix(adata, layer)
    sparse.save_npz(filename, X)


def to_csv_or_mtx(
    adata: AnnData,
    filename: PathInput,
    layer: Optional[str] = None,
) -> None:
    """
    Write an AnnData matrix as CSV if dense, or Matrix Market if sparse.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix to export.
    filename: str or Path
        Output path. The suffix is added by the selected writer if missing.
    layer: str, optional
        Layer to export instead of `adata.X`.
    """

    if layer is None:
        if not sparse.issparse(adata.X):
            to_csv(adata, filename)
        else:
            to_mtx(adata, filename)
    else:
        if not sparse.issparse(adata.layers[layer]):
            to_csv(adata, filename, layer)
        else:
            to_mtx(adata, filename, layer)


def to_csv_or_npz(
    adata: AnnData,
    filename: PathInput,
    layer: Optional[str] = None,
) -> None:
    """
    Write an AnnData matrix as CSV if dense, or NumPy `.npz` if sparse.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix to export.
    filename: str or Path
        Output path. The suffix is added by the selected writer if missing.
    layer: str, optional
        Layer to export instead of `adata.X`.
    """

    if layer is None:
        if not sparse.issparse(adata.X):
            to_csv(adata, filename)
        else:
            to_npz(adata, filename)
    else:
        if not sparse.issparse(adata.layers[layer]):
            to_csv(adata, filename, layer)
        else:
            to_npz(adata, filename, layer)
