#!/usr/bin/env python

from __future__ import annotations

import math
import warnings
from typing import Any, Optional, Tuple, Union, cast, overload

import numpy as np
from anndata import AnnData
from scipy import sparse
from scipy.sparse import csr_matrix

from ..._compat import Literal
from ..._validation import (
    _as_boolean,
    _as_memory_size,
    _as_positive_number,
    _as_string,
)
from .._stats import _column_mean_variance
from .._typing import anndata_checker

ScaleMatrixSlot = Literal["X", "layer", "obsm"]


@overload
def normalize(
    adata: AnnData,
    target_sum: float = 1e4,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def normalize(
    adata: AnnData,
    target_sum: float = 1e4,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    *,
    copy: Literal[True],
) -> AnnData: ...


@overload
def normalize(
    adata: AnnData,
    target_sum: float = 1e4,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def normalize(
    adata: AnnData,
    target_sum: float = 1e4,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Normalize observations by library size.

    This function provides stable library-size normalization for AnnData
    matrices and sparse count matrices. It is intended to avoid small numerical
    differences caused by Scanpy backend changes across versions.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    target_sum: float (default: 1e4)
        Target total expression value per observation after normalization.
    expression: str, optional
        Expression matrix to normalize. If `None`, use `adata.X`; otherwise, use
        `adata.layers[expression]`.
    key_added: str, optional
        Layer key used to store normalized values. If `None`, overwrite the
        selected input expression matrix.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with normalized values added.
        Otherwise, updates `adata` in place and returns None.

        Normalized values are stored in:

        - `adata.X`: if `expression is None` and `key_added is None`;
        - `adata.layers[expression]`: if `expression` is provided and
          `key_added is None`;
        - `adata.layers[key_added]`: if `key_added` is provided.

    """

    target_sum = _as_positive_number(target_sum, "target_sum")
    expression = None if expression is None else _as_string(expression, "expression")
    key_added = None if key_added is None else _as_string(key_added, "key_added")
    copy = _as_boolean(copy, "copy")

    adata = adata.copy() if copy else adata
    matrix = _read_expression_matrix(adata, expression)

    if sparse.issparse(matrix):
        normalized_matrix = cast(csr_matrix, cast(Any, matrix).tocsr(copy=True))
        if not np.issubdtype(normalized_matrix.dtype, np.floating):
            normalized_matrix = normalized_matrix.astype(np.float32)

        counts = np.asarray(normalized_matrix.sum(axis=1)).ravel()
        counts = counts / target_sum
        counts = counts + (counts == 0)
        np.true_divide(
            normalized_matrix.data,
            np.repeat(counts, np.diff(normalized_matrix.indptr)),
            out=normalized_matrix.data,
        )
    else:
        normalized_matrix = np.asarray(matrix)
        if key_added is not None:
            normalized_matrix = normalized_matrix.copy()
        if not np.issubdtype(normalized_matrix.dtype, np.floating):
            normalized_matrix = normalized_matrix.astype(np.float32)

        counts = np.asarray(normalized_matrix.sum(axis=1)).ravel()
        counts = counts / target_sum
        counts = counts + (counts == 0)
        np.true_divide(normalized_matrix, counts[:, None], out=normalized_matrix)

    _write_expression_matrix(
        adata,
        normalized_matrix,
        expression=expression,
        key_added=key_added,
    )

    return adata if copy else None


@overload
def log1p(
    adata: AnnData,
    base: float = math.e,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    max_memory: Optional[Union[int, str]] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def log1p(
    adata: AnnData,
    base: float = math.e,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    max_memory: Optional[Union[int, str]] = None,
    *,
    copy: Literal[True],
) -> AnnData: ...


@overload
def log1p(
    adata: AnnData,
    base: float = math.e,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    max_memory: Optional[Union[int, str]] = None,
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def log1p(
    adata: AnnData,
    base: float = math.e,
    expression: Optional[str] = None,
    key_added: Optional[str] = None,
    max_memory: Optional[Union[int, str]] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Log-transform an expression matrix using log(X + 1).

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    base: float (default: math.e)
        Logarithm base. The transformation is `log(X + 1) / log(base)`.
    expression: str, optional
        Expression matrix to transform. If `None`, use `adata.X`; otherwise, use
        `adata.layers[expression]`.
    key_added: str, optional
        Layer key used to store transformed values. If `None`, overwrite the
        selected input expression matrix.
    max_memory: int or str, optional
        Approximate maximum memory allocated to dense working arrays. If `None`,
        no chunking is used. Integers are interpreted as bytes. Human-readable
        strings such as `"512MB"`, `"2GB"` or `"1GiB"` are accepted.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with transformed values added.
        Otherwise, updates `adata` in place and returns None.

        Transformed values are stored in:

        - `adata.X`: if `expression is None` and `key_added is None`;
        - `adata.layers[expression]`: if `expression` is provided and
          `key_added is None`;
        - `adata.layers[key_added]`: if `key_added` is provided.

    """

    base = _as_positive_number(base, "base")
    if base == 1:
        raise ValueError(
            f"invalid argument value for 'base': expected value != 1 "
            f"but received {base!r}"
        )
    expression = None if expression is None else _as_string(expression, "expression")
    key_added = None if key_added is None else _as_string(key_added, "key_added")
    max_memory_bytes = (
        None if max_memory is None else _as_memory_size(max_memory, "max_memory")
    )
    copy = _as_boolean(copy, "copy")

    adata = adata.copy() if copy else adata
    matrix = _read_expression_matrix(adata, expression)
    log_base = math.log(base)

    if sparse.issparse(matrix):
        if key_added is not None:
            transformed_matrix = cast(Any, matrix).copy()
        else:
            transformed_matrix = cast(Any, matrix)

        if not np.issubdtype(transformed_matrix.dtype, np.floating):
            transformed_matrix = transformed_matrix.astype(np.float32)

        np.log1p(transformed_matrix.data, out=transformed_matrix.data)
        if log_base != 1.0:
            np.true_divide(
                transformed_matrix.data,
                log_base,
                out=transformed_matrix.data,
            )
    else:
        transformed_matrix = np.asarray(matrix)
        if key_added is not None:
            transformed_matrix = transformed_matrix.copy()
        if not np.issubdtype(transformed_matrix.dtype, np.floating):
            transformed_matrix = transformed_matrix.astype(np.float32)

        if max_memory_bytes is None:
            row_chunk_size = max(1, int(transformed_matrix.shape[0]))
        else:
            bytes_per_row = max(
                1,
                int(transformed_matrix.shape[1] * transformed_matrix.dtype.itemsize),
            )
            if max_memory_bytes < bytes_per_row:
                warnings.warn(
                    "Requested memory budget is smaller than the memory required "
                    "for a single row. Computation will proceed with row_chunk_size=1.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            row_chunk_size = max(1, int(max_memory_bytes // bytes_per_row))

        for start in range(0, transformed_matrix.shape[0], row_chunk_size):
            end = min(start + row_chunk_size, transformed_matrix.shape[0])
            chunk = transformed_matrix[start:end, :]
            np.log1p(chunk, out=chunk)
            if log_base != 1.0:
                np.true_divide(chunk, log_base, out=chunk)

    _write_expression_matrix(
        adata,
        transformed_matrix,
        expression=expression,
        key_added=key_added,
    )

    return adata if copy else None


@overload
def scale(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    representation: Optional[str] = None,
    key_added: Optional[str] = None,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def scale(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    representation: Optional[str] = None,
    key_added: Optional[str] = None,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: Literal[True],
) -> AnnData: ...


@overload
def scale(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    representation: Optional[str] = None,
    key_added: Optional[str] = None,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def scale(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    representation: Optional[str] = None,
    key_added: Optional[str] = None,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Scale an expression matrix or representation to unit variance.

    Variables or representation dimensions are scaled column-wise. If
    `zero_center=True`, feature means are subtracted before dividing by the
    standard deviation. Sparse matrices are densified when `zero_center=True`.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    expression: str, optional
        Expression matrix to scale. If `None` and `representation` is None, use
        `adata.X`; otherwise use `adata.layers[expression]`.
    representation: str, optional
        Representation matrix in `adata.obsm` to scale.
    key_added: str, optional
        Key used to store scaled values. If `None`, overwrite the selected input
        matrix.
    zero_center: bool (default: True)
        If `True`, subtract feature means before variance scaling.
    max_value: float, optional
        Clip scaled values to this absolute value. If `None`, do not clip.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with scaled values added.
        Otherwise, updates `adata` in place and returns None.

        Scaled values are stored in:

        - `adata.X`: if `expression is None`, `representation is None`,
          and `key_added is None`;
        - `adata.layers[expression]`: if `expression` is provided and
          `key_added is None`;
        - `adata.layers[key_added]`: if scaling expression values and
          `key_added` is provided;
        - `adata.obsm[representation]`: if `representation` is provided and
          `key_added is None`;
        - `adata.obsm[key_added]`: if scaling a representation and `key_added`
          is provided.

    """

    expression = None if expression is None else _as_string(expression, "expression")
    representation = (
        None if representation is None else _as_string(representation, "representation")
    )
    if expression is not None and representation is not None:
        raise ValueError(
            "invalid argument combination: "
            "'expression' and 'representation' cannot be both specified"
        )

    key_added = None if key_added is None else _as_string(key_added, "key_added")
    zero_center = _as_boolean(zero_center, "zero_center")
    max_value = (
        None if max_value is None else _as_positive_number(max_value, "max_value")
    )
    copy = _as_boolean(copy, "copy")

    adata = adata.copy() if copy else adata
    matrix, slot, source_key = _read_scale_matrix(
        adata,
        expression=expression,
        representation=representation,
    )

    if sparse.issparse(matrix):
        sparse_matrix = cast(Any, matrix)
        if not np.issubdtype(sparse_matrix.dtype, np.floating):
            sparse_matrix = sparse_matrix.astype(np.float64)

        means, variances = _column_mean_variance(sparse_matrix)
        stds = np.sqrt(variances)
        stds[(stds == 0) | ~np.isfinite(stds)] = 1

        if zero_center:
            scaled_matrix = sparse_matrix.toarray().astype(np.float64, copy=False)
            np.subtract(scaled_matrix, means, out=scaled_matrix)
            np.true_divide(scaled_matrix, stds, out=scaled_matrix)
            if max_value is not None:
                np.clip(scaled_matrix, -max_value, max_value, out=scaled_matrix)
        else:
            scaled_sparse_matrix = cast(csr_matrix, sparse_matrix.tocsr(copy=True))
            np.true_divide(
                scaled_sparse_matrix.data,
                stds[scaled_sparse_matrix.indices],
                out=scaled_sparse_matrix.data,
            )
            if max_value is not None:
                np.clip(
                    scaled_sparse_matrix.data,
                    -max_value,
                    max_value,
                    out=scaled_sparse_matrix.data,
                )
            scaled_matrix = scaled_sparse_matrix
    else:
        scaled_matrix = np.asarray(matrix)
        if key_added is not None:
            scaled_matrix = scaled_matrix.copy()
        if not np.issubdtype(scaled_matrix.dtype, np.floating):
            scaled_matrix = scaled_matrix.astype(np.float64)

        dense_scaled_matrix = cast(np.ndarray, scaled_matrix)
        means, variances = _column_mean_variance(dense_scaled_matrix)
        stds = np.sqrt(variances)
        stds[(stds == 0) | ~np.isfinite(stds)] = 1

        if zero_center:
            np.subtract(dense_scaled_matrix, means, out=dense_scaled_matrix)
        np.true_divide(dense_scaled_matrix, stds, out=dense_scaled_matrix)
        if max_value is not None:
            np.clip(dense_scaled_matrix, -max_value, max_value, out=dense_scaled_matrix)
        scaled_matrix = dense_scaled_matrix

    matrix_shape = cast(Tuple[int, int], scaled_matrix.shape)
    if slot != "obsm" and matrix_shape[1] == adata.n_vars:
        adata.var["mean"] = means
        adata.var["std"] = stds

    _write_scale_matrix(
        adata,
        scaled_matrix,
        slot=slot,
        source_key=source_key,
        key_added=key_added,
    )

    return adata if copy else None


def _read_expression_matrix(
    adata: AnnData,
    expression: Optional[str],
) -> Any:

    return adata.X if expression is None else adata.layers[expression]


def _write_expression_matrix(
    adata: AnnData,
    matrix: Any,
    *,
    expression: Optional[str],
    key_added: Optional[str],
) -> None:

    if key_added is not None:
        adata.layers[key_added] = matrix
    elif expression is None:
        adata.X = matrix
    else:
        adata.layers[expression] = matrix


def _read_scale_matrix(
    adata: AnnData,
    *,
    expression: Optional[str],
    representation: Optional[str],
) -> Tuple[Any, ScaleMatrixSlot, Optional[str]]:

    if representation is not None:
        return adata.obsm[representation], "obsm", representation
    if expression is not None:
        return adata.layers[expression], "layer", expression

    return adata.X, "X", None


def _write_scale_matrix(
    adata: AnnData,
    matrix: Any,
    *,
    slot: ScaleMatrixSlot,
    source_key: Optional[str],
    key_added: Optional[str],
) -> None:

    if key_added is not None:
        if slot == "obsm":
            adata.obsm[key_added] = matrix
        else:
            adata.layers[key_added] = matrix
    elif slot == "obsm":
        adata.obsm[cast(str, source_key)] = matrix
    elif slot == "layer":
        adata.layers[cast(str, source_key)] = matrix
    else:
        adata.X = matrix
