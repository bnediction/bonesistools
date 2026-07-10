#!/usr/bin/env python

from __future__ import annotations

import warnings
from typing import Any, Optional, Tuple, cast

import numpy as np
from scipy import sparse


def _column_mean_variance(
    matrix: Any,
    max_memory: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:

    n_obs = int(matrix.shape[0])
    n_vars = int(matrix.shape[1])
    if sparse.issparse(matrix):
        sparse_format = getattr(matrix, "format", None)
        if sparse_format == "csr":
            sums = np.bincount(
                matrix.indices,
                weights=matrix.data,
                minlength=n_vars,
            )
            squared_sums = np.bincount(
                matrix.indices,
                weights=matrix.data * matrix.data,
                minlength=n_vars,
            )
        elif sparse_format == "csc":
            sums = np.zeros(n_vars, dtype=np.float64)
            squared_sums = np.zeros(n_vars, dtype=np.float64)
            nonempty_columns = np.diff(matrix.indptr) > 0
            if np.any(nonempty_columns):
                starts = matrix.indptr[:-1][nonempty_columns]
                data = np.asarray(matrix.data, dtype=np.float64)
                sums[nonempty_columns] = np.add.reduceat(data, starts)
                squared_sums[nonempty_columns] = np.add.reduceat(data * data, starts)
        else:
            return _column_mean_variance(
                matrix.tocsr(),
                max_memory=max_memory,
            )

        means = sums / n_obs
        mean_squares = squared_sums / n_obs
    else:
        dense_matrix = cast(np.ndarray, matrix)
        means = np.empty(n_vars, dtype=np.float64)
        mean_squares = np.empty(n_vars, dtype=np.float64)
        chunk_size = _dense_column_chunk_size(dense_matrix, max_memory)
        for start in range(0, n_vars, chunk_size):
            end = min(start + chunk_size, n_vars)
            dense_chunk = dense_matrix[:, start:end]
            means[start:end] = np.asarray(
                dense_chunk.mean(axis=0, dtype=np.float64)
            ).ravel()
            squared_chunk = np.multiply(dense_chunk, dense_chunk)
            mean_squares[start:end] = np.asarray(
                squared_chunk.mean(axis=0, dtype=np.float64)
            ).ravel()

    variances = np.maximum(mean_squares - means**2, 0)
    if n_obs > 1:
        variances *= n_obs / (n_obs - 1)
    return means, variances


def _dense_column_chunk_size(
    matrix: np.ndarray,
    max_memory: Optional[int],
) -> int:

    n_obs = int(matrix.shape[0])
    n_vars = int(matrix.shape[1])
    if max_memory is None:
        return max(1, n_vars)

    bytes_per_value = np.dtype(matrix.dtype).itemsize
    bytes_per_column = n_obs * bytes_per_value
    if max_memory < bytes_per_column:
        warnings.warn(
            "Requested memory budget is smaller than the memory required for a "
            "single variable. Computation will proceed with chunk_size=1.",
            RuntimeWarning,
            stacklevel=2,
        )
    return max(1, int(max_memory // bytes_per_column))
