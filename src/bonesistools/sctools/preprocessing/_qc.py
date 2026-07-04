#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Collection as CollectionInstance
from typing import Any, Dict, Optional, Sequence, Tuple, Union, cast, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ..._compat import Literal
from ..._validation import _as_boolean, _as_positive_integer
from .._typing import Matrix, anndata_checker
from ..tools._utils import get_expression

_numba_njit: Any
_numba_prange: Any
try:
    from numba import njit as _numba_njit
    from numba import prange as _numba_prange
except ImportError:  # pragma: no cover - optional acceleration dependency
    _numba_njit = None
    _numba_prange = range


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    *,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    backend: Literal["auto", "python", "numba"] = "auto",
    copy: Literal[False] = False,
) -> None: ...


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    *,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    backend: Literal["auto", "python", "numba"] = "auto",
    copy: Literal[True],
) -> AnnData: ...


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    *,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    backend: Literal["auto", "python", "numba"] = "auto",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    *,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    backend: Literal["auto", "python", "numba"] = "auto",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Calculate quality-control metrics for observations and variables.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    expression: str, optional
        Expression matrix used for quality-control metrics. If None, use
        `adata.X`; otherwise, use `adata.layers[expression]`.
    qc_vars: str or sequence of str (default: ())
        Boolean columns in `adata.var` identifying variable groups to summarize,
        such as mitochondrial genes.
    percent_top: sequence of int, optional (default: (50, 100, 200, 500))
        Ranks at which cumulative expression among the most expressed variables
        is reported. If None, top-rank percentages are not computed.
    log1p: bool (default: True)
        Whether to add log1p-transformed metric columns.
    backend: {"auto", "python", "numba"} (default: "auto")
        Backend used for sparse top-feature percentages and sparse median/MAD
        metrics. If "auto", use the Numba backend when available and otherwise
        use the reference Python backend. If "python", force the reference
        implementation. If "numba", require the accelerated Numba backend.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with quality-control metrics
        added. Otherwise, updates `adata` in place and returns None.

        Metric names follow the project naming convention: `n_` indicates a
        cardinal, `total` and `total_` indicate sums, `mean` indicates a mean,
        `variance` indicates a variance, `pct_` indicates a percentage, and
        `log1p_` indicates a log1p transformation.

        Observation-level metrics are stored in `adata.obs`:

        - `adata.obs["n_features"]`: number of detected features;
        - `adata.obs["total"]`: sum of matrix values per observation;
        - `adata.obs[f"pct_top{n}_features"]`: percentage of total expression
          contained in the top `n` most expressed features, for each value in
          `percent_top`;
        - `adata.obs[f"total_{qc_var}"]`: sum of matrix values assigned to
          each variable group listed in `qc_vars`;
        - `adata.obs[f"pct_{qc_var}"]`: percentage of expression assigned to each
          variable group listed in `qc_vars`;

        Variable-level metrics are stored in `adata.var`:

        - `adata.var["n_barcodes"]`: number of observations where the feature
          is detected;
        - `adata.var["mean"]`: mean expression per feature;
        - `adata.var["median"]`: median expression per feature;
        - `adata.var["variance"]`: sample variance per feature, normalized by
          `n_barcodes - 1`;
        - `adata.var["mad"]`: raw median absolute deviation per feature,
          without normal-consistency scaling;
        - `adata.var["pct_dropout"]`: percentage of observations where the
          feature is not detected;
        - `adata.var["total"]`: sum of matrix values per feature.

        If `log1p=True`, matching `log1p_` transformation columns are added,
        such as `log1p_n_features`, `log1p_total`, `log1p_mean`,
        `log1p_median`, `log1p_variance` and `log1p_mad`.

    Notes
    -----
    The `numba` backend uses JIT-compiled implementations and may produce small
    floating-point differences relative to the reference backend. Use
    `backend="python"` for reproducible debugging, and `backend="numba"` to
    require the accelerated implementation.
    """

    if isinstance(qc_vars, str):
        resolved_qc_vars = (qc_vars,)
    elif not isinstance(qc_vars, CollectionInstance):
        raise TypeError(
            f"unsupported argument type for 'qc_vars': "
            f"expected {str} or a sequence of {str} but received {type(qc_vars)}"
        )
    else:
        resolved_qc_vars = tuple(qc_vars)
        for qc_var in resolved_qc_vars:
            if not isinstance(qc_var, str):
                raise TypeError(
                    f"unsupported element type in 'qc_vars': "
                    f"expected {str} but received {type(qc_var)}"
                )

    if percent_top is None:
        resolved_percent_top = None
    elif isinstance(percent_top, str) or not isinstance(
        percent_top,
        CollectionInstance,
    ):
        raise TypeError(
            f"unsupported argument type for 'percent_top': "
            f"expected a sequence of {int} or None but received {type(percent_top)}"
        )
    else:
        resolved_percent_top = tuple(
            _as_positive_integer(value, "percent_top") for value in percent_top
        )

    log1p = _as_boolean(log1p, "log1p")
    copy = _as_boolean(copy, "copy")
    if backend not in {"auto", "python", "numba"}:
        raise ValueError(
            f"unsupported backend: {backend!r}; " "expected 'auto', 'python' or 'numba'"
        )
    if backend == "numba" and _numba_njit is None:
        raise ImportError(
            "backend='numba' requires the optional dependency 'numba'. "
            "Install it with:\n\n"
            "pip install numba\n\n"
            "or use backend='python'."
        )

    adata = adata.copy() if copy else adata
    expression_mtx = get_expression(adata, layer=expression, copy=False)
    matrix = (
        cast(Matrix, cast(Any, expression_mtx).tocsr())
        if sparse.issparse(expression_mtx)
        else expression_mtx
    )

    n_barcodes, n_features = matrix.shape
    total_per_barcode = _sum_axis(matrix, axis=1)
    total_per_feature = _sum_axis(matrix, axis=0)
    variance_per_feature = _variance_per_feature(matrix)
    median_per_feature, mad_per_feature = _median_mad_per_feature(
        matrix,
        backend=backend,
    )
    n_features_per_barcode, n_barcodes_per_feature = _positive_counts(matrix)

    obs_metrics = pd.DataFrame(index=adata.obs_names)
    obs_metrics["n_features"] = n_features_per_barcode
    if log1p:
        obs_metrics["log1p_n_features"] = np.log1p(n_features_per_barcode)
    obs_metrics["total"] = total_per_barcode
    if log1p:
        obs_metrics["log1p_total"] = np.log1p(total_per_barcode)

    if resolved_percent_top is not None:
        if any(top > n_features for top in resolved_percent_top):
            raise IndexError("Positions outside range of features.")
        top_percentages = _percent_top_by_rank(
            matrix,
            resolved_percent_top,
            total_per_barcode,
            backend=backend,
        )
        for top, values in top_percentages.items():
            obs_metrics[f"pct_top{top}_features"] = values

    var = cast(pd.DataFrame, adata.var)
    for qc_var in resolved_qc_vars:
        qc_mask = _qc_var_mask(var, qc_var)
        total_qc = _sum_axis(cast(Any, matrix)[:, qc_mask], axis=1)
        obs_metrics[f"total_{qc_var}"] = total_qc
        if log1p:
            obs_metrics[f"log1p_total_{qc_var}"] = np.log1p(total_qc)
        obs_metrics[f"pct_{qc_var}"] = _percentage(total_qc, total_per_barcode)

    var_metrics = pd.DataFrame(index=adata.var_names)
    var_metrics["n_barcodes"] = n_barcodes_per_feature
    var_metrics["mean"] = total_per_feature / n_barcodes
    var_metrics["variance"] = variance_per_feature
    var_metrics["median"] = median_per_feature
    var_metrics["mad"] = mad_per_feature
    var_metrics["pct_dropout"] = 100 * (1 - n_barcodes_per_feature / n_barcodes)
    var_metrics["total"] = total_per_feature

    if log1p:
        for column in ("mean", "variance", "median", "mad", "total"):
            var_metrics[f"log1p_{column}"] = np.log1p(
                var_metrics[column].to_numpy()
            )

    obs_metrics.index = adata.obs_names
    var_metrics.index = adata.var_names

    for column in obs_metrics:
        adata.obs[column] = obs_metrics[column].to_numpy()
    for column in var_metrics:
        adata.var[column] = var_metrics[column].to_numpy()

    return adata if copy else None


def _sum_axis(matrix: Any, axis: int) -> np.ndarray:

    return np.asarray(matrix.sum(axis=axis)).ravel()


def _variance_per_feature(matrix: Any) -> np.ndarray:

    n_obs = matrix.shape[0]
    if n_obs <= 1:
        return np.full(matrix.shape[1], np.nan)

    if sparse.issparse(matrix):
        total = _sum_axis(matrix, axis=0)
        squared_total = _sum_axis(matrix.multiply(matrix), axis=0)
    else:
        array = np.asarray(matrix)
        total = array.sum(axis=0)
        squared_total = np.einsum("ij,ij->j", array, array)

    return (squared_total - total * total / n_obs) / (n_obs - 1)


def _median_mad_per_feature(
    matrix: Any,
    backend: Literal["auto", "python", "numba"] = "auto",
) -> Tuple[np.ndarray, np.ndarray]:

    if sparse.issparse(matrix):
        csc_matrix = cast(Any, matrix).tocsc()
        if backend in {"auto", "numba"} and _median_mad_sparse_csc_numba is not None:
            return _median_mad_sparse_csc_numba(
                csc_matrix.data.astype(np.float64, copy=False),
                csc_matrix.indptr.astype(np.int64, copy=False),
                csc_matrix.shape[0],
            )

        n_obs, n_features = csc_matrix.shape
        medians = np.empty(n_features, dtype=np.float64)
        mads = np.empty(n_features, dtype=np.float64)

        for column in range(n_features):
            start = csc_matrix.indptr[column]
            stop = csc_matrix.indptr[column + 1]
            data = np.asarray(csc_matrix.data[start:stop])
            implicit_zeros = n_obs - data.size
            median = _median_with_implicit_value(data, n_obs, implicit_zeros, 0.0)
            deviations = np.abs(data - median)

            medians[column] = median
            mads[column] = _median_with_implicit_value(
                deviations,
                n_obs,
                implicit_zeros,
                abs(median),
            )

        return medians, mads

    array = np.asarray(matrix)
    medians = np.median(array, axis=0)
    mads = np.median(np.abs(array - medians), axis=0)
    return medians, mads


def _median_with_implicit_value(
    values: np.ndarray,
    size: int,
    implicit_count: int,
    implicit_value: float,
) -> float:

    if size == 0:
        return np.nan

    lower = (size - 1) // 2
    upper = size // 2
    return (
        _kth_with_implicit_value(values, lower, implicit_count, implicit_value)
        + _kth_with_implicit_value(values, upper, implicit_count, implicit_value)
    ) / 2


def _kth_with_implicit_value(
    values: np.ndarray,
    k: int,
    implicit_count: int,
    implicit_value: float,
) -> float:

    lower_values = values[values < implicit_value]
    equal_count = np.count_nonzero(values == implicit_value) + implicit_count
    if k < lower_values.size:
        return float(np.partition(lower_values, k)[k])
    if k < lower_values.size + equal_count:
        return float(implicit_value)

    upper_values = values[values > implicit_value]
    upper_k = k - lower_values.size - equal_count
    return float(np.partition(upper_values, upper_k)[upper_k])


if _numba_njit is None:  # pragma: no cover - import-time optional dependency branch
    _median_mad_sparse_csc_numba = None
else:  # pragma: no cover - Numba-compiled code is exercised but not line-traced

    @_numba_njit(cache=True)
    def _kth_with_implicit_sorted_numba(
        sorted_values: np.ndarray,
        k: int,
        implicit_count: int,
        implicit_value: float,
    ) -> float:

        lower_count = 0
        equal_count_data = 0
        for value in sorted_values:
            if value < implicit_value:
                lower_count += 1
            elif value == implicit_value:
                equal_count_data += 1

        equal_count = equal_count_data + implicit_count
        if k < lower_count:
            return float(sorted_values[k])
        if k < lower_count + equal_count:
            return float(implicit_value)

        upper_k = k - lower_count - equal_count
        upper_start = lower_count + equal_count_data
        return float(sorted_values[upper_start + upper_k])

    @_numba_njit(cache=True)
    def _median_with_implicit_value_numba(
        values: np.ndarray,
        size: int,
        implicit_count: int,
        implicit_value: float,
    ) -> float:

        if size == 0:
            return np.nan

        sorted_values = np.sort(values.copy())
        lower = (size - 1) // 2
        upper = size // 2
        return (
            _kth_with_implicit_sorted_numba(
                sorted_values,
                lower,
                implicit_count,
                implicit_value,
            )
            + _kth_with_implicit_sorted_numba(
                sorted_values,
                upper,
                implicit_count,
                implicit_value,
            )
        ) / 2

    @_numba_njit(cache=True, parallel=True)
    def _median_mad_sparse_csc_numba(
        data: np.ndarray,
        indptr: np.ndarray,
        n_obs: int,
    ) -> Tuple[np.ndarray, np.ndarray]:

        n_features = indptr.size - 1
        medians = np.empty(n_features, dtype=np.float64)
        mads = np.empty(n_features, dtype=np.float64)

        for column in _numba_prange(n_features):
            start = indptr[column]
            stop = indptr[column + 1]
            size = stop - start
            values = np.empty(size, dtype=np.float64)
            for index in range(size):
                values[index] = data[start + index]

            implicit_zeros = n_obs - size
            median = _median_with_implicit_value_numba(
                values,
                n_obs,
                implicit_zeros,
                0.0,
            )
            deviations = np.empty(size, dtype=np.float64)
            for index in range(size):
                deviations[index] = abs(values[index] - median)

            medians[column] = median
            mads[column] = _median_with_implicit_value_numba(
                deviations,
                n_obs,
                implicit_zeros,
                abs(median),
            )

        return medians, mads


def _positive_counts(matrix: Any) -> Tuple[np.ndarray, np.ndarray]:

    positive = matrix > 0
    if sparse.issparse(positive):
        return (
            np.asarray(positive.getnnz(axis=1)).ravel(),
            np.asarray(positive.getnnz(axis=0)).ravel(),
        )

    return (
        np.asarray(positive.sum(axis=1)).ravel(),
        np.asarray(positive.sum(axis=0)).ravel(),
    )


def _percent_top(
    matrix: Any,
    top: int,
    total_per_barcode: np.ndarray,
    backend: Literal["auto", "python", "numba"] = "auto",
) -> np.ndarray:

    return _percent_top_by_rank(matrix, (top,), total_per_barcode, backend=backend)[
        top
    ]


def _percent_top_by_rank(
    matrix: Any,
    tops: Sequence[int],
    total_per_barcode: np.ndarray,
    backend: Literal["auto", "python", "numba"] = "auto",
) -> Dict[int, np.ndarray]:

    sorted_tops = np.asarray(sorted(tops), dtype=np.int64)
    if not sparse.issparse(matrix):
        values = _percent_top_dense(np.asarray(matrix), sorted_tops, total_per_barcode)
        values_by_top = {
            int(top): values[:, index] * 100 for index, top in enumerate(sorted_tops)
        }
        return {top: values_by_top[top] for top in tops}

    csr_matrix = cast(Any, matrix).tocsr()
    if backend in {"auto", "numba"} and _percent_top_sparse_csr_numba is not None:
        values = _percent_top_sparse_csr_numba(
            csr_matrix.data,
            csr_matrix.indptr,
            sorted_tops,
        )
        values_by_top = {
            int(top): values[:, index] * 100 for index, top in enumerate(sorted_tops)
        }
        return {top: values_by_top[top] for top in tops}

    values = _percent_top_sparse_csr(
        csr_matrix.data,
        csr_matrix.indptr,
        sorted_tops,
    )
    values_by_top = {
        int(top): values[:, index] * 100 for index, top in enumerate(sorted_tops)
    }
    return {top: values_by_top[top] for top in tops}


def _percent_top_dense(
    matrix: np.ndarray,
    tops: np.ndarray,
    total_per_barcode: np.ndarray,
) -> np.ndarray:

    max_top = tops[-1]
    partitioned = np.partition(matrix, matrix.shape[1] - tops, axis=1)
    partitioned = partitioned[:, ::-1][:, :max_top]
    values = np.zeros((matrix.shape[0], tops.size))
    accumulated = np.zeros(matrix.shape[0])
    previous = 0

    for index, top in enumerate(tops):
        accumulated += partitioned[:, previous:top].sum(axis=1)
        values[:, index] = accumulated
        previous = top

    with np.errstate(divide="ignore", invalid="ignore"):
        return values / total_per_barcode[:, None]


def _percent_top_sparse_csr(
    data: np.ndarray,
    indptr: np.ndarray,
    tops: np.ndarray,
) -> np.ndarray:

    tops = np.sort(tops.astype(np.int64))
    max_top = tops[-1]
    n_obs = indptr.size - 1
    sums = np.zeros(n_obs, dtype=data.dtype)
    partitioned = np.zeros((n_obs, max_top), dtype=data.dtype)

    for row in range(n_obs):
        start = indptr[row]
        stop = indptr[row + 1]
        size = stop - start
        row_values = data[start:stop]
        sums[row] = row_values.sum()
        if size <= max_top:
            partitioned[row, :size] = row_values
        else:
            partitioned[row, :] = -np.partition(-row_values, max_top)[:max_top]
        partitioned[row, :] = np.partition(partitioned[row, :], max_top - tops)

    partitioned = partitioned[:, ::-1][:, :max_top]
    values = np.zeros((n_obs, tops.size))
    accumulated = np.zeros(n_obs, dtype=data.dtype)
    previous = 0

    for index, top in enumerate(tops):
        accumulated += partitioned[:, previous:top].sum(axis=1)
        values[:, index] = accumulated
        previous = top

    with np.errstate(divide="ignore", invalid="ignore"):
        return values / sums.reshape((n_obs, 1))


if _numba_njit is None:  # pragma: no cover - import-time optional dependency branch
    _percent_top_sparse_csr_numba = None
else:  # pragma: no cover - Numba-compiled code is exercised but not line-traced

    @_numba_njit(cache=True, parallel=True)
    def _percent_top_sparse_csr_numba(
        data: np.ndarray,
        indptr: np.ndarray,
        tops: np.ndarray,
    ) -> np.ndarray:

        tops = np.sort(tops.astype(np.int64))
        max_top = tops[-1]
        n_obs = indptr.size - 1
        sums = np.zeros(n_obs, dtype=data.dtype)
        partitioned = np.zeros((n_obs, max_top), dtype=data.dtype)
        values = np.empty((n_obs, tops.size), dtype=np.float64)

        for row in _numba_prange(n_obs):
            start = indptr[row]
            stop = indptr[row + 1]
            size = stop - start
            sums[row] = np.sum(data[start:stop])

            if size <= max_top:
                for index in range(size):
                    partitioned[row, index] = data[start + index]
            else:
                row_values = data[start:stop]
                top_values = -np.partition(-row_values, max_top)[:max_top]
                for index in range(max_top):
                    partitioned[row, index] = top_values[index]

            partitioned[row, :] = np.partition(
                partitioned[row, :],
                max_top - tops,
            )

        partitioned = partitioned[:, ::-1][:, :max_top]
        accumulated = np.zeros(n_obs, dtype=data.dtype)
        previous = 0

        for top_index in range(tops.size):
            top = tops[top_index]
            for row in _numba_prange(n_obs):
                segment_sum = np.sum(partitioned[row, previous:top])
                accumulated[row] += segment_sum
                values[row, top_index] = accumulated[row]
            previous = top

        return values / sums.reshape((n_obs, 1))


def _percentage(
    numerator: np.ndarray,
    denominator: np.ndarray,
) -> np.ndarray:

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.true_divide(numerator, denominator) * 100


def _qc_var_mask(var: pd.DataFrame, qc_var: str) -> np.ndarray:

    if qc_var not in var:
        raise KeyError(f"column {qc_var!r} not found in adata.var")

    values = var[qc_var]
    if not hasattr(values, "dtype") or values.dtype != bool:
        raise TypeError(
            f"unsupported column dtype for 'qc_vars': "
            f"expected boolean values in adata.var[{qc_var!r}]"
        )

    return values.to_numpy(dtype=bool)
