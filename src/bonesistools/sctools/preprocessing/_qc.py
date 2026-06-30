#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Collection as CollectionInstance
from typing import Any, Optional, Sequence, Tuple, Union, cast, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ..._compat import Literal
from ..._validation import _as_boolean, _as_positive_integer, _as_string
from .._typing import anndata_checker


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    copy: Literal[False] = False,
) -> None: ...


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    *,
    copy: Literal[True],
) -> AnnData: ...


@overload
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def qc(
    adata: AnnData,
    expression: Optional[str] = None,
    qc_vars: Union[str, Sequence[str]] = (),
    percent_top: Optional[Sequence[int]] = (50, 100, 200, 500),
    log1p: bool = True,
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
    """

    expression = None if expression is None else _as_string(expression, "expression")
    qc_vars = _as_qc_vars(qc_vars)
    percent_top = _as_percent_top(percent_top)
    log1p = _as_boolean(log1p, "log1p")
    copy = _as_boolean(copy, "copy")

    adata = adata.copy() if copy else adata
    expression_mtx = adata.X if expression is None else adata.layers[expression]
    obs_metrics, var_metrics = _calculate_qc_metrics(
        expression_mtx,
        cast(pd.DataFrame, adata.var),
        qc_vars=qc_vars,
        percent_top=percent_top,
        log1p=log1p,
    )
    obs_metrics.index = adata.obs_names
    var_metrics.index = adata.var_names

    for column in obs_metrics:
        adata.obs[column] = obs_metrics[column].to_numpy()
    for column in var_metrics:
        adata.var[column] = var_metrics[column].to_numpy()

    return adata if copy else None


def _calculate_qc_metrics(
    expression_mtx: Any,
    var: pd.DataFrame,
    *,
    qc_vars: Tuple[str, ...],
    percent_top: Optional[Tuple[int, ...]],
    log1p: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    matrix = (
        expression_mtx.tocsr() if sparse.issparse(expression_mtx) else expression_mtx
    )
    n_obs, n_vars = matrix.shape
    total_per_obs = _sum_axis(matrix, axis=1)
    total_per_var = _sum_axis(matrix, axis=0)
    variance_per_var = _variance_axis0(matrix)
    median_per_var, mad_per_var = _median_mad_axis0(matrix)
    n_vars_per_obs, n_cells_per_var = _positive_counts(matrix)

    obs_metrics = pd.DataFrame()
    obs_metrics["n_features"] = n_vars_per_obs
    if log1p:
        obs_metrics["log1p_n_features"] = np.log1p(n_vars_per_obs)
    obs_metrics["total"] = total_per_obs
    if log1p:
        obs_metrics["log1p_total"] = np.log1p(total_per_obs)

    if percent_top is not None:
        _check_percent_top(percent_top, n_vars)
        for top in percent_top:
            obs_metrics[f"pct_top{top}_features"] = _percent_top(
                matrix, top, total_per_obs
            )

    for qc_var in qc_vars:
        qc_mask = _qc_var_mask(var, qc_var)
        total_qc = _sum_axis(matrix[:, qc_mask], axis=1)
        obs_metrics[f"total_{qc_var}"] = total_qc
        if log1p:
            obs_metrics[f"log1p_total_{qc_var}"] = np.log1p(total_qc)
        obs_metrics[f"pct_{qc_var}"] = _percentage(total_qc, total_per_obs)

    var_metrics = pd.DataFrame()
    var_metrics["n_barcodes"] = n_cells_per_var
    var_metrics["mean"] = total_per_var / n_obs
    if log1p:
        var_metrics["log1p_mean"] = np.log1p(var_metrics["mean"].to_numpy())
    var_metrics["variance"] = variance_per_var
    if log1p:
        var_metrics["log1p_variance"] = np.log1p(variance_per_var)
    var_metrics["median"] = median_per_var
    if log1p:
        var_metrics["log1p_median"] = np.log1p(median_per_var)
    var_metrics["mad"] = mad_per_var
    if log1p:
        var_metrics["log1p_mad"] = np.log1p(mad_per_var)
    var_metrics["pct_dropout"] = 100 * (1 - n_cells_per_var / n_obs)
    var_metrics["total"] = total_per_var
    if log1p:
        var_metrics["log1p_total"] = np.log1p(total_per_var)

    return obs_metrics, var_metrics


def _as_qc_vars(qc_vars: Union[str, Sequence[str]]) -> Tuple[str, ...]:

    if isinstance(qc_vars, str):
        return (qc_vars,)
    if not isinstance(qc_vars, CollectionInstance):
        raise TypeError(
            f"unsupported argument type for 'qc_vars': "
            f"expected {str} or a sequence of {str} but received {type(qc_vars)}"
        )

    resolved = tuple(qc_vars)
    for qc_var in resolved:
        if not isinstance(qc_var, str):
            raise TypeError(
                f"unsupported element type in 'qc_vars': "
                f"expected {str} but received {type(qc_var)}"
            )

    return resolved


def _as_percent_top(
    percent_top: Optional[Sequence[int]],
) -> Optional[Tuple[int, ...]]:

    if percent_top is None:
        return None
    if isinstance(percent_top, str) or not isinstance(
        percent_top,
        CollectionInstance,
    ):
        raise TypeError(
            f"unsupported argument type for 'percent_top': "
            f"expected a sequence of {int} or None but received {type(percent_top)}"
        )

    return tuple(_as_positive_integer(value, "percent_top") for value in percent_top)


def _check_percent_top(percent_top: Tuple[int, ...], n_vars: int) -> None:

    if any(top > n_vars for top in percent_top):
        raise IndexError("Positions outside range of features.")


def _sum_axis(matrix: Any, axis: int) -> np.ndarray:

    return np.asarray(matrix.sum(axis=axis)).ravel()


def _variance_axis0(matrix: Any) -> np.ndarray:

    n_obs = matrix.shape[0]
    if n_obs <= 1:
        return np.full(matrix.shape[1], np.nan)

    if sparse.issparse(matrix):
        total = _sum_axis(matrix, axis=0)
        squared_total = _sum_axis(matrix.multiply(matrix), axis=0)
    else:
        array = np.asarray(matrix)
        total = array.sum(axis=0)
        squared_total = np.square(array).sum(axis=0)

    return (squared_total - total * total / n_obs) / (n_obs - 1)


def _median_mad_axis0(matrix: Any) -> Tuple[np.ndarray, np.ndarray]:

    if sparse.issparse(matrix):
        return _sparse_median_mad_axis0(matrix)

    array = np.asarray(matrix)
    medians = np.median(array, axis=0)
    mads = np.median(np.abs(array - medians), axis=0)
    return medians, mads


def _sparse_median_mad_axis0(matrix: Any) -> Tuple[np.ndarray, np.ndarray]:

    csc_matrix = matrix.tocsc()
    n_obs, n_vars = csc_matrix.shape
    medians = np.empty(n_vars, dtype=np.float64)
    mads = np.empty(n_vars, dtype=np.float64)

    for column in range(n_vars):
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
    total_per_obs: np.ndarray,
) -> np.ndarray:

    if sparse.issparse(matrix):
        top_sums = _sparse_top_sums(matrix, top)
    else:
        array = np.asarray(matrix)
        partitioned = np.partition(array, array.shape[1] - top, axis=1)
        top_sums = partitioned[:, -top:].sum(axis=1)

    return _percentage(top_sums, total_per_obs)


def _sparse_top_sums(matrix: Any, top: int) -> np.ndarray:

    matrix = matrix.tocsr()
    top_sums = np.zeros(matrix.shape[0], dtype=matrix.dtype)
    for row in range(matrix.shape[0]):
        start = matrix.indptr[row]
        stop = matrix.indptr[row + 1]
        data = matrix.data[start:stop]
        if data.size <= top:
            top_sums[row] = data.sum()
        else:
            top_sums[row] = np.partition(data, data.size - top)[-top:].sum()

    return top_sums


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
