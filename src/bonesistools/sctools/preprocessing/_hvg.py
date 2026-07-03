#!/usr/bin/env python

from __future__ import annotations

import importlib
import warnings
from dataclasses import dataclass
from typing import Any, Optional, Tuple, cast, overload

import numpy as np
from anndata import AnnData
from scipy import sparse

from ..._compat import Literal
from ..._validation import (
    _as_boolean,
    _as_literal,
    _as_positive_integer,
    _as_positive_number,
    _as_string,
)
from .._typing import anndata_checker
from ..tools._utils import get_expression

HVGMethod = Literal["loess"]
HVG_METHODS: Tuple[HVGMethod, ...] = ("loess",)


@dataclass(frozen=True)
class _StandardizedExpression:

    n_obs: int
    means: np.ndarray
    regularized_std: np.ndarray
    clipped_sum: np.ndarray
    squared_sum: np.ndarray


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    n_features: int = 2000,
    method: HVGMethod = "loess",
    span: float = 0.3,
    key_added: str = "highly_variable",
    copy: Literal[False] = False,
) -> None: ...


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    n_features: int = 2000,
    method: HVGMethod = "loess",
    span: float = 0.3,
    key_added: str = "highly_variable",
    copy: Literal[True],
) -> AnnData: ...


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    n_features: int = 2000,
    method: HVGMethod = "loess",
    span: float = 0.3,
    key_added: str = "highly_variable",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    n_features: int = 2000,
    method: HVGMethod = "loess",
    span: float = 0.3,
    key_added: str = "highly_variable",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Select highly variable features.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    expression: str, optional
        Expression matrix used to estimate variability. If None, use
        `adata.X`; otherwise, use `adata.layers[expression]`.
    n_features: int (default: 2000)
        Maximum number of highly variable features to select.
    method: {'loess'} (default: 'loess')
        Method used to score highly variable features.

        The `loess` method reproduces the Seurat v3 highly variable
        feature selection strategy: LOESS-regularized mean-variance trend
        estimation, followed by clipping, standardization of expression
        values using the regularized standard deviation, and ranking by
        normalized variance scores [1].
    span: float (default: 0.3)
        LOESS smoothing span. Must be greater than 0 and smaller than or equal
        to 1.
    key_added: str (default: 'highly_variable')
        Column in `adata.var` where the selected feature mask is stored.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with highly variable feature
        results added. Otherwise, updates `adata` in place and returns None.

        HVG results are stored in:

        - `adata.var[key_added]`: selected highly variable features;
        - `adata.var[f"{key_added}_rank"]`: selected feature rank;
        - `adata.var[f"{key_added}_score"]`: normalized variance score;
        - `adata.uns[key_added]`: HVG metadata.

    References
    ----------
    Stuart et al. (2019). Comprehensive integration of single-cell data. Cell,
    177(7), 1888-1902.
    """

    expression = None if expression is None else _as_string(expression, "expression")
    n_features = _as_positive_integer(n_features, "n_features")
    method = cast(HVGMethod, _as_literal(method, choices=HVG_METHODS, name="method"))
    span = _as_positive_number(span, "span")
    if span > 1:
        raise ValueError(
            f"invalid argument value for 'span': "
            f"expected value smaller than or equal to 1 but received {span!r}"
        )
    key_added = _as_string(key_added, "key_added")
    copy = _as_boolean(copy, "copy")

    adata = adata.copy() if copy else adata
    expression_mtx = get_expression(
        adata,
        layer=expression,
        copy=False,
    )

    if len(expression_mtx.shape) != 2:
        raise ValueError("invalid expression matrix: expected a two-dimensional matrix")

    scores = _loess_fit(expression_mtx, span=span)
    ordered_indices = _order_by_score(scores)

    n_eligible = int(ordered_indices.size)
    if n_eligible < n_features:
        warnings.warn(
            f"Requested {n_features} highly variable features, but only "
            f"{n_eligible} features have finite variability scores.",
            RuntimeWarning,
            stacklevel=2,
        )

    selected = _select_features(
        ordered_indices,
        n_features=n_features,
        n_vars=adata.n_vars,
    )
    ranks = _assign_feature_ranks(
        ordered_indices,
        n_features=n_features,
        n_vars=adata.n_vars,
    )

    adata.var[key_added] = selected
    adata.var[f"{key_added}_rank"] = ranks
    adata.var[f"{key_added}_score"] = scores
    adata.uns[key_added] = {
        "method": method,
        "n_features": n_features,
        "n_selected": int(selected.sum()),
        "expression": expression,
        "params": {
            "fit": "loess",
            "span": span,
            "score": "normalized_variance_score",
        },
    }

    return adata if copy else None


def _loess_fit(
    matrix: Any,
    *,
    span: float,
) -> np.ndarray:

    mean = _compute_mean(matrix)
    variance = _compute_variance(matrix, mean)
    trend = _compute_loess_trend(mean, variance, span=span)
    std = _compute_regularized_std(trend)
    standardized = _standardize(matrix, mean, std)
    score = _compute_normalized_variance(standardized)

    return score


def _compute_mean(matrix: Any) -> np.ndarray:

    n_obs = int(matrix.shape[0])
    if sparse.issparse(matrix):
        sums = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64, copy=False)

        return sums / n_obs

    dense_matrix = cast(np.ndarray, matrix)

    return np.asarray(dense_matrix.mean(axis=0, dtype=np.float64)).ravel()


def _compute_variance(
    matrix: Any,
    mean: Optional[np.ndarray] = None,
) -> np.ndarray:

    if mean is None:
        mean = _compute_mean(matrix)

    n_obs = int(matrix.shape[0])
    n_vars = int(matrix.shape[1])
    if sparse.issparse(matrix):
        sparse_format = getattr(matrix, "format", None)
        if sparse_format == "csr":
            squared_sums = np.bincount(
                matrix.indices,
                weights=matrix.data * matrix.data,
                minlength=n_vars,
            )
        elif sparse_format == "csc":
            squared_sums = np.zeros(n_vars, dtype=np.float64)
            nonempty_columns = np.diff(matrix.indptr) > 0
            if np.any(nonempty_columns):
                starts = matrix.indptr[:-1][nonempty_columns]
                data = np.asarray(matrix.data, dtype=np.float64)
                squared_sums[nonempty_columns] = np.add.reduceat(
                    data * data,
                    starts,
                )
        else:
            return _compute_variance(
                matrix.tocsr(),
                mean=mean,
            )
        mean_squares = squared_sums / n_obs
    else:
        dense_matrix = cast(np.ndarray, matrix)
        squared_matrix = np.multiply(dense_matrix, dense_matrix)
        mean_squares = np.asarray(
            squared_matrix.mean(axis=0, dtype=np.float64)
        ).ravel()

    variances = np.maximum(mean_squares - mean**2, 0)
    if n_obs > 1:
        variances *= n_obs / (n_obs - 1)
    return variances


def _compute_dispersion(
    mean: np.ndarray,
    variance: np.ndarray,
) -> np.ndarray:

    with np.errstate(divide="ignore", invalid="ignore"):
        dispersion = variance / mean
    dispersion[~np.isfinite(dispersion)] = np.nan

    return dispersion


def _compute_loess_trend(
    means: np.ndarray,
    variances: np.ndarray,
    *,
    span: float,
) -> np.ndarray:

    try:
        loess_module = importlib.import_module("skmisc.loess")
    except ImportError as error:
        raise ImportError(
            "method='loess' requires the optional dependency "
            "'scikit-misc'. Install it with:\n\n"
            "pip install scikit-misc"
        ) from error
    loess = cast(Any, getattr(loess_module, "loess"))

    log_expected_variances = np.zeros(means.shape[0], dtype=np.float64)
    not_constant = (variances > 0) & (means > 0)
    if np.any(not_constant):
        x = np.log10(means[not_constant])
        y = np.log10(variances[not_constant])
        model = loess(x, y, span=span, degree=2)
        model.fit()
        log_expected_variances[not_constant] = model.outputs.fitted_values
    return cast(np.ndarray, 10**log_expected_variances)


def _compute_regularized_std(trend: np.ndarray) -> np.ndarray:

    trend = np.maximum(
        np.asarray(trend, dtype=np.float64),
        np.finfo(np.float64).eps,
    )

    return np.sqrt(trend)


def _standardize(
    matrix: Any,
    means: np.ndarray,
    regularized_std: np.ndarray,
) -> _StandardizedExpression:

    n_obs = int(matrix.shape[0])
    clip_values = means + regularized_std * np.sqrt(n_obs)

    if sparse.issparse(matrix):
        counts = sparse.csr_matrix(matrix.astype(np.float64, copy=True))
        mask = counts.data > clip_values[counts.indices]
        counts.data[mask] = clip_values[counts.indices[mask]]

        clipped_sum = np.asarray(counts.sum(axis=0)).ravel()
        np.square(counts.data, out=counts.data)
        squared_sum = np.asarray(counts.sum(axis=0)).ravel()
    else:
        counts = np.asarray(matrix, dtype=np.float64).copy()
        np.minimum(counts, clip_values.reshape(1, -1), out=counts)
        clipped_sum = counts.sum(axis=0)
        np.square(counts, out=counts)
        squared_sum = counts.sum(axis=0)

    return _StandardizedExpression(
        n_obs=n_obs,
        means=means,
        regularized_std=regularized_std,
        clipped_sum=clipped_sum,
        squared_sum=squared_sum,
    )


def _compute_normalized_variance(
    standardized: _StandardizedExpression,
) -> np.ndarray:

    n_obs = standardized.n_obs
    means = standardized.means
    regularized_std = standardized.regularized_std
    clipped_sum = standardized.clipped_sum
    squared_sum = standardized.squared_sum

    denominator = (n_obs - 1) * np.square(regularized_std)
    numerator = n_obs * np.square(means) + squared_sum - 2 * clipped_sum * means
    scores = numerator / denominator
    scores[denominator == 0] = np.nan
    return scores


def _order_by_score(scores: np.ndarray) -> np.ndarray:

    eligible = np.flatnonzero(np.isfinite(scores))

    return eligible[np.argsort(-scores[eligible], kind="mergesort")]


def _select_features(
    ordered_indices: np.ndarray,
    *,
    n_features: int,
    n_vars: int,
) -> np.ndarray:

    selected = np.zeros(n_vars, dtype=bool)
    selected_indices = ordered_indices[:n_features]
    selected[selected_indices] = True

    return selected


def _assign_feature_ranks(
    ordered_indices: np.ndarray,
    *,
    n_features: int,
    n_vars: int,
) -> np.ndarray:

    ranks = np.full(n_vars, np.nan, dtype=float)
    selected_indices = ordered_indices[:n_features]
    ranks[selected_indices] = np.arange(1, selected_indices.size + 1, dtype=float)

    return ranks
