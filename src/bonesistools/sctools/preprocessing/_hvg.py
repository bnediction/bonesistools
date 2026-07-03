#!/usr/bin/env python

from __future__ import annotations

import importlib
import warnings
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
from .._stats import _column_mean_variance
from .._typing import anndata_checker
from ..tools._utils import get_expression

HVGMethod = Literal["loess"]
HVG_METHODS: Tuple[HVGMethod, ...] = ("loess",)


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

    means, variances = _column_mean_variance(expression_mtx)
    expected_variances = _loess_fit(means, variances, span=span)
    scores = _normalized_variance_score(
        expression_mtx,
        means,
        expected_variances,
    )

    selected = np.zeros(adata.n_vars, dtype=bool)
    ranks = np.full(adata.n_vars, np.nan, dtype=float)

    eligible = np.flatnonzero(np.isfinite(scores))
    n_eligible = int(eligible.size)
    if n_eligible < n_features:
        warnings.warn(
            f"Requested {n_features} highly variable features, but only "
            f"{n_eligible} features have finite variability scores.",
            RuntimeWarning,
            stacklevel=2,
        )

    ordered_indices = eligible[np.argsort(-scores[eligible], kind="mergesort")]
    selected_indices = ordered_indices[:n_features]
    selected[selected_indices] = True
    ranks[selected_indices] = np.arange(1, selected_indices.size + 1, dtype=float)

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
    means: np.ndarray,
    variances: np.ndarray,
    span: float = 0.3,
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


def _normalized_variance_score(
    matrix: Any,
    means: np.ndarray,
    expected_variances: np.ndarray,
) -> np.ndarray:

    n_obs = int(matrix.shape[0])
    expected_variances = np.maximum(
        np.asarray(expected_variances, dtype=np.float64),
        np.finfo(np.float64).eps,
    )
    regularized_std = np.sqrt(expected_variances)
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

    denominator = (n_obs - 1) * np.square(regularized_std)
    numerator = n_obs * np.square(means) + squared_sum - 2 * clipped_sum * means
    scores = numerator / denominator
    scores[denominator == 0] = np.nan
    return scores
