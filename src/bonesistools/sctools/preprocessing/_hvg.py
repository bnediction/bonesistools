#!/usr/bin/env python

from __future__ import annotations

import importlib
import numbers
import warnings
from collections.abc import Sequence as SequenceInstance
from dataclasses import dataclass
from typing import Any, Optional, Tuple, Union, cast, overload

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
from .._typing import Matrix, anndata_checker
from ..tools._utils import get_expression

HVGMethod = Literal["loess", "binning"]
HVG_METHODS: Tuple[HVGMethod, ...] = ("loess", "binning")
HVG_BATCH_WEIGHTINGS: Tuple[Literal["equal", "cell_count"], ...] = (
    "equal",
    "cell_count",
)
HVG_BATCH_SELECTIONS: Tuple[Literal["consensus", "rank"], ...] = (
    "consensus",
    "rank",
)
_DEFAULT_SCORE_CUTOFFS = {
    "loess": (2.0, np.inf),
    "binning": (0.5, np.inf),
}
_NORMAL_CONSISTENCY_MAD_SCALE = 0.6744897501960817


@dataclass(frozen=True)
class _HVGScore:

    means: np.ndarray
    scores: np.ndarray


@dataclass(frozen=True)
class _BatchHVGScore:

    means: np.ndarray
    scores: np.ndarray
    nbatches: np.ndarray
    rank_summary: np.ndarray


@dataclass(frozen=True)
class _HVGSelection:

    selected: np.ndarray
    ranks: np.ndarray
    score_cutoff: Optional[Tuple[float, float]]
    mean_cutoff: Optional[Tuple[float, float]]


@dataclass(frozen=True)
class _BatchGroups:

    masks: Tuple[np.ndarray, ...]
    sizes: np.ndarray


@dataclass(frozen=True)
class _ClippedMoments:

    n_obs: int
    means: np.ndarray
    regularized_variance: np.ndarray
    clipped_sum: np.ndarray
    squared_sum: np.ndarray


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    method: Literal["loess", "binning"] = "loess",
    n_features: Optional[int] = 2000,
    mean: Tuple[float, float] = (0.0125, 3.0),
    score: Union[Literal["auto"], Tuple[float, float]] = "auto",
    span: float = 0.3,
    n_bins: int = 20,
    batch_key: Optional[str] = None,
    batch_selection: Literal["consensus", "rank"] = "consensus",
    batch_weighting: Literal["equal", "cell_count"] = "equal",
    key_added: str = "highly_variable",
    copy: Literal[False] = False,
) -> None: ...


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    method: Literal["loess", "binning"] = "loess",
    n_features: Optional[int] = 2000,
    mean: Tuple[float, float] = (0.0125, 3.0),
    score: Union[Literal["auto"], Tuple[float, float]] = "auto",
    span: float = 0.3,
    n_bins: int = 20,
    batch_key: Optional[str] = None,
    batch_selection: Literal["consensus", "rank"] = "consensus",
    batch_weighting: Literal["equal", "cell_count"] = "equal",
    key_added: str = "highly_variable",
    copy: Literal[True],
) -> AnnData: ...


@overload
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    method: Literal["loess", "binning"] = "loess",
    n_features: Optional[int] = 2000,
    mean: Tuple[float, float] = (0.0125, 3.0),
    score: Union[Literal["auto"], Tuple[float, float]] = "auto",
    span: float = 0.3,
    n_bins: int = 20,
    batch_key: Optional[str] = None,
    batch_selection: Literal["consensus", "rank"] = "consensus",
    batch_weighting: Literal["equal", "cell_count"] = "equal",
    key_added: str = "highly_variable",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def hvg(
    adata: AnnData,
    *,
    expression: Optional[str] = None,
    method: Literal["loess", "binning"] = "loess",
    n_features: Optional[int] = 2000,
    mean: Tuple[float, float] = (0.0125, 3.0),
    score: Union[Literal["auto"], Tuple[float, float]] = "auto",
    span: float = 0.3,
    n_bins: int = 20,
    batch_key: Optional[str] = None,
    batch_selection: Literal["consensus", "rank"] = "consensus",
    batch_weighting: Literal["equal", "cell_count"] = "equal",
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
        `adata.X`; otherwise, use `adata.layers[expression]`. Use raw counts
        for `method='loess'` and log-normalized expression values for
        `method='binning'`.
    method: {'loess', 'binning'} (default: 'loess')
        Method used to score highly variable features.

        The `loess` method reproduces the Seurat v3
        LOESS-regularized mean-variance strategy: trend estimation,
        followed by clipping, standardization of expression values
        using the regularized variance, and ranking by  normalized
        variance scores [1].

        The `binning` method reproduces the Cell Ranger
        mean-binned normalized dispersion strategy: features are
        grouped by mean expression, dispersions are normalized within
        each bin using the median and median absolute deviation, and
        features are ranked by normalized dispersion scores [2, 3].
    n_features: int or None (default: 2000)
        Number of highly variable features to select. If None, features are
        selected using `score` and `mean` cutoffs instead.
    mean: tuple of float (default: (0.0125, 3.0))
        Mean-expression interval used when `n_features=None`. Features are
        selected when `mean[0] < feature_mean < mean[1]`. For
        `method='loess'`, the interval is applied to the mean of log1p
        library-size-normalized counts with target sum 1e4. For
        `method='binning'`, the interval is applied to the mean of the
        provided expression matrix.
    score: {'auto'} or tuple of float (default: 'auto')
        Score interval used when `n_features=None`. Features are selected when
        `score[0] < feature_score < score[1]`. If 'auto', use
        method-specific defaults: `(2.0, np.inf)` for `method='loess'` and
        `(0.5, np.inf)` for `method='binning'`. For `method='loess'`, scores
        are normalized variance scores. For `method='binning'`, scores are
        normalized dispersion scores.
    span: float (default: 0.3)
        LOESS smoothing span. Must be greater than 0 and smaller than or equal
        to 1. Used by `method='loess'`.
    n_bins: int (default: 20)
        Number of mean-expression bins used by `method='binning'`.
    batch_key: str, optional
        Column in `adata.obs` defining batches. If provided, features are
        scored independently in each batch and scores are averaged across
        batches.
    batch_selection: {'consensus', 'rank'} (default: 'consensus')
        Batch-aware top-feature selection strategy used when both `batch_key`
        and `n_features` are provided. If 'consensus', prioritize features
        selected in the largest number of batches. If 'rank', prioritize
        features with the best summarized within-batch rank.
    batch_weighting: {'equal', 'cell_count'} (default: 'equal')
        Batch score aggregation strategy. If 'equal', each batch has the same
        weight. If 'cell_count', batches are weighted by their number of
        observations. When `n_features` is provided, features selected in more
        batches are prioritized before applying method-specific tie-breaks. If
        `n_features` is None, cutoffs are applied to the aggregated scores and
        means. Ignored when `batch_key` is None.
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
        - `adata.var[f"{key_added}_score"]`: method-specific variability score;
        - `adata.uns[key_added]`: HVG metadata.

    References
    ----------
    [1] Stuart et al. (2019). Comprehensive integration of single-cell data. Cell,
    177(7), 1888-1902.

    [2] Macosko et al. (2015). Highly parallel genome-wide expression
    profiling of individual cells using nanoliter droplets. Cell, 161(5),
    1202-1214.

    [3] Zheng et al. (2017). Massively parallel digital transcriptional
    profiling of single cells. Nature Communications, 8, 14049.
    """

    expression = None if expression is None else _as_string(expression, "expression")
    n_features = (
        None
        if n_features is None
        else _as_positive_integer(n_features, "n_features")
    )
    method = cast(HVGMethod, _as_literal(method, choices=HVG_METHODS, name="method"))
    batch_key = None if batch_key is None else _as_string(batch_key, "batch_key")
    batch_weighting = cast(
        Literal["equal", "cell_count"],
        _as_literal(
            batch_weighting,
            choices=HVG_BATCH_WEIGHTINGS,
            name="batch_weighting",
        ),
    )
    batch_selection = cast(
        Literal["consensus", "rank"],
        _as_literal(
            batch_selection,
            choices=HVG_BATCH_SELECTIONS,
            name="batch_selection",
        ),
    )
    span = _as_positive_number(span, "span")
    if span > 1:
        raise ValueError(
            f"invalid argument value for 'span': "
            f"expected value smaller than or equal to 1 but received {span!r}"
        )
    n_bins = _as_positive_integer(n_bins, "n_bins")
    key_added = _as_string(key_added, "key_added")
    score_cutoff = _as_optional_cutoff_interval(score, "score")
    mean_cutoff = _as_cutoff_interval(mean, "mean")
    copy = _as_boolean(copy, "copy")

    if (
        n_features is not None
        and (
            score_cutoff is not None
            or mean_cutoff != (0.0125, 3.0)
        )
    ):
        warnings.warn(
            "`score` and `mean` cutoffs are ignored when `n_features` is "
            "provided.",
            UserWarning,
            stacklevel=2,
        )

    adata = adata.copy() if copy else adata
    expression_mtx = get_expression(
        adata,
        layer=expression,
        copy=False,
    )

    if len(expression_mtx.shape) != 2:
        raise ValueError("invalid expression matrix: expected a two-dimensional matrix")

    batch_groups = _resolve_batch_groups(adata, batch_key)
    scoring = (
        _score_hvg_features(
            expression_mtx,
            method=method,
            span=span,
            n_bins=n_bins,
        )
        if batch_groups is None
        else _score_hvg_features_by_batch(
            expression_mtx,
            batch_groups=batch_groups,
            batch_weighting=batch_weighting,
            n_features=n_features,
            method=method,
            span=span,
            n_bins=n_bins,
        )
    )
    ordered_indices = (
        _order_by_score(scoring.scores)
        if batch_groups is None or n_features is None
        else _order_batch_hvg_features(
            cast(_BatchHVGScore, scoring),
            method=method,
            batch_selection=batch_selection,
        )
    )
    selection_means = (
        _selection_means(
            expression_mtx,
            method=method,
            scoring_means=scoring.means,
            batch_groups=batch_groups,
            batch_weighting=batch_weighting,
        )
        if n_features is None
        else scoring.means
    )
    selection = _select_hvg_features(
        ordered_indices=ordered_indices,
        scores=scoring.scores,
        means=selection_means,
        n_features=n_features,
        score=score_cutoff,
        mean=mean_cutoff,
        method=method,
        n_vars=adata.n_vars,
    )

    adata.var[key_added] = selection.selected
    adata.var[f"{key_added}_rank"] = selection.ranks
    adata.var[f"{key_added}_score"] = scoring.scores
    adata.uns[key_added] = {
        "method": method,
        "n_features": n_features,
        "n_selected": int(selection.selected.sum()),
        "expression": expression,
        "params": _hvg_params(
            method=method,
            span=span,
            n_bins=n_bins,
            n_features=n_features,
            selection=selection,
            batch_key=batch_key,
            batch_weighting=batch_weighting,
            batch_selection=batch_selection,
        ),
    }

    return adata if copy else None


def _score_hvg_features(
    matrix: Matrix,
    *,
    method: HVGMethod,
    span: float,
    n_bins: int,
) -> _HVGScore:

    if method == "loess":
        return _loess_fit(
            matrix,
            span=span,
        )
    if method == "binning":
        return _binning_fit(
            matrix,
            n_bins=n_bins,
        )

    raise ValueError(f"unsupported HVG method: {method!r}")


def _score_hvg_features_by_batch(
    matrix: Matrix,
    *,
    batch_groups: _BatchGroups,
    batch_weighting: Literal["equal", "cell_count"],
    n_features: Optional[int],
    method: HVGMethod,
    span: float,
    n_bins: int,
) -> _BatchHVGScore:

    batch_means = []
    batch_scores = []
    batch_ranks = []
    nbatches = np.zeros(matrix.shape[1], dtype=int)
    for mask in batch_groups.masks:
        batch_matrix = _subset_obs_matrix(matrix, mask)
        if int(mask.sum()) < 2:
            batch_means.append(_compute_mean(batch_matrix))
            batch_scores.append(np.full(matrix.shape[1], np.nan, dtype=np.float64))
            batch_ranks.append(np.full(matrix.shape[1], np.nan, dtype=np.float64))
            continue

        scoring = _score_hvg_features(
            batch_matrix,
            method=method,
            span=span,
            n_bins=n_bins,
        )
        batch_means.append(scoring.means)
        batch_scores.append(scoring.scores)
        batch_rank = np.full(matrix.shape[1], np.nan, dtype=np.float64)
        if n_features is not None:
            selected_indices = _order_by_score(scoring.scores)[:n_features]
            nbatches[selected_indices] += 1
            batch_rank[selected_indices] = np.arange(
                1,
                selected_indices.size + 1,
                dtype=np.float64,
            )
        batch_ranks.append(batch_rank)

    weights = _batch_weights(batch_groups, batch_weighting=batch_weighting)

    return _BatchHVGScore(
        means=_weighted_nanmean(np.vstack(batch_means), weights),
        scores=_weighted_nanmean(np.vstack(batch_scores), weights),
        nbatches=nbatches,
        rank_summary=_batch_rank_summary(
            np.vstack(batch_ranks),
            weights,
            batch_weighting=batch_weighting,
        ),
    )


def _hvg_params(
    *,
    method: HVGMethod,
    span: float,
    n_bins: int,
    n_features: Optional[int],
    selection: _HVGSelection,
    batch_key: Optional[str],
    batch_weighting: Literal["equal", "cell_count"],
    batch_selection: Literal["consensus", "rank"],
) -> dict:

    params = {
        "selection": "top_n" if n_features is not None else "cutoff",
        "score_cutoff": selection.score_cutoff,
        "mean_cutoff": selection.mean_cutoff,
        "batch_key": batch_key,
        "batch_weighting": batch_weighting,
        "batch_selection": batch_selection,
    }
    if method == "loess":
        return {
            "fit": "loess",
            "span": span,
            "score": "normalized_variance_score",
            **params,
        }

    return {
        "fit": "mean_binned_dispersion",
        "n_bins": n_bins,
        "score": "normalized_dispersion_score",
        **params,
    }


def _loess_fit(
    matrix: Matrix,
    *,
    span: float,
) -> _HVGScore:

    mean = _compute_mean(matrix)
    variance = _compute_variance(matrix, mean)
    trend = _compute_loess_trend(mean, variance, span=span)
    regularized_variance = _compute_regularized_variance(trend)
    clipped_moments = _compute_clipped_moments(
        matrix,
        mean,
        regularized_variance,
    )
    score = _compute_normalized_variance(clipped_moments)

    return _HVGScore(
        means=mean,
        scores=score,
    )


def _binning_fit(
    matrix: Matrix,
    *,
    n_bins: int,
) -> _HVGScore:

    mean = _compute_mean(matrix)
    variance = _compute_variance(matrix, mean)
    dispersion = _compute_dispersion(mean, variance)
    bins = _bin_by_mean(mean, n_bins=n_bins)
    score = _normalize_dispersion(dispersion, bins)

    return _HVGScore(
        means=mean,
        scores=score,
    )


def _compute_mean(matrix: Matrix) -> np.ndarray:

    n_obs = int(matrix.shape[0])
    if sparse.issparse(matrix):
        sparse_matrix = cast(Any, matrix)
        sums = (
            np.asarray(sparse_matrix.sum(axis=0))
            .ravel()
            .astype(np.float64, copy=False)
        )

        return sums / n_obs

    dense_matrix = cast(np.ndarray, matrix)

    return np.asarray(dense_matrix.mean(axis=0, dtype=np.float64)).ravel()


def _selection_means(
    matrix: Matrix,
    *,
    method: HVGMethod,
    scoring_means: np.ndarray,
    batch_groups: Optional[_BatchGroups],
    batch_weighting: Literal["equal", "cell_count"],
) -> np.ndarray:

    if method == "loess":
        if batch_groups is not None:
            return _compute_log_normalized_mean_by_batch(
                matrix,
                batch_groups=batch_groups,
                batch_weighting=batch_weighting,
            )
        return _compute_log_normalized_mean(matrix)
    return scoring_means


def _compute_log_normalized_mean(
    matrix: Matrix,
    *,
    target_sum: float = 1e4,
) -> np.ndarray:

    n_obs = int(matrix.shape[0])
    n_vars = int(matrix.shape[1])
    if sparse.issparse(matrix):
        sparse_matrix = cast(Any, matrix)
        counts = (
            sparse_matrix
            if sparse.isspmatrix_csr(sparse_matrix)
            else sparse.csr_matrix(sparse_matrix)
        )
        totals = np.asarray(counts.sum(axis=1)).ravel()
        scale = np.divide(
            target_sum,
            totals,
            out=np.zeros(totals.shape[0], dtype=np.float64),
            where=totals != 0,
        )
        values = np.asarray(counts.data, dtype=np.float64) * np.repeat(
            scale,
            np.diff(counts.indptr),
        )
        np.log1p(values, out=values)
        sums = np.bincount(
            counts.indices,
            weights=values,
            minlength=n_vars,
        )

        return sums / n_obs

    counts = np.asarray(matrix, dtype=np.float64)
    totals = counts.sum(axis=1)
    scale = np.divide(
        target_sum,
        totals,
        out=np.zeros(totals.shape[0], dtype=np.float64),
        where=totals != 0,
    )
    normalized = counts * scale.reshape(-1, 1)
    np.log1p(normalized, out=normalized)

    return np.asarray(normalized.mean(axis=0, dtype=np.float64)).ravel()


def _compute_log_normalized_mean_by_batch(
    matrix: Matrix,
    *,
    batch_groups: _BatchGroups,
    batch_weighting: Literal["equal", "cell_count"],
) -> np.ndarray:

    batch_means = [
        _compute_log_normalized_mean(_subset_obs_matrix(matrix, mask))
        for mask in batch_groups.masks
    ]
    weights = _batch_weights(batch_groups, batch_weighting=batch_weighting)

    return _weighted_nanmean(np.vstack(batch_means), weights)


def _compute_variance(
    matrix: Matrix,
    mean: Optional[np.ndarray] = None,
) -> np.ndarray:

    if mean is None:
        mean = _compute_mean(matrix)

    n_obs = int(matrix.shape[0])
    n_vars = int(matrix.shape[1])
    if sparse.issparse(matrix):
        sparse_matrix = cast(Any, matrix)
        sparse_format = getattr(sparse_matrix, "format", None)
        if sparse_format == "csr":
            squared_sums = np.bincount(
                sparse_matrix.indices,
                weights=sparse_matrix.data * sparse_matrix.data,
                minlength=n_vars,
            )
        elif sparse_format == "csc":
            squared_sums = np.zeros(n_vars, dtype=np.float64)
            nonempty_columns = np.diff(sparse_matrix.indptr) > 0
            if np.any(nonempty_columns):
                starts = sparse_matrix.indptr[:-1][nonempty_columns]
                data = np.asarray(sparse_matrix.data, dtype=np.float64)
                squared_sums[nonempty_columns] = np.add.reduceat(
                    data * data,
                    starts,
                )
        else:
            return _compute_variance(
                sparse_matrix.tocsr(),
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


def _bin_by_mean(
    means: np.ndarray,
    *,
    n_bins: int,
) -> np.ndarray:

    bins = np.full(means.shape[0], -1, dtype=int)
    finite_indices = np.flatnonzero(np.isfinite(means) & (means > 0))
    if finite_indices.size == 0:
        return bins

    finite_means = means[finite_indices]
    if n_bins == 20:
        percentiles = np.arange(10, 105, 5)
    else:
        percentiles = np.linspace(
            100 / n_bins,
            100,
            n_bins,
        )
    bin_edges = np.r_[-np.inf, np.percentile(finite_means, percentiles), np.inf]
    bin_edges = np.unique(bin_edges)
    bins[finite_indices] = np.searchsorted(
        bin_edges,
        finite_means,
        side="left",
    ) - 1

    return bins


def _normalize_dispersion(
    dispersion: np.ndarray,
    bins: np.ndarray,
) -> np.ndarray:

    scores = np.full(dispersion.shape[0], np.nan, dtype=np.float64)
    valid_bins = np.unique(bins[bins >= 0])
    for bin_index in valid_bins:
        bin_mask = bins == bin_index
        finite_mask = bin_mask & np.isfinite(dispersion)
        bin_dispersion = dispersion[finite_mask]
        if bin_dispersion.size == 0:
            continue

        median = np.median(bin_dispersion)
        mad = np.median(np.abs(bin_dispersion - median))
        mad = mad / _NORMAL_CONSISTENCY_MAD_SCALE

        with np.errstate(divide="ignore", invalid="ignore"):
            scores[finite_mask] = (bin_dispersion - median) / mad

    return scores


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


def _compute_regularized_variance(trend: np.ndarray) -> np.ndarray:

    return np.maximum(
        np.asarray(trend, dtype=np.float64),
        np.finfo(np.float64).eps,
    )


def _compute_clipped_moments(
    matrix: Matrix,
    means: np.ndarray,
    regularized_variance: np.ndarray,
) -> _ClippedMoments:

    n_obs = int(matrix.shape[0])
    clip_values = means + np.sqrt(regularized_variance * n_obs)

    if sparse.issparse(matrix):
        sparse_matrix = cast(Any, matrix)
        counts = sparse.csr_matrix(sparse_matrix.astype(np.float64, copy=True))
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

    return _ClippedMoments(
        n_obs=n_obs,
        means=means,
        regularized_variance=regularized_variance,
        clipped_sum=clipped_sum,
        squared_sum=squared_sum,
    )


def _compute_normalized_variance(
    clipped_moments: _ClippedMoments,
) -> np.ndarray:

    n_obs = clipped_moments.n_obs
    means = clipped_moments.means
    regularized_variance = clipped_moments.regularized_variance
    clipped_sum = clipped_moments.clipped_sum
    squared_sum = clipped_moments.squared_sum

    empirical_clipped_variance = (
        n_obs * np.square(means) + squared_sum - 2 * clipped_sum * means
    ) / (n_obs - 1)
    scores = empirical_clipped_variance / regularized_variance
    scores[regularized_variance == 0] = np.nan
    return scores


def _order_by_score(scores: np.ndarray) -> np.ndarray:

    eligible = np.flatnonzero(~np.isnan(scores))

    return eligible[np.argsort(-scores[eligible], kind="mergesort")]


def _order_batch_hvg_features(
    scoring: _BatchHVGScore,
    *,
    method: HVGMethod,
    batch_selection: Literal["consensus", "rank"],
) -> np.ndarray:

    eligible = np.flatnonzero(np.isfinite(scoring.scores))
    ranks = np.where(
        np.isnan(scoring.rank_summary),
        np.inf,
        scoring.rank_summary,
    )
    if batch_selection == "rank":
        return eligible[
            np.lexsort(
                (
                    -scoring.scores[eligible],
                    -scoring.nbatches[eligible],
                    ranks[eligible],
                )
            )
        ]

    if batch_selection != "consensus":
        raise ValueError(f"unsupported batch selection: {batch_selection!r}")

    if method == "loess":
        return eligible[
            np.lexsort(
                (
                    -scoring.scores[eligible],
                    ranks[eligible],
                    -scoring.nbatches[eligible],
                )
            )
        ]

    if method == "binning":
        return eligible[
            np.lexsort(
                (
                    -scoring.scores[eligible],
                    -scoring.nbatches[eligible],
                )
            )
        ]

    raise ValueError(f"unsupported HVG method: {method!r}")


def _resolve_batch_groups(
    adata: AnnData,
    batch_key: Optional[str],
) -> Optional[_BatchGroups]:

    if batch_key is None:
        return None
    if batch_key not in adata.obs:
        raise KeyError(f"column {batch_key!r} not found in adata.obs")

    labels = cast(Any, adata.obs[batch_key])
    valid_mask = ~np.asarray(labels.isna(), dtype=bool)
    values = list(labels[valid_mask].unique())
    masks = tuple(
        np.asarray((labels == value) & valid_mask, dtype=bool)
        for value in values
    )
    masks = tuple(mask for mask in masks if np.any(mask))
    if not masks:
        raise ValueError(
            f"invalid observation batches in adata.obs[{batch_key!r}]: "
            "expected at least one non-missing batch"
        )

    return _BatchGroups(
        masks=masks,
        sizes=np.asarray([mask.sum() for mask in masks], dtype=np.float64),
    )


def _subset_obs_matrix(matrix: Matrix, obs_mask: np.ndarray) -> Matrix:

    return cast(Matrix, cast(Any, matrix)[obs_mask, :])


def _batch_weights(
    batch_groups: _BatchGroups,
    *,
    batch_weighting: Literal["equal", "cell_count"],
) -> np.ndarray:

    if batch_weighting == "equal":
        return np.ones(batch_groups.sizes.shape[0], dtype=np.float64)
    if batch_weighting == "cell_count":
        return batch_groups.sizes.astype(np.float64, copy=False)

    raise ValueError(f"unsupported batch weighting: {batch_weighting!r}")


def _weighted_nanmean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:

    valid = ~np.isnan(values)
    weighted_values = np.where(valid, values, 0.0) * weights.reshape(-1, 1)
    denominators = np.sum(valid * weights.reshape(-1, 1), axis=0)

    return np.divide(
        weighted_values.sum(axis=0),
        denominators,
        out=np.full(values.shape[1], np.nan, dtype=np.float64),
        where=denominators != 0,
    )


def _batch_rank_summary(
    ranks: np.ndarray,
    weights: np.ndarray,
    *,
    batch_weighting: Literal["equal", "cell_count"],
) -> np.ndarray:

    if batch_weighting == "equal":
        return cast(
            np.ndarray,
            np.ma.median(np.ma.masked_invalid(ranks), axis=0).filled(np.nan),
        )

    return _weighted_nanmedian(ranks, weights)


def _weighted_nanmedian(values: np.ndarray, weights: np.ndarray) -> np.ndarray:

    medians = np.full(values.shape[1], np.nan, dtype=np.float64)
    for index in range(values.shape[1]):
        valid = ~np.isnan(values[:, index])
        if not np.any(valid):
            continue

        valid_values = values[valid, index]
        valid_weights = weights[valid]
        order = np.argsort(valid_values, kind="mergesort")
        sorted_values = valid_values[order]
        cumulative_weights = np.cumsum(valid_weights[order])
        cutoff = 0.5 * cumulative_weights[-1]
        median_index = np.searchsorted(cumulative_weights, cutoff, side="left")
        medians[index] = sorted_values[median_index]

    return medians


def _as_cutoff_interval(
    value: Tuple[float, float],
    name: str,
) -> Tuple[float, float]:

    if isinstance(value, (str, bytes)) or not isinstance(value, SequenceInstance):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected tuple of two floats but received {type(value)}"
        )
    if len(value) != 2:
        raise ValueError(
            f"invalid argument value for '{name}': expected two bounds "
            f"but received {value!r}"
        )

    interval = cast(Tuple[float, float], value)
    lower = _as_cutoff_bound(interval[0], f"{name}[0]")
    upper = _as_cutoff_bound(interval[1], f"{name}[1]")
    if lower >= upper:
        raise ValueError(
            f"invalid argument value for '{name}': expected lower bound "
            f"smaller than upper bound but received {value!r}"
        )

    return lower, upper


def _as_optional_cutoff_interval(
    value: Union[Literal["auto"], Tuple[float, float]],
    name: str,
) -> Optional[Tuple[float, float]]:

    if isinstance(value, str):
        if value == "auto":
            return None
        raise ValueError(
            f"invalid argument value for '{name}': expected 'auto' "
            f"but received {value!r}"
        )

    return _as_cutoff_interval(value, name)


def _as_cutoff_bound(
    value: Any,
    name: str,
) -> float:

    if not isinstance(value, numbers.Real) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {float} but received {type(value)}"
        )

    value = float(value)
    if np.isnan(value):
        raise ValueError(
            f"invalid argument value for '{name}': expected real value "
            f"but received {value!r}"
        )

    return value


def _select_hvg_features(
    *,
    ordered_indices: np.ndarray,
    scores: np.ndarray,
    means: np.ndarray,
    n_features: Optional[int],
    score: Optional[Tuple[float, float]],
    mean: Tuple[float, float],
    method: HVGMethod,
    n_vars: int,
) -> _HVGSelection:

    if n_features is not None:
        return _select_hvg_by_rank(
            ordered_indices,
            n_features=n_features,
            n_vars=n_vars,
        )

    return _select_hvg_by_cutoff(
        ordered_indices=ordered_indices,
        scores=scores,
        means=means,
        score=score,
        mean=mean,
        method=method,
        n_vars=n_vars,
    )


def _select_hvg_by_rank(
    ordered_indices: np.ndarray,
    *,
    n_features: int,
    n_vars: int,
) -> _HVGSelection:

    n_eligible = int(ordered_indices.size)
    if n_eligible < n_features:
        warnings.warn(
            f"Requested {n_features} highly variable features, but only "
            f"{n_eligible} features have finite variability scores.",
            RuntimeWarning,
            stacklevel=2,
        )

    selected_indices = ordered_indices[:n_features]

    return _hvg_selection_from_ordered_indices(
        selected_indices,
        n_vars=n_vars,
        score_cutoff=None,
        mean_cutoff=None,
    )


def _select_hvg_by_cutoff(
    *,
    ordered_indices: np.ndarray,
    scores: np.ndarray,
    means: np.ndarray,
    score: Optional[Tuple[float, float]],
    mean: Tuple[float, float],
    method: HVGMethod,
    n_vars: int,
) -> _HVGSelection:

    score_cutoff = score if score is not None else _DEFAULT_SCORE_CUTOFFS[method]
    mean_cutoff = mean
    selected = (
        np.isfinite(scores)
        & np.isfinite(means)
        & (scores > score_cutoff[0])
        & (scores < score_cutoff[1])
        & (means > mean_cutoff[0])
        & (means < mean_cutoff[1])
    )
    selected_indices = ordered_indices[selected[ordered_indices]]

    return _hvg_selection_from_ordered_indices(
        selected_indices,
        n_vars=n_vars,
        score_cutoff=score_cutoff,
        mean_cutoff=mean_cutoff,
    )


def _hvg_selection_from_ordered_indices(
    selected_indices: np.ndarray,
    *,
    n_vars: int,
    score_cutoff: Optional[Tuple[float, float]],
    mean_cutoff: Optional[Tuple[float, float]],
) -> _HVGSelection:

    selected = np.zeros(n_vars, dtype=bool)
    selected[selected_indices] = True
    ranks = np.full(n_vars, np.nan, dtype=float)
    ranks[selected_indices] = np.arange(1, selected_indices.size + 1, dtype=float)

    return _HVGSelection(
        selected=selected,
        ranks=ranks,
        score_cutoff=score_cutoff,
        mean_cutoff=mean_cutoff,
    )
