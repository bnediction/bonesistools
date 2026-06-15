#!/usr/bin/env python

from __future__ import annotations

import warnings
from collections.abc import Sequence as SequenceInstance
from typing import Any, Optional, Sequence, Tuple, Union, cast

import numpy as np
import pandas as pd
from anndata import AnnData
from pandas import Index

from ..._compat import Literal
from ..._validation import (
    _as_boolean,
    _as_literal,
    _as_memory_size,
    _as_string,
)
from .._typing import Matrix, anndata_checker
from ._utils import _as_dense_matrix_chunk, _get_expression_with_gene_names

CorrectionMethod = Literal["benjamini-hochberg", "bonferroni"]

CORRECTION_METHODS: Tuple[CorrectionMethod, ...] = (
    "benjamini-hochberg",
    "bonferroni",
)


def _adjust_pvalues(
    pvals: np.ndarray,
    correction: Optional[CorrectionMethod],
) -> np.ndarray:

    adjusted_pvals = np.full_like(pvals, np.nan, dtype=float)
    valid = ~np.isnan(pvals)
    valid_pvals = pvals[valid]
    if valid_pvals.size == 0:
        return adjusted_pvals

    if correction is None:
        adjusted_pvals[valid] = valid_pvals
    elif correction == "bonferroni":
        adjusted_pvals[valid] = np.minimum(valid_pvals * valid_pvals.size, 1.0)
    elif correction == "benjamini-hochberg":
        order = np.argsort(valid_pvals)
        ranks = np.arange(1, valid_pvals.size + 1)
        adjusted = valid_pvals[order] * valid_pvals.size / ranks
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        valid_adjusted = np.empty_like(adjusted)
        valid_adjusted[order] = np.minimum(adjusted, 1.0)
        adjusted_pvals[valid] = valid_adjusted
    else:
        raise ValueError(
            f"invalid argument value for 'correction': "
            f"expected 'benjamini-hochberg', 'bonferroni' or None "
            f"but received {correction!r}"
        )

    return adjusted_pvals


@anndata_checker
def wilcoxon_tests(
    adata: AnnData,
    obs: str,
    groups: Union[Literal["all"], Sequence[Any]] = "all",
    background: Union[Literal["rest"], Sequence[Any]] = "rest",
    expression: Optional[str] = None,
    var_subset: Optional[Sequence[str]] = None,
    correction: Optional[CorrectionMethod] = "benjamini-hochberg",
    tie_correction: bool = True,
    max_memory: Union[int, str] = "10GB",
) -> pd.DataFrame:
    """
    Compute vectorized Wilcoxon rank-sum tests for groups of observations.

    For each target group, observations in `adata.obs[obs] == group` are
    compared with either the remaining observations or a named background group.
    Ranks are computed gene-wise in a single matrix operation. Ties are ranked
    with average ranks and tie correction is applied by default.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    obs: str
        Observation column in `adata.obs` defining groups.
    groups: 'all' or sequence (default: "all")
        Groups to test. If `"all"`, test all groups in `adata.obs[obs]`;
        when `background` is a sequence, background groups themselves are not
        tested.
    background: 'rest' or sequence (default: "rest")
        Background population. If `"rest"`, each target group is compared
        against all observations outside that group. If a sequence is provided,
        each target group is compared against the union of observations
        belonging to those background labels. To compare against an actual group
        named `"rest"`, pass `background=["rest"]`.
    expression: str, optional
        Expression source. If None or `"X"`, use `adata.X`. If `"raw.X"`, use
        `adata.raw.X`. Otherwise, interpret as a layer key in `adata.layers`.
    var_subset: sequence of str, optional
        Variable names to test. If None, test all variables in the selected
        expression matrix.
    correction: {'benjamini-hochberg', 'bonferroni'} or None
        Multiple-testing correction applied independently for each group. If
        None, p-values are not adjusted.
    tie_correction: bool (default: True)
        Whether to apply the rank-sum variance correction for tied values.
    max_memory: int or str (default: "10GB")
        Approximate maximum memory allocated to dense rank matrices. Integers
        are interpreted as bytes. Human-readable strings such as `"512MB"`,
        `"2GB"` or `"1GiB"` are accepted.

    Returns
    -------
    DataFrame
        Long-form table indexed by gene name. The table contains:

        - `group`: tested group;
        - `statistics`: normal-approximation z statistics;
        - `pvals`: two-sided p-values;
        - `pvals_adj`: adjusted p-values;
        - `u_statistics`: Mann-Whitney U statistics for the target group.
        - `sum_ranks`: sum of target-group ranks.
    """

    obs = _as_string(obs, "obs")
    correction = _as_literal(
        correction,
        choices=CORRECTION_METHODS,
        name="correction",
        allow_none=True,
    )
    tie_correction = _as_boolean(tie_correction, "tie_correction")
    max_memory = _as_memory_size(max_memory, "max_memory")

    expression_mtx, gene_names = _get_expression_with_gene_names(
        adata,
        expression,
        var_subset,
    )
    if len(expression_mtx.shape) != 2:
        raise ValueError("invalid expression matrix: expected a two-dimensional matrix")

    labels = cast(pd.Series, adata.obs[obs])
    background_groups = _resolve_background_groups(labels, background)
    target_groups = _resolve_groups(labels, groups, background_groups)
    if background_groups is not None:
        overlap = set(target_groups) & set(background_groups)
        if overlap:
            raise ValueError(
                "invalid argument combination: target groups and background groups "
                f"must be disjoint, but overlap on {sorted(overlap)!r}"
            )

    if background_groups is None:
        results = _wilcoxon_tests_rest(
            expression_mtx,
            labels,
            target_groups,
            gene_names,
            tie_correction=tie_correction,
            max_memory=max_memory,
            correction=correction,
        )
    else:
        results = _wilcoxon_tests_fixed_background(
            expression_mtx,
            labels,
            target_groups,
            background_groups,
            gene_names,
            tie_correction=tie_correction,
            max_memory=max_memory,
            correction=correction,
        )

    if not results:
        return pd.DataFrame(
            columns=cast(
                Any,
                [
                    "group",
                    "statistics",
                    "pvals",
                    "pvals_adj",
                    "u_statistics",
                    "sum_ranks",
                ],
            )
        )

    df = pd.concat(results, axis=0)
    df = df[["group", "statistics", "pvals", "pvals_adj", "u_statistics", "sum_ranks"]]
    df.index.name = "names"

    return cast(
        pd.DataFrame,
        cast(Any, df).sort_values(by=["group", "pvals"], kind="mergesort"),
    )


def _wilcoxon_tests_rest(
    expression_mtx: Matrix,
    labels: pd.Series,
    target_groups: Sequence[Any],
    gene_names: Index,
    tie_correction: bool,
    max_memory: int,
    correction: Optional[CorrectionMethod],
) -> Sequence[pd.DataFrame]:

    from scipy.stats import norm, rankdata

    n_genes = expression_mtx.shape[1]
    group_masks = [np.asarray(labels == group, dtype=bool) for group in target_groups]
    n_groups = np.asarray([mask.sum() for mask in group_masks], dtype=int)
    n_backgrounds = expression_mtx.shape[0] - n_groups
    if np.any(n_groups == 0):
        raise ValueError("invalid group selection: target group is empty")
    if np.any(n_backgrounds == 0):
        raise ValueError("invalid group selection: background group is empty")

    group_indicators = np.column_stack(group_masks).astype(np.uint8)
    sum_ranks = np.empty((len(target_groups), n_genes), dtype=float)
    sigma = np.empty((len(target_groups), n_genes), dtype=float)
    chunk_size = _rank_chunk_size(expression_mtx, max_memory)

    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        expression_chunk_mtx = _as_dense_matrix_chunk(expression_mtx, start, end)
        ranks = rankdata(expression_chunk_mtx, axis=0, method="average")
        sum_ranks[:, start:end] = group_indicators.T @ ranks
        tie_sums = _rank_sum_tie_sums(expression_chunk_mtx) if tie_correction else None
        for group_index, (n_group, n_background) in enumerate(
            zip(n_groups, n_backgrounds)
        ):
            sigma[group_index, start:end] = _rank_sum_sigma_from_tie_sums(
                expression_chunk_mtx.shape[1],
                int(n_group),
                int(n_background),
                tie_sums,
            )

    results = []
    for group_index, group in enumerate(target_groups):
        n_group = int(n_groups[group_index])
        n_background = int(n_backgrounds[group_index])
        u_statistics = sum_ranks[group_index, :] - n_group * (n_group + 1) / 2
        mu = n_group * n_background / 2
        group_sigma = sigma[group_index, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            statistics = (u_statistics - mu) / group_sigma
        u_statistics = np.where(group_sigma == 0, np.nan, u_statistics)
        statistics = np.where(group_sigma == 0, np.nan, statistics)
        pvals = 2 * norm.sf(np.abs(statistics))
        group_result = pd.DataFrame(
            {
                "group": group,
                "statistics": statistics,
                "pvals": pvals,
                "u_statistics": u_statistics,
                "sum_ranks": sum_ranks[group_index, :],
            },
            index=gene_names,
        )
        group_result["pvals_adj"] = _adjust_pvalues(
            group_result["pvals"].to_numpy(dtype=float),
            correction,
        )
        results.append(group_result)

    return results


def _wilcoxon_tests_fixed_background(
    expression_mtx: Matrix,
    labels: pd.Series,
    target_groups: Sequence[Any],
    background_groups: Sequence[Any],
    gene_names: Index,
    tie_correction: bool,
    max_memory: int,
    correction: Optional[CorrectionMethod],
) -> Sequence[pd.DataFrame]:

    background_mask = np.asarray(labels.isin(background_groups), dtype=bool)
    n_background = int(background_mask.sum())
    if n_background == 0:
        raise ValueError("invalid group selection: background group is empty")

    group_masks = {
        group: np.asarray(labels == group, dtype=bool) for group in target_groups
    }

    results = []
    for group, group_mask in group_masks.items():
        if int(group_mask.sum()) == 0:
            raise ValueError("invalid group selection: target group is empty")

        combined_mask = group_mask | background_mask
        group_result = _wilcoxon_test_matrix(
            cast(Any, expression_mtx)[combined_mask, :],
            group_mask[combined_mask],
            gene_names,
            tie_correction=tie_correction,
            max_memory=max_memory,
        )
        group_result.insert(0, "group", group)
        group_result["pvals_adj"] = _adjust_pvalues(
            group_result["pvals"].to_numpy(dtype=float),
            correction,
        )
        results.append(group_result)

    return results


def _wilcoxon_test_matrix(
    expression_mtx: Matrix,
    mask: np.ndarray,
    gene_names: Index,
    tie_correction: bool,
    max_memory: int,
) -> pd.DataFrame:

    from scipy.stats import norm, rankdata

    group_mask = np.asarray(mask, dtype=bool)
    n_group = int(group_mask.sum())
    n_background = int((~group_mask).sum())
    if n_group == 0:
        raise ValueError("invalid group selection: target group is empty")
    if n_background == 0:
        raise ValueError("invalid group selection: background group is empty")

    n_genes = expression_mtx.shape[1]
    chunk_size = _rank_chunk_size(expression_mtx, max_memory)
    sum_ranks = np.empty(n_genes, dtype=float)
    u_statistics = np.empty(n_genes, dtype=float)
    sigma = np.empty(n_genes, dtype=float)
    mu = n_group * n_background / 2
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        expression_chunk_mtx = _as_dense_matrix_chunk(expression_mtx, start, end)
        ranks = rankdata(expression_chunk_mtx, axis=0, method="average")
        sum_ranks[start:end] = ranks[group_mask, :].sum(axis=0)
        u_statistics[start:end] = sum_ranks[start:end] - n_group * (n_group + 1) / 2
        tie_sums = _rank_sum_tie_sums(expression_chunk_mtx) if tie_correction else None
        sigma[start:end] = _rank_sum_sigma_from_tie_sums(
            expression_chunk_mtx.shape[1],
            n_group=n_group,
            n_background=n_background,
            tie_sums=tie_sums,
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        statistics = (u_statistics - mu) / sigma
    u_statistics = np.where(sigma == 0, np.nan, u_statistics)
    statistics = np.where(sigma == 0, np.nan, statistics)
    pvals = 2 * norm.sf(np.abs(statistics))

    return pd.DataFrame(
        {
            "statistics": statistics,
            "pvals": pvals,
            "u_statistics": u_statistics,
            "sum_ranks": sum_ranks,
        },
        index=gene_names,
    )


def _resolve_groups(
    labels: pd.Series,
    groups: Union[Literal["all"], Sequence[Any]],
    background_groups: Optional[Sequence[Any]],
) -> Sequence[Any]:

    if isinstance(groups, str) and groups == "all":
        if isinstance(labels.dtype, pd.CategoricalDtype):
            resolved = list(labels.cat.categories)
        else:
            resolved = list(labels.dropna().unique())
        if background_groups is not None:
            background_group_set = set(background_groups)
            resolved = [
                group for group in resolved if group not in background_group_set
            ]
        return resolved

    if isinstance(groups, str) or not isinstance(groups, SequenceInstance):
        raise TypeError(
            f"unsupported argument type for 'groups': "
            f"expected 'all' or a sequence but received {type(groups)}"
        )

    resolved = list(groups)
    available_groups = set(labels.dropna().unique())
    missing = [group for group in resolved if group not in available_groups]
    if missing:
        formatted_missing = ", ".join(repr(group) for group in missing)
        raise ValueError(
            f"invalid argument value for 'groups': unknown group(s) {formatted_missing}"
        )

    return resolved


def _resolve_background_groups(
    labels: pd.Series,
    background: Union[Literal["rest"], Sequence[Any]],
) -> Optional[Sequence[Any]]:

    if isinstance(background, str) and background == "rest":
        return None

    if isinstance(background, str) or not isinstance(background, SequenceInstance):
        raise TypeError(
            f"unsupported argument type for 'background': "
            f"expected 'rest' or a sequence but received {type(background)}"
        )

    background_groups = list(background)
    if not background_groups:
        raise ValueError(
            "invalid argument value for 'background': "
            "expected at least one background group"
        )

    available_groups = set(labels.dropna().unique())
    missing = [group for group in background_groups if group not in available_groups]
    if missing:
        formatted_missing = ", ".join(repr(group) for group in missing)
        raise ValueError(
            f"invalid argument value for 'background': "
            f"unknown background group(s) {formatted_missing}"
        )

    return background_groups


def _rank_chunk_size(expression_mtx: Matrix, max_memory: int) -> int:

    bytes_per_rank = np.dtype(np.float64).itemsize
    bytes_per_rank_column = expression_mtx.shape[0] * bytes_per_rank
    if max_memory < bytes_per_rank_column:
        warnings.warn(
            "Requested memory budget is smaller than the memory required for a "
            "single gene. Computation will proceed with chunk_size=1.",
            RuntimeWarning,
            stacklevel=2,
        )
    return max(1, int(max_memory // bytes_per_rank_column))


def _rank_sum_tie_sums(expression_mtx: np.ndarray) -> np.ndarray:

    return np.asarray(
        [
            np.sum(counts**3 - counts)
            for counts in (
                np.unique(expression_mtx[:, gene_index], return_counts=True)[1]
                for gene_index in range(expression_mtx.shape[1])
            )
        ],
        dtype=float,
    )


def _rank_sum_sigma_from_tie_sums(
    n_genes: int,
    n_group: int,
    n_background: int,
    tie_sums: Optional[np.ndarray],
) -> np.ndarray:

    n_total = n_group + n_background
    if tie_sums is None:
        return np.full(
            n_genes,
            np.sqrt(n_group * n_background * (n_total + 1) / 12),
            dtype=float,
        )

    variance = (
        n_group
        * n_background
        / 12
        * ((n_total + 1) - tie_sums / (n_total * (n_total - 1)))
    )
    return np.sqrt(np.maximum(variance, 0))
