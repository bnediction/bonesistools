#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import bonesistools as bt
from bonesistools.omics._typing import Matrix
from bonesistools.omics.preprocessing import _hvg


def _hvg_adata(sparse_input=False):

    rng = np.random.default_rng(0)
    n_obs = 80
    n_vars = 120
    means = np.linspace(0.2, 5.0, n_vars)
    matrix = rng.poisson(means, size=(n_obs, n_vars)).astype(float)
    matrix[:, 12] += rng.choice([0, 20], size=n_obs, p=[0.85, 0.15])
    matrix[:, 37] += rng.choice([0, 12], size=n_obs, p=[0.75, 0.25])
    matrix[:, 75] += rng.choice([0, 18], size=n_obs, p=[0.90, 0.10])
    matrix[:, 104] += rng.choice([0, 10], size=n_obs, p=[0.80, 0.20])

    X = sparse.csr_matrix(matrix) if sparse_input else matrix
    adata = ad.AnnData(
        X=X,
        var=pd.DataFrame(index=["g{}".format(index) for index in range(n_vars)]),
    )
    adata.layers["counts"] = (
        sparse.csr_matrix(matrix) if sparse_input else matrix.copy()
    )
    return adata


def _weighted_nanmean(values, weights):

    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    valid = ~np.isnan(values)
    weighted_values = np.where(valid, values, 0.0) * weights.reshape(-1, 1)
    denominators = np.sum(valid * weights.reshape(-1, 1), axis=0)

    return np.divide(
        weighted_values.sum(axis=0),
        denominators,
        out=np.full(values.shape[1], np.nan, dtype=np.float64),
        where=denominators != 0,
    )


def _top_score_indices(scores, n_features):

    ordered = np.flatnonzero(~np.isnan(scores))

    return ordered[np.argsort(-scores[ordered], kind="mergesort")][:n_features]


def test_hvg_dense_matrix_adds_expected_outputs():
    adata = _hvg_adata()

    result = bt.omics.pp.hvg(adata, n_features=5)
    var = cast(pd.DataFrame, adata.var)

    assert result is None
    assert "highly_variable" in var
    assert "highly_variable_rank" in var
    assert "highly_variable_score" in var
    assert var["highly_variable"].dtype == bool
    assert int(var["highly_variable"].sum()) == 5
    assert adata.uns["highly_variable"] == {
        "method": "loess",
        "n_features": 5,
        "n_selected": 5,
        "expression": None,
        "params": {
            "fit": "loess",
            "span": 0.3,
            "score": "normalized_variance_score",
            "selection": "top_n",
            "score_cutoff": None,
            "mean_cutoff": None,
            "batch_key": None,
            "batch_weighting": "equal",
            "batch_selection": "consensus",
        },
    }


def test_hvg_sparse_matrix_matches_dense_scores_and_selection():
    dense = _hvg_adata(sparse_input=False)
    sparse_adata = _hvg_adata(sparse_input=True)

    bt.omics.pp.hvg(dense, n_features=5)
    bt.omics.pp.hvg(sparse_adata, n_features=5)

    dense_var = cast(pd.DataFrame, dense.var)
    sparse_var = cast(pd.DataFrame, sparse_adata.var)
    assert np.array_equal(
        dense_var["highly_variable"].to_numpy(),
        sparse_var["highly_variable"].to_numpy(),
    )
    np.testing.assert_allclose(
        dense_var["highly_variable_rank"].to_numpy(),
        sparse_var["highly_variable_rank"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        dense_var["highly_variable_score"].to_numpy(),
        sparse_var["highly_variable_score"].to_numpy(),
    )


def test_hvg_binning_sparse_matrix_matches_dense_scores_and_selection():
    dense = _hvg_adata(sparse_input=False)
    sparse_adata = _hvg_adata(sparse_input=True)

    bt.omics.pp.hvg(dense, method="binning", n_features=5)
    bt.omics.pp.hvg(sparse_adata, method="binning", n_features=5)

    dense_var = cast(pd.DataFrame, dense.var)
    sparse_var = cast(pd.DataFrame, sparse_adata.var)
    assert np.array_equal(
        dense_var["highly_variable"].to_numpy(),
        sparse_var["highly_variable"].to_numpy(),
    )
    np.testing.assert_allclose(
        dense_var["highly_variable_rank"].to_numpy(),
        sparse_var["highly_variable_rank"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        dense_var["highly_variable_score"].to_numpy(),
        sparse_var["highly_variable_score"].to_numpy(),
        equal_nan=True,
    )


def test_hvg_batch_key_uses_equal_batch_weighting_by_default():
    adata = _hvg_adata()
    matrix = np.asarray(adata.X)
    batch = np.array(["small"] * 20 + ["large"] * 60)
    adata.obs["batch"] = batch

    bt.omics.pp.hvg(
        adata,
        method="binning",
        n_features=5,
        n_bins=10,
        batch_key="batch",
    )

    small_scores = _hvg._score_hvg_features(
        matrix[batch == "small", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    large_scores = _hvg._score_hvg_features(
        matrix[batch == "large", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    expected_scores = _weighted_nanmean(
        np.vstack([small_scores, large_scores]),
        np.ones(2),
    )
    nbatches = np.zeros(matrix.shape[1], dtype=int)
    nbatches[_top_score_indices(small_scores, 5)] += 1
    nbatches[_top_score_indices(large_scores, 5)] += 1
    eligible = np.flatnonzero(np.isfinite(expected_scores))
    expected_order = eligible[
        np.lexsort((-expected_scores[eligible], -nbatches[eligible]))
    ]
    expected_selected = np.zeros(matrix.shape[1], dtype=bool)
    expected_selected[expected_order[:5]] = True
    var = cast(pd.DataFrame, adata.var)

    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        expected_scores,
        equal_nan=True,
    )
    assert np.array_equal(var["highly_variable"].to_numpy(), expected_selected)
    np.testing.assert_allclose(
        var["highly_variable_rank"].to_numpy()[expected_order[:5]],
        np.arange(1, 6, dtype=float),
    )
    assert adata.uns["highly_variable"]["params"]["batch_key"] == "batch"
    assert adata.uns["highly_variable"]["params"]["batch_weighting"] == "equal"
    assert adata.uns["highly_variable"]["params"]["batch_selection"] == "consensus"


def test_hvg_batch_key_can_weight_batches_by_cell_count():
    adata = _hvg_adata()
    matrix = np.asarray(adata.X)
    batch = np.array(["small"] * 20 + ["large"] * 60)
    adata.obs["batch"] = batch

    bt.omics.pp.hvg(
        adata,
        method="binning",
        n_features=5,
        n_bins=10,
        batch_key="batch",
        batch_weighting="cell_count",
    )

    small_scores = _hvg._score_hvg_features(
        matrix[batch == "small", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    large_scores = _hvg._score_hvg_features(
        matrix[batch == "large", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    batch_scores = np.vstack([small_scores, large_scores])
    expected_scores = _weighted_nanmean(batch_scores, np.array([20, 60]))
    equal_scores = _weighted_nanmean(batch_scores, np.ones(2))
    nbatches = np.zeros(matrix.shape[1], dtype=int)
    nbatches[_top_score_indices(small_scores, 5)] += 1
    nbatches[_top_score_indices(large_scores, 5)] += 1
    eligible = np.flatnonzero(np.isfinite(expected_scores))
    expected_order = eligible[
        np.lexsort((-expected_scores[eligible], -nbatches[eligible]))
    ]
    expected_selected = np.zeros(matrix.shape[1], dtype=bool)
    expected_selected[expected_order[:5]] = True
    var = cast(pd.DataFrame, adata.var)
    finite = np.isfinite(expected_scores) & np.isfinite(equal_scores)

    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        expected_scores,
        equal_nan=True,
    )
    assert np.array_equal(var["highly_variable"].to_numpy(), expected_selected)
    assert np.max(np.abs(expected_scores[finite] - equal_scores[finite])) > 1e-12
    assert adata.uns["highly_variable"]["params"]["batch_key"] == "batch"
    assert adata.uns["highly_variable"]["params"]["batch_weighting"] == "cell_count"
    assert adata.uns["highly_variable"]["params"]["batch_selection"] == "consensus"


def test_hvg_batch_selection_rank_prioritizes_summarized_within_batch_rank():
    adata = _hvg_adata()
    matrix = np.asarray(adata.X)
    batch = np.array(["small"] * 20 + ["large"] * 60)
    adata.obs["batch"] = batch

    bt.omics.pp.hvg(
        adata,
        method="binning",
        n_features=5,
        n_bins=10,
        batch_key="batch",
        batch_selection="rank",
    )

    small_scores = _hvg._score_hvg_features(
        matrix[batch == "small", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    large_scores = _hvg._score_hvg_features(
        matrix[batch == "large", :],
        method="binning",
        span=0.3,
        n_bins=10,
    ).scores
    expected_scores = _weighted_nanmean(
        np.vstack([small_scores, large_scores]),
        np.ones(2),
    )
    ranks = np.full((2, matrix.shape[1]), np.nan, dtype=float)
    nbatches = np.zeros(matrix.shape[1], dtype=int)
    for batch_index, scores in enumerate([small_scores, large_scores]):
        selected_indices = _top_score_indices(scores, 5)
        nbatches[selected_indices] += 1
        ranks[batch_index, selected_indices] = np.arange(
            1,
            selected_indices.size + 1,
            dtype=float,
        )
    rank_summary = np.ma.median(
        np.ma.masked_invalid(ranks),
        axis=0,
    ).filled(np.nan)
    rank_sort_values = np.where(np.isnan(rank_summary), np.inf, rank_summary)
    eligible = np.flatnonzero(np.isfinite(expected_scores))
    expected_order = eligible[
        np.lexsort(
            (
                -expected_scores[eligible],
                -nbatches[eligible],
                rank_sort_values[eligible],
            )
        )
    ]
    consensus_order = eligible[
        np.lexsort(
            (
                -expected_scores[eligible],
                -nbatches[eligible],
            )
        )
    ]
    expected_selected = np.zeros(matrix.shape[1], dtype=bool)
    expected_selected[expected_order[:5]] = True
    consensus_selected = np.zeros(matrix.shape[1], dtype=bool)
    consensus_selected[consensus_order[:5]] = True
    var = cast(pd.DataFrame, adata.var)

    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        expected_scores,
    )
    assert np.array_equal(var["highly_variable"].to_numpy(), expected_selected)
    assert not np.array_equal(expected_selected, consensus_selected)
    np.testing.assert_allclose(
        var["highly_variable_rank"].to_numpy()[expected_order[:5]],
        np.arange(1, 6, dtype=float),
    )
    assert adata.uns["highly_variable"]["params"]["batch_selection"] == "rank"


def test_hvg_loess_batch_key_prioritizes_batches_then_median_rank():
    pytest.importorskip("skmisc")

    adata = _hvg_adata()
    matrix = np.asarray(adata.X)
    batch = np.array(["a"] * 40 + ["b"] * 40)
    adata.obs["batch"] = batch

    bt.omics.pp.hvg(
        adata,
        method="loess",
        n_features=8,
        batch_key="batch",
    )

    a_scores = _hvg._score_hvg_features(
        matrix[batch == "a", :],
        method="loess",
        span=0.3,
        n_bins=20,
    ).scores
    b_scores = _hvg._score_hvg_features(
        matrix[batch == "b", :],
        method="loess",
        span=0.3,
        n_bins=20,
    ).scores
    expected_scores = _weighted_nanmean(
        np.vstack([a_scores, b_scores]),
        np.ones(2),
    )
    ranks = np.full((2, matrix.shape[1]), np.nan, dtype=float)
    nbatches = np.zeros(matrix.shape[1], dtype=int)
    for batch_index, scores in enumerate([a_scores, b_scores]):
        selected_indices = _top_score_indices(scores, 8)
        nbatches[selected_indices] += 1
        ranks[batch_index, selected_indices] = np.arange(
            1,
            selected_indices.size + 1,
            dtype=float,
        )
    rank_summary = np.ma.median(
        np.ma.masked_invalid(ranks),
        axis=0,
    ).filled(np.nan)
    rank_sort_values = np.where(np.isnan(rank_summary), np.inf, rank_summary)
    eligible = np.flatnonzero(np.isfinite(expected_scores))
    expected_order = eligible[
        np.lexsort(
            (
                -expected_scores[eligible],
                rank_sort_values[eligible],
                -nbatches[eligible],
            )
        )
    ]
    expected_selected = np.zeros(matrix.shape[1], dtype=bool)
    expected_selected[expected_order[:8]] = True
    var = cast(pd.DataFrame, adata.var)

    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        expected_scores,
    )
    assert np.array_equal(var["highly_variable"].to_numpy(), expected_selected)
    np.testing.assert_allclose(
        var["highly_variable_rank"].to_numpy()[expected_order[:8]],
        np.arange(1, 9, dtype=float),
    )


def test_hvg_batch_rank_summary_supports_cell_count_weighting():
    ranks = np.array(
        [
            [1.0, 4.0, np.nan],
            [10.0, 2.0, 3.0],
            [20.0, np.nan, 1.0],
        ]
    )
    weights = np.array([1.0, 3.0, 6.0])

    summary = _hvg._batch_rank_summary(
        ranks,
        weights,
        batch_weighting="cell_count",
    )

    np.testing.assert_allclose(summary, np.array([20.0, 2.0, 1.0]))


def test_hvg_weighted_nanmedian_preserves_left_cutoff_and_missing_values():
    values = np.array(
        [
            [3.0, np.nan, 5.0, 0.0],
            [1.0, np.nan, 5.0, -np.inf],
            [2.0, np.nan, 7.0, np.inf],
            [np.nan, np.nan, 9.0, np.nan],
        ]
    )
    weights = np.array([2.0, 1.0, 3.0, 4.0])

    summary = _hvg._weighted_nanmedian(values, weights)

    assert np.array_equal(
        summary,
        np.array([2.0, np.nan, 7.0, 0.0]),
        equal_nan=True,
    )


def test_hvg_binning_scores_match_mean_binned_dispersion_formula():
    matrix = np.array(
        [
            [0, 1, 1, 2, 3, 3, 4, 5],
            [0, 1, 2, 2, 3, 5, 4, 8],
            [0, 1, 1, 4, 3, 3, 8, 5],
            [0, 1, 3, 2, 6, 3, 4, 5],
            [0, 2, 1, 2, 3, 6, 4, 9],
            [0, 1, 1, 5, 3, 3, 9, 5],
        ],
        dtype=float,
    )
    adata = ad.AnnData(
        X=matrix,
        var=pd.DataFrame(index=["g{}".format(index) for index in range(8)]),
    )

    bt.omics.pp.hvg(
        adata,
        method="binning",
        n_features=3,
        n_bins=2,
    )
    var = cast(pd.DataFrame, adata.var)

    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        np.array(
            [
                np.nan,
                -2.108749563831,
                0.131796847739,
                1.217182652653,
                -0.131796847739,
                -0.674489750196,
                3.599762976185,
                0.0,
            ]
        ),
        equal_nan=True,
    )
    assert var["highly_variable"].to_numpy().tolist() == [
        False,
        False,
        True,
        True,
        False,
        False,
        True,
        False,
    ]
    np.testing.assert_allclose(
        var["highly_variable_rank"].to_numpy(),
        np.array([np.nan, np.nan, 3.0, 2.0, np.nan, np.nan, 1.0, np.nan]),
        equal_nan=True,
    )
    assert adata.uns["highly_variable"] == {
        "method": "binning",
        "n_features": 3,
        "n_selected": 3,
        "expression": None,
        "params": {
            "fit": "mean_binned_dispersion",
            "n_bins": 2,
            "score": "normalized_dispersion_score",
            "selection": "top_n",
            "score_cutoff": None,
            "mean_cutoff": None,
            "batch_key": None,
            "batch_weighting": "equal",
            "batch_selection": "consensus",
        },
    }


def test_hvg_ranking_follows_descending_positive_scores():
    adata = _hvg_adata()

    bt.omics.pp.hvg(adata, n_features=5)
    var = cast(pd.DataFrame, adata.var)
    selected = cast(pd.DataFrame, var.loc[var["highly_variable"]])
    ordered = selected.sort_values("highly_variable_rank")

    assert ordered.index[0] == var["highly_variable_score"].idxmax()
    assert ordered["highly_variable_rank"].tolist() == [1.0, 2.0, 3.0, 4.0, 5.0]
    assert ordered["highly_variable_score"].is_monotonic_decreasing


def test_hvg_non_selected_features_have_nan_rank_and_non_positive_scores_are_excluded():
    adata = _hvg_adata()

    bt.omics.pp.hvg(adata, n_features=5)
    var = cast(pd.DataFrame, adata.var)
    non_selected = cast(pd.DataFrame, var.loc[~var["highly_variable"]])

    assert np.isnan(non_selected["highly_variable_rank"].to_numpy()).all()
    assert (var.loc[var["highly_variable"], "highly_variable_score"] > 0).all()
    assert not (var["highly_variable"] & (var["highly_variable_score"] <= 0)).any()


def test_hvg_cutoffs_select_features_by_score_and_mean():
    adata = _hvg_adata()

    bt.omics.pp.hvg(
        adata,
        n_features=None,
        score=(2.0, np.inf),
        mean=(0.0, 3.5),
    )
    var = cast(pd.DataFrame, adata.var)
    scores = var["highly_variable_score"].to_numpy()
    means = _hvg._compute_log_normalized_mean(cast(Matrix, adata.X))
    expected = np.isfinite(scores) & (scores > 2.0) & (means > 0.0) & (means < 3.5)

    assert np.array_equal(var["highly_variable"].to_numpy(), expected)
    selected = cast(pd.DataFrame, var.loc[var["highly_variable"]])
    ordered = selected.sort_values("highly_variable_rank")
    assert ordered["highly_variable_score"].is_monotonic_decreasing
    assert adata.uns["highly_variable"]["n_features"] is None
    assert adata.uns["highly_variable"]["params"]["selection"] == "cutoff"
    assert adata.uns["highly_variable"]["params"]["score_cutoff"] == (2.0, np.inf)
    assert adata.uns["highly_variable"]["params"]["mean_cutoff"] == (0.0, 3.5)


def test_hvg_uses_method_default_cutoffs_when_n_features_is_none():
    adata = _hvg_adata()

    bt.omics.pp.hvg(
        adata,
        n_features=None,
        score="auto",
    )
    var = cast(pd.DataFrame, adata.var)
    scores = var["highly_variable_score"].to_numpy()
    means = _hvg._compute_log_normalized_mean(cast(Matrix, adata.X))
    expected = np.isfinite(scores) & (scores > 2.0) & (means > 0.0125) & (means < 3.0)

    assert np.array_equal(var["highly_variable"].to_numpy(), expected)
    assert adata.uns["highly_variable"]["params"]["score_cutoff"] == (2.0, np.inf)
    assert adata.uns["highly_variable"]["params"]["mean_cutoff"] == (0.0125, 3.0)


def test_hvg_loess_mean_cutoff_uses_log_normalized_mean():
    adata = _hvg_adata()
    bt.omics.pp.normalize(adata, expression="counts", key_added="normalized")
    bt.omics.pp.log1p(adata, expression="normalized", key_added="log1p")

    bt.omics.pp.hvg(
        adata,
        expression="counts",
        n_features=None,
        score=(-np.inf, np.inf),
        mean=(0.5, 3.5),
    )
    var = cast(pd.DataFrame, adata.var)
    scores = var["highly_variable_score"].to_numpy()
    log_means = np.asarray(adata.layers["log1p"]).mean(axis=0)
    expected = np.isfinite(scores) & (log_means > 0.5)
    expected &= log_means < 3.5

    assert np.array_equal(var["highly_variable"].to_numpy(), expected)


def test_hvg_n_features_takes_priority_over_cutoffs():
    adata = _hvg_adata()

    with pytest.warns(UserWarning, match="cutoffs are ignored"):
        bt.omics.pp.hvg(
            adata,
            n_features=3,
            score=(1e9, np.inf),
            mean=(1e9, np.inf),
        )
    var = cast(pd.DataFrame, adata.var)

    assert int(var["highly_variable"].sum()) == 3
    assert adata.uns["highly_variable"]["params"]["selection"] == "top_n"
    assert adata.uns["highly_variable"]["params"]["score_cutoff"] is None
    assert adata.uns["highly_variable"]["params"]["mean_cutoff"] is None


def test_hvg_binning_cutoffs_select_features_by_score_and_mean():
    adata = _hvg_adata()

    bt.omics.pp.hvg(
        adata,
        method="binning",
        n_features=None,
        score=(1.0, np.inf),
        mean=(0.0, 4.0),
    )
    var = cast(pd.DataFrame, adata.var)
    scores = var["highly_variable_score"].to_numpy()
    means = np.asarray(adata.X).mean(axis=0)
    expected = np.isfinite(scores) & (scores > 1.0) & (means > 0.0) & (means < 4.0)

    assert np.array_equal(var["highly_variable"].to_numpy(), expected)
    assert adata.uns["highly_variable"]["params"]["selection"] == "cutoff"
    assert adata.uns["highly_variable"]["params"]["score_cutoff"] == (1.0, np.inf)
    assert adata.uns["highly_variable"]["params"]["mean_cutoff"] == (0.0, 4.0)


def test_hvg_binning_uses_method_default_cutoffs_when_n_features_is_none():
    adata = _hvg_adata()

    bt.omics.pp.hvg(adata, method="binning", n_features=None)
    var = cast(pd.DataFrame, adata.var)
    scores = var["highly_variable_score"].to_numpy()
    means = np.asarray(adata.X).mean(axis=0)
    expected = np.isfinite(scores) & (scores > 0.5) & (means > 0.0125) & (means < 3.0)

    assert np.array_equal(var["highly_variable"].to_numpy(), expected)
    assert adata.uns["highly_variable"]["params"]["score_cutoff"] == (0.5, np.inf)
    assert adata.uns["highly_variable"]["params"]["mean_cutoff"] == (0.0125, 3.0)


def test_hvg_binning_ignores_unexpressed_features_when_binning():
    rng = np.random.default_rng(0)
    matrix = rng.poisson(
        np.linspace(0.2, 5.0, 120),
        size=(80, 120),
    ).astype(float)
    with_zeros = np.column_stack(
        [
            np.zeros((matrix.shape[0], 3), dtype=float),
            matrix,
        ]
    )
    expressed = ad.AnnData(
        X=matrix,
        var=pd.DataFrame(index=["g{}".format(index) for index in range(120)]),
    )
    extended = ad.AnnData(
        X=with_zeros,
        var=pd.DataFrame(index=["z0", "z1", "z2"] + list(expressed.var_names)),
    )

    bt.omics.pp.hvg(expressed, method="binning", n_features=30)
    bt.omics.pp.hvg(extended, method="binning", n_features=30)
    expressed_var = cast(pd.DataFrame, expressed.var)
    extended_var = cast(pd.DataFrame, extended.var)

    assert not extended_var.iloc[:3]["highly_variable"].any()
    assert np.isnan(extended_var.iloc[:3]["highly_variable_score"].to_numpy()).all()
    assert np.array_equal(
        expressed_var["highly_variable"].to_numpy(),
        extended_var.iloc[3:]["highly_variable"].to_numpy(),
    )
    np.testing.assert_allclose(
        expressed_var["highly_variable_score"].to_numpy(),
        extended_var.iloc[3:]["highly_variable_score"].to_numpy(),
    )


def test_hvg_warns_when_fewer_features_have_finite_scores_than_requested():
    adata = ad.AnnData(
        X=np.ones((1, 5), dtype=float),
        var=pd.DataFrame(index=["g{}".format(index) for index in range(5)]),
    )

    with pytest.warns(RuntimeWarning):
        bt.omics.pp.hvg(adata, n_features=3)

    var = cast(pd.DataFrame, adata.var)
    assert not var["highly_variable"].any()
    assert np.isnan(var["highly_variable_rank"].to_numpy()).all()
    assert np.isnan(var["highly_variable_score"].to_numpy()).all()
    assert adata.uns["highly_variable"]["n_selected"] == 0


def test_hvg_inplace_false_returns_results_without_modifying_adata():
    adata = _hvg_adata()
    expected = adata.copy()

    result = bt.omics.pp.hvg(adata, n_features=4, inplace=False)
    bt.omics.pp.hvg(expected, n_features=4)
    expected_var = cast(pd.DataFrame, expected.var)

    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["selected", "rank", "score"]
    assert result.index.equals(adata.var_names)
    assert result["selected"].dtype == bool
    assert int(result["selected"].sum()) == 4
    np.testing.assert_array_equal(
        result["selected"].to_numpy(),
        expected_var["highly_variable"].to_numpy(),
    )
    np.testing.assert_allclose(
        result["rank"].to_numpy(),
        expected_var["highly_variable_rank"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result["score"].to_numpy(),
        expected_var["highly_variable_score"].to_numpy(),
    )
    assert "highly_variable" not in adata.var
    assert "highly_variable" not in adata.uns


def test_hvg_expression_layer_and_custom_key():
    adata = _hvg_adata()

    bt.omics.pp.hvg(
        adata,
        expression="counts",
        n_features=4,
        span=0.5,
        key_added="counts_hvg",
    )
    var = cast(pd.DataFrame, adata.var)

    assert "counts_hvg" in var
    assert "counts_hvg_rank" in var
    assert "counts_hvg_score" in var
    assert int(var["counts_hvg"].sum()) == 4
    assert adata.uns["counts_hvg"]["expression"] == "counts"
    assert adata.uns["counts_hvg"]["params"]["span"] == 0.5


def test_hvg_validates_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, method=cast(Any, "seurat_v3"))

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, n_features=0)

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, span=0)

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, span=1.5)

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, n_bins=0)

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, batch_key=cast(Any, 1))

    with pytest.raises(KeyError):
        bt.omics.pp.hvg(mini_adata, batch_key="unknown")

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, batch_weighting=cast(Any, "weighted"))

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, batch_selection=cast(Any, "score"))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, span=cast(Any, "0.3"))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, n_bins=cast(Any, "20"))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, score=cast(Any, 1.0))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, score=cast(Any, None))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, mean=cast(Any, "default"))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, mean=cast(Any, "auto"))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, mean=cast(Any, (0.0, "max")))

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, score=(1.0, 1.0))

    with pytest.raises(ValueError):
        bt.omics.pp.hvg(mini_adata, mean=(np.nan, 1.0))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, key_added=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.omics.pp.hvg(mini_adata, inplace=cast(Any, "yes"))


def test_hvg_loess_reports_missing_skmisc(mini_adata, monkeypatch):
    original_import_module = _hvg.importlib.import_module

    def fake_import_module(name):
        if name == "skmisc.loess":
            raise ImportError("missing skmisc")
        return original_import_module(name)

    monkeypatch.setattr(_hvg.importlib, "import_module", fake_import_module)

    with pytest.raises(ImportError, match="scikit-misc"):
        bt.omics.pp.hvg(mini_adata, method="loess")


def test_hvg_loess_matches_scanpy_seurat_v3():
    sc = pytest.importorskip("scanpy")
    pytest.importorskip("skmisc")

    rng = np.random.default_rng(0)
    means = np.linspace(0.1, 5.0, 120)
    theta = 8.0
    probability = theta / (theta + means.reshape(1, -1))
    matrix = rng.negative_binomial(
        theta,
        probability,
        size=(80, 120),
    ).astype(np.float32)
    adata = ad.AnnData(
        X=matrix,
        var=pd.DataFrame(index=["g{}".format(index) for index in range(120)]),
    )
    scanpy_adata = adata.copy()

    bt.omics.pp.hvg(adata, n_features=30, method="loess")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(
            scanpy_adata,
            n_top_genes=30,
            flavor="seurat_v3",
            inplace=True,
        )

    var = cast(pd.DataFrame, adata.var)
    scanpy_var = cast(pd.DataFrame, scanpy_adata.var)
    assert np.array_equal(
        var["highly_variable"].to_numpy(),
        scanpy_var["highly_variable"].to_numpy(),
    )
    np.testing.assert_allclose(
        var["highly_variable_score"].to_numpy(),
        scanpy_var["variances_norm"].to_numpy(),
    )
    np.testing.assert_allclose(
        var["highly_variable_rank"].to_numpy() - 1,
        scanpy_var["highly_variable_rank"].to_numpy(),
        equal_nan=True,
    )
