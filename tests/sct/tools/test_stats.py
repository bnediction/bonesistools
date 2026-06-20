#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools._stats import _column_mean_variance
from bonesistools.sctools._typing import Matrix
from bonesistools.sctools.tools import _stats


def _import_scanpy():
    try:
        import scanpy as sc
    except Exception as error:
        pytest.skip(f"scanpy is not importable: {error}")
    return sc


def _toy_wilcoxon_adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 0.0, 2.0],
                [3.0, 0.0, 2.0],
                [2.0, 5.0, 2.0],
                [4.0, 10.0, 9.0],
                [5.0, 20.0, 8.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B", "C"])},
            index=["c1", "c2", "c3", "c4", "c5"],
        ),
        var=pd.DataFrame(index=["Gata1", "Klf1", "Tal1"]),
    )
    adata.layers["counts"] = csr_matrix(adata.X)
    return adata


def test_column_mean_variance_supports_sparse_arrays():
    matrix = sparse.csr_array(
        np.array(
            [
                [1.0, 0.0, 2.0],
                [3.0, 0.0, 4.0],
                [5.0, 6.0, 0.0],
            ]
        )
    )

    means, variances = _column_mean_variance(matrix)

    dense_matrix = matrix.toarray()
    np.testing.assert_allclose(means, dense_matrix.mean(axis=0))
    np.testing.assert_allclose(variances, dense_matrix.var(axis=0, ddof=1))


def _previous_fixed_background_loop(
    adata: ad.AnnData,
    groups,
    background: str,
) -> pd.DataFrame:
    labels = cast(pd.Series, adata.obs["cluster"])
    results = []
    for group in groups:
        group_mask = np.asarray(labels == group, dtype=bool)
        background_mask = np.asarray(labels == background, dtype=bool)
        combined_mask = group_mask | background_mask
        group_result = _stats._wilcoxon_test_matrix(
            cast(Any, adata.X)[combined_mask, :],
            group_mask[combined_mask],
            adata.var_names,
            tie_correction=True,
            max_memory=10_000_000_000,
        )
        group_result.insert(0, "group", group)
        group_result["pvals_adj"] = _stats._adjust_pvalues(
            group_result["pvals"].to_numpy(dtype=float),
            None,
        )
        results.append(group_result)

    df = pd.concat(results, axis=0)
    df.index.name = "names"
    return df[
        ["group", "statistics", "pvals", "pvals_adj", "u_statistics", "sum_ranks"]
    ]


def test_wilcoxon_tests_returns_expected_group_statistics():
    adata = _toy_wilcoxon_adata()

    result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background="rest",
        correction="bonferroni",
        tie_correction=False,
    )

    sigma = np.sqrt(2 * 3 * 6 / 12)
    expected_statistics = {
        "Gata1": (1.0 - 3.0) / sigma,
        "Klf1": (0.0 - 3.0) / sigma,
        "Tal1": (1.0 - 3.0) / sigma,
    }

    assert result.index.tolist() == ["Klf1", "Gata1", "Tal1"]
    assert result["group"].unique().tolist() == ["A"]
    assert result.loc["Gata1", "u_statistics"] == 1.0
    assert result.loc["Klf1", "u_statistics"] == 0.0
    assert result.loc["Tal1", "u_statistics"] == 1.0
    assert result.loc["Gata1", "sum_ranks"] == 4.0
    assert result.loc["Klf1", "sum_ranks"] == 3.0
    assert result.loc["Tal1", "sum_ranks"] == 4.0
    for gene, statistic in expected_statistics.items():
        assert result.loc[gene, "statistics"] == pytest.approx(statistic)
    assert (result["pvals_adj"] >= result["pvals"]).all()


def test_wilcoxon_tests_exposes_primary_statistics_column():
    result = bt.sct.tl.wilcoxon_tests(
        _toy_wilcoxon_adata(),
        groupby="cluster",
        groups=["A"],
    )

    assert "statistics" in result.columns


def test_wilcoxon_tests_matches_expected_tie_correction_values():
    group_a = np.array(
        [
            [5, 1, 0, 0, 1, 6, 0, 0, 0, 0],
            [4, 1, 0, 0, 2, 7, 1, 2, 0, 0],
            [4, 1, 1, 0, 2, 8, 2, 4, 1, 0],
            [3, 1, 1, 5, 3, 9, 3, 6, 1, 1],
            [3, 1, 2, 5, 3, 10, 4, 8, 1, 1],
        ],
        dtype=float,
    )
    group_b = np.array(
        [
            [0, 1, 5, 0, 1, 0, 6, 1, 0, 0],
            [0, 1, 4, 0, 1, 1, 7, 3, 0, 1],
            [1, 1, 4, 0, 2, 2, 8, 5, 0, 1],
            [1, 1, 3, 0, 2, 3, 9, 7, 0, 1],
            [2, 1, 3, 0, 3, 4, 10, 9, 1, 1],
        ],
        dtype=float,
    )
    adata = ad.AnnData(
        X=np.vstack([group_a, group_b]),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 5 + ["B"] * 5)},
            index=[f"c{i}" for i in range(10)],
        ),
        var=pd.DataFrame(index=[f"G{i}" for i in range(1, 11)]),
    )

    result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction=None,
    )

    expected = pd.DataFrame(
        {
            "u_statistics": [
                25.0,
                np.nan,
                0.0,
                17.5,
                16.0,
                25.0,
                0.0,
                10.0,
                17.5,
                7.5,
            ],
            "sum_ranks": [
                40.0,
                27.5,
                15.0,
                32.5,
                31.0,
                40.0,
                15.0,
                25.0,
                32.5,
                22.5,
            ],
            "statistics": [
                2.643,
                np.nan,
                -2.643,
                1.500,
                0.775,
                2.611,
                -2.611,
                -0.522,
                1.225,
                -1.225,
            ],
            "pvals": [
                0.0082,
                np.nan,
                0.0082,
                0.1336,
                0.4386,
                0.0090,
                0.0090,
                0.6015,
                0.2207,
                0.2207,
            ],
        },
        index=[f"G{i}" for i in range(1, 11)],
    )

    result = result.loc[expected.index]
    assert result["u_statistics"].to_numpy() == pytest.approx(
        expected["u_statistics"].to_numpy(),
        nan_ok=True,
    )
    assert result["sum_ranks"].to_numpy() == pytest.approx(
        expected["sum_ranks"].to_numpy(),
        nan_ok=True,
    )
    assert result["statistics"].to_numpy() == pytest.approx(
        expected["statistics"].to_numpy(),
        abs=5e-4,
        nan_ok=True,
    )
    assert result["pvals"].to_numpy() == pytest.approx(
        expected["pvals"].to_numpy(),
        abs=5e-5,
        nan_ok=True,
    )
    assert result["pvals_adj"].to_numpy() == pytest.approx(
        result["pvals"].to_numpy(),
        nan_ok=True,
    )


def test_wilcoxon_tests_warns_when_memory_is_smaller_than_one_rank_column():
    adata = _toy_wilcoxon_adata()

    with pytest.warns(RuntimeWarning):
        result = bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            groups=["A"],
            correction=None,
            max_memory=39,
        )

    expected = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        correction=None,
    )
    pd.testing.assert_frame_equal(result, expected)


def test_wilcoxon_tests_groups_all_background_and_var_subset():
    adata = _toy_wilcoxon_adata()

    result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups="all",
        background=["B"],
        var_subset=["Tal1", "Gata1"],
        correction=None,
    )

    assert result["group"].unique().tolist() == ["A", "C"]
    assert sorted(result.loc[result["group"] == "A"].index.tolist()) == [
        "Gata1",
        "Tal1",
    ]
    assert sorted(result.loc[result["group"] == "C"].index.tolist()) == [
        "Gata1",
        "Tal1",
    ]
    gata1_pvals_adj = cast(pd.Series, result.loc["Gata1", "pvals_adj"])
    gata1_pvals = cast(pd.Series, result.loc["Gata1", "pvals"])
    tal1_pvals_adj = cast(pd.Series, result.loc["Tal1", "pvals_adj"])
    tal1_pvals = cast(pd.Series, result.loc["Tal1", "pvals"])

    assert gata1_pvals_adj.tolist() == gata1_pvals.tolist()
    assert tal1_pvals_adj.tolist() == tal1_pvals.tolist()


def test_wilcoxon_tests_accepts_boolean_var_subset_column():
    adata = _toy_wilcoxon_adata()
    adata.var["selected"] = [True, False, True]

    from_column = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        var_subset="selected",
        correction=None,
    )
    from_names = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        var_subset=["Gata1", "Tal1"],
        correction=None,
    )

    pd.testing.assert_frame_equal(from_column, from_names)


def test_wilcoxon_tests_fixed_background_matches_previous_loop():
    adata = _toy_wilcoxon_adata()
    labels = cast(pd.Series, adata.obs["cluster"])

    old = _previous_fixed_background_loop(adata, groups=["A", "C"], background="B")
    new = pd.concat(
        _stats._wilcoxon_tests_fixed_background(
            cast(Matrix, adata.X),
            labels,
            ["A", "C"],
            ["B"],
            adata.var_names,
            tie_correction=True,
            max_memory=10_000_000_000,
            correction=None,
        ),
        axis=0,
    )
    new.index.name = "names"
    new = new[old.columns]

    pd.testing.assert_series_equal(new["group"], old["group"])
    np.testing.assert_allclose(
        new["statistics"].to_numpy(),
        old["statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        new["pvals"].to_numpy(),
        old["pvals"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        new["u_statistics"].to_numpy(),
        old["u_statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        new["sum_ranks"].to_numpy(),
        old["sum_ranks"].to_numpy(),
        equal_nan=True,
    )


def test_wilcoxon_tests_fixed_background_cannot_use_global_ranks():
    from scipy.stats import rankdata

    adata = ad.AnnData(
        X=np.array([[1.0], [10.0], [2.0], [3.0], [4.0], [5.0]]),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B", "C", "C"])},
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction=None,
    )

    labels = cast(pd.Series, adata.obs["cluster"])
    group_mask = np.asarray(labels == "A", dtype=bool)
    background_mask = np.asarray(labels == "B", dtype=bool)
    combined_mask = group_mask | background_mask
    expression_mtx = cast(np.ndarray, adata.X)
    subset_ranks = rankdata(
        expression_mtx[combined_mask, :],
        axis=0,
        method="average",
    )
    global_ranks = rankdata(expression_mtx, axis=0, method="average")

    expected_sum_ranks = subset_ranks[group_mask[combined_mask], :].sum(axis=0)[0]
    invalid_global_sum_ranks = global_ranks[group_mask, :].sum(axis=0)[0]

    assert result.loc["g1", "sum_ranks"] == expected_sum_ranks
    assert invalid_global_sum_ranks != expected_sum_ranks


def test_wilcoxon_tests_fixed_background_chunking_preserves_default_results():
    base = _toy_wilcoxon_adata()
    repeats = 50
    base_mtx = cast(np.ndarray, base.X)
    obs = cast(pd.DataFrame, base.obs).copy()
    adata = ad.AnnData(
        X=np.tile(base_mtx, (1, repeats)),
        obs=obs,
        var=pd.DataFrame(
            index=[
                f"{gene}_{repeat}"
                for repeat in range(repeats)
                for gene in base.var_names
            ]
        ),
    )

    unchunked = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A", "C"],
        background=["B"],
        correction=None,
    )
    small_memory = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A", "C"],
        background=["B"],
        correction=None,
        max_memory="1KB",
    )

    np.testing.assert_allclose(
        small_memory["statistics"].to_numpy(),
        unchunked["statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        small_memory["pvals"].to_numpy(),
        unchunked["pvals"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        small_memory["u_statistics"].to_numpy(),
        unchunked["u_statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        small_memory["sum_ranks"].to_numpy(),
        unchunked["sum_ranks"].to_numpy(),
        equal_nan=True,
    )


def test_wilcoxon_tests_supports_composite_background_groups():
    adata = _toy_wilcoxon_adata()

    composite = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B", "C"],
        correction=None,
    )
    rest = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background="rest",
        correction=None,
    )

    pd.testing.assert_frame_equal(composite, rest)


def test_wilcoxon_tests_distinguishes_rest_keyword_from_rest_group():
    adata = ad.AnnData(
        X=np.array([[10.0], [11.0], [0.0], [1.0], [-10.0], [-9.0]]),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "rest", "rest", "B", "B"])},
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    keyword_result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background="rest",
        correction=None,
    )
    group_result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["rest"],
        correction=None,
    )

    assert keyword_result.loc["g1", "sum_ranks"] != group_result.loc["g1", "sum_ranks"]


def test_wilcoxon_tests_supports_numeric_group_labels():
    adata = ad.AnnData(
        X=np.array([[5.0], [4.0], [0.0], [1.0]]),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical([0, 0, 1, 1])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=[0],
        background=[1],
        correction=None,
    )

    assert result.loc["g1", "group"] == 0
    assert result.loc["g1", "u_statistics"] == 4.0


def test_wilcoxon_tests_uses_sparse_layer():
    adata = _toy_wilcoxon_adata()

    dense_result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        expression="X",
    )
    sparse_result = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        expression="counts",
    )

    pd.testing.assert_frame_equal(sparse_result, dense_result)


def test_wilcoxon_tests_validates_arguments():
    adata = _toy_wilcoxon_adata()

    with pytest.raises(ValueError, match="invalid argument value for 'groups'"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", groups=["missing"])

    with pytest.raises(ValueError, match="invalid argument value for 'background'"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", background=["missing"])

    with pytest.raises(ValueError, match="must be disjoint"):
        bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            groups=["A"],
            background=["A"],
        )

    with pytest.raises(ValueError, match="must be disjoint"):
        bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            groups=["B", "A"],
            background=["B", "C"],
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'background'"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", background=cast(Any, "B"))

    with pytest.raises(ValueError, match="invalid argument value for 'background'"):
        bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            background=cast(Any, ["B", 1]),
        )

    with pytest.raises(ValueError, match="expected at least one background group"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", background=[])

    with pytest.raises(KeyError, match="variable"):
        bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            var_subset=["missing"],
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'groups'"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", groups=cast(Any, 1))

    with pytest.raises(TypeError, match="unsupported element type in 'var_subset'"):
        bt.sct.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            var_subset=cast(Any, ["Gata1", 1]),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'max_memory'"):
        bt.sct.tl.wilcoxon_tests(adata, groupby="cluster", max_memory="1XB")


def _toy_welch_test_adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 2.0, 0.0, 5.0],
                [2.0, 3.0, 0.0, 5.0],
                [3.0, 4.0, 1.0, 5.0],
                [4.0, 8.0, 2.0, 5.0],
                [5.0, 9.0, 2.0, 5.0],
                [6.0, 10.0, 3.0, 5.0],
            ],
            dtype=np.float64,
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "A", "B", "B", "B"])},
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4"]),
    )
    adata.layers["counts"] = csr_matrix(adata.X)
    return adata


def test_welch_tests_returns_expected_fixed_background_statistics():
    from scipy import stats

    adata = _toy_welch_test_adata()

    result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction=None,
    )

    group_mtx = np.asarray(adata.X[:3, :])
    background_mtx = np.asarray(adata.X[3:, :])
    mean_group = group_mtx.mean(axis=0, dtype=np.float64)
    mean_background = background_mtx.mean(axis=0, dtype=np.float64)
    variance_group = group_mtx.var(axis=0, ddof=1)
    variance_background = background_mtx.var(axis=0, ddof=1)
    statistics, pvals = stats.ttest_ind_from_stats(
        mean1=mean_group,
        std1=np.sqrt(variance_group),
        nobs1=group_mtx.shape[0],
        mean2=mean_background,
        std2=np.sqrt(variance_background),
        nobs2=background_mtx.shape[0],
        equal_var=False,
    )
    statistics[np.isnan(statistics)] = 0
    pvals[np.isnan(pvals)] = 1

    result = result.loc[adata.var_names]
    assert result["group"].unique().tolist() == ["A"]
    np.testing.assert_allclose(result["statistics"].to_numpy(), statistics)
    np.testing.assert_allclose(result["pvals"].to_numpy(), pvals)
    np.testing.assert_allclose(result["pvals_adj"].to_numpy(), pvals)
    np.testing.assert_allclose(result["mean_group"].to_numpy(), mean_group)
    np.testing.assert_allclose(
        result["mean_background"].to_numpy(),
        mean_background,
    )
    np.testing.assert_allclose(result["variance_group"].to_numpy(), variance_group)
    np.testing.assert_allclose(
        result["variance_background"].to_numpy(),
        variance_background,
    )


def test_welch_tests_exposes_primary_statistics_column():
    result = bt.sct.tl.welch_tests(
        _toy_welch_test_adata(),
        groupby="cluster",
        groups=["A"],
    )

    assert "statistics" in result.columns


def test_welch_tests_overestimates_background_variance_like_scanpy():
    from scipy import stats

    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 2.0],
                [2.0, 4.0],
                [5.0, 5.0],
                [6.0, 7.0],
                [7.0, 8.0],
                [8.0, 10.0],
            ],
            dtype=np.float64,
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B", "B", "B"])},
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction=None,
        overestimate_variance=True,
    ).loc[adata.var_names]

    group_mtx = np.asarray(adata.X[:2, :])
    background_mtx = np.asarray(adata.X[2:, :])
    statistics, pvals = stats.ttest_ind_from_stats(
        mean1=group_mtx.mean(axis=0, dtype=np.float64),
        std1=np.sqrt(group_mtx.var(axis=0, ddof=1)),
        nobs1=2,
        mean2=background_mtx.mean(axis=0, dtype=np.float64),
        std2=np.sqrt(background_mtx.var(axis=0, ddof=1)),
        nobs2=2,
        equal_var=False,
    )

    np.testing.assert_allclose(result["statistics"].to_numpy(), statistics)
    np.testing.assert_allclose(result["pvals"].to_numpy(), pvals)


def test_welch_tests_uses_sparse_layer():
    adata = _toy_welch_test_adata()

    dense_result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        expression="X",
    )
    sparse_result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        expression="counts",
    )

    pd.testing.assert_frame_equal(sparse_result, dense_result)


def test_welch_tests_accepts_boolean_var_subset_column():
    adata = _toy_welch_test_adata()
    adata.var["selected"] = [True, False, True, False]

    from_column = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        var_subset="selected",
        correction=None,
    )
    from_names = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        var_subset=["G1", "G3"],
        correction=None,
    )

    pd.testing.assert_frame_equal(from_column, from_names)


def test_welch_tests_chunking_preserves_default_results():
    base = _toy_welch_test_adata()
    repeats = 40
    base_mtx = cast(np.ndarray, base.X)
    adata = ad.AnnData(
        X=np.tile(base_mtx, (1, repeats)),
        obs=cast(pd.DataFrame, base.obs).copy(),
        var=pd.DataFrame(
            index=[
                f"{gene}_{repeat}"
                for repeat in range(repeats)
                for gene in base.var_names
            ]
        ),
    )

    unchunked = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction="benjamini-hochberg",
    )
    chunked = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction="benjamini-hochberg",
        max_memory=128,
    )

    pd.testing.assert_frame_equal(chunked, unchunked)


def test_welch_tests_matches_scanpy_rank_genes_groups_t_test():
    sc = _import_scanpy()
    adata = _toy_welch_test_adata()
    scanpy_adata = adata.copy()

    sc.tl.rank_genes_groups(
        scanpy_adata,
        "cluster",
        groups=["A"],
        reference="rest",
        method="t-test",
        corr_method="benjamini-hochberg",
        use_raw=False,
        n_genes=adata.n_vars,
    )
    result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background="rest",
        correction="benjamini-hochberg",
    )

    scanpy_result = pd.DataFrame(
        {
            "statistics": scanpy_adata.uns["rank_genes_groups"]["scores"]["A"],
            "pvals": scanpy_adata.uns["rank_genes_groups"]["pvals"]["A"],
            "pvals_adj": scanpy_adata.uns["rank_genes_groups"]["pvals_adj"]["A"],
        },
        index=scanpy_adata.uns["rank_genes_groups"]["names"]["A"],
    )
    result = result.loc[scanpy_result.index]

    np.testing.assert_array_equal(
        result["statistics"].to_numpy(dtype=np.float32),
        scanpy_result["statistics"].to_numpy(),
    )
    np.testing.assert_allclose(
        result["pvals"].to_numpy(),
        scanpy_result["pvals"].to_numpy(),
    )
    np.testing.assert_allclose(
        result["pvals_adj"].to_numpy(),
        scanpy_result["pvals_adj"].to_numpy(),
    )


def test_welch_tests_matches_scanpy_rank_genes_groups_overestimated_t_test():
    sc = _import_scanpy()
    adata = _toy_welch_test_adata()
    scanpy_adata = adata.copy()

    sc.tl.rank_genes_groups(
        scanpy_adata,
        "cluster",
        groups=["A"],
        reference="B",
        method="t-test_overestim_var",
        corr_method="bonferroni",
        use_raw=False,
        n_genes=adata.n_vars,
    )
    result = bt.sct.tl.welch_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        correction="bonferroni",
        overestimate_variance=True,
    )

    scanpy_result = pd.DataFrame(
        {
            "statistics": scanpy_adata.uns["rank_genes_groups"]["scores"]["A"],
            "pvals": scanpy_adata.uns["rank_genes_groups"]["pvals"]["A"],
            "pvals_adj": scanpy_adata.uns["rank_genes_groups"]["pvals_adj"]["A"],
        },
        index=scanpy_adata.uns["rank_genes_groups"]["names"]["A"],
    )
    result = result.loc[scanpy_result.index]

    np.testing.assert_array_equal(
        result["statistics"].to_numpy(dtype=np.float32),
        scanpy_result["statistics"].to_numpy(),
    )
    np.testing.assert_allclose(
        result["pvals"].to_numpy(),
        scanpy_result["pvals"].to_numpy(),
    )
    np.testing.assert_allclose(
        result["pvals_adj"].to_numpy(),
        scanpy_result["pvals_adj"].to_numpy(),
    )


def test_welch_tests_validates_arguments():
    adata = _toy_welch_test_adata()

    with pytest.raises(ValueError, match="invalid argument value for 'groups'"):
        bt.sct.tl.welch_tests(adata, groupby="cluster", groups=["missing"])

    with pytest.raises(ValueError, match="invalid argument value for 'background'"):
        bt.sct.tl.welch_tests(adata, groupby="cluster", background=["missing"])

    with pytest.raises(ValueError, match="must be disjoint"):
        bt.sct.tl.welch_tests(adata, groupby="cluster", groups=["A"], background=["A"])

    with pytest.raises(TypeError, match="unsupported argument type for 'background'"):
        bt.sct.tl.welch_tests(adata, groupby="cluster", background=cast(Any, "B"))

    with pytest.raises(ValueError, match="invalid argument value for 'max_memory'"):
        bt.sct.tl.welch_tests(adata, groupby="cluster", max_memory="1XB")

    with pytest.raises(
        ValueError,
        match="target group must contain at least two observations",
    ):
        small_group = adata.copy()
        small_group.obs["cluster"] = pd.Categorical(["A", "B", "B", "B", "B", "B"])
        bt.sct.tl.welch_tests(small_group, groupby="cluster", groups=["A"])
