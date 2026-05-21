#!/usr/bin/env python

import numpy as np
import pandas as pd
import pytest

import bonesistools as bt
from bonesistools.sctools.tools import _utils


def test_choose_mtx_representation_selects_x_layer_and_raw(mini_adata):
    mini_adata.raw = mini_adata

    assert np.array_equal(
        bt.sct.tl.choose_mtx_representation(mini_adata, copy=False),
        mini_adata.X,
    )
    assert np.array_equal(
        bt.sct.tl.choose_mtx_representation(
            mini_adata,
            layer="counts",
            copy=False,
        ),
        mini_adata.layers["counts"],
    )
    assert np.array_equal(
        bt.sct.tl.choose_mtx_representation(
            mini_adata,
            use_raw=True,
            copy=False,
        ),
        mini_adata.raw.X,
    )


def test_choose_mtx_representation_rejects_raw_and_layer(mini_adata):
    with pytest.raises(ValueError, match="invalid argument combination"):
        bt.sct.tl.choose_mtx_representation(
            mini_adata,
            use_raw=True,
            layer="counts",
        )


def test_choose_representation_truncates_dimensions(mini_adata):
    rep = bt.sct.tl.choose_representation(
        mini_adata,
        use_rep="X_pca",
        n_components=2,
    )

    assert rep.shape == (4, 2)
    assert np.array_equal(rep, mini_adata.obsm["X_pca"][:, :2])


def test_choose_representation_reports_missing_key(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found in scdata.obsm"):
        bt.sct.tl.choose_representation(mini_adata, use_rep="missing")


def test_get_distances_and_connectivities_from_default_neighbors(mini_adata):
    assert _utils._get_distances(mini_adata) is mini_adata.obsp["distances"]
    assert _utils._get_connectivities(mini_adata) is mini_adata.obsp["connectivities"]


def test_anndata_to_dataframe_uses_layer_obs_and_log_base(mini_adata):
    mini_adata.uns["log1p"] = {"base": 2}
    mini_adata.layers["log2_counts"] = np.log1p(mini_adata.layers["counts"]) / np.log(2)

    df = bt.sct.tl.anndata_to_dataframe(
        mini_adata,
        obs=["cluster", "batch"],
        layer="log2_counts",
        is_log=True,
    )

    assert df.columns.tolist() == ["g1", "g2", "g3", "cluster", "batch"]
    assert np.allclose(df.loc[:, ["g1", "g2", "g3"]], mini_adata.layers["counts"])
    assert df["batch"].tolist() == ["b1", "b2", "b1", "b2"]


def test_pairwise_distances_returns_or_stores_matrix(mini_adata):
    distances = bt.sct.tl.pairwise_distances(
        mini_adata,
        use_rep="X_pca",
        n_components=2,
    )

    assert distances.shape == (4, 4)
    assert np.allclose(np.diag(distances), 0.0)

    result = bt.sct.tl.pairwise_distances(
        mini_adata,
        use_rep="X_pca",
        key_added="custom_distances",
    )
    assert result is None
    assert "custom_distances" in mini_adata.obsp
    assert mini_adata.uns["custom_distances"]["use_rep"] == "X_pca"


def test_barycenters_returns_cluster_means(mini_adata):
    barycenters = bt.sct.tl.barycenters(mini_adata, obs="cluster", use_rep="X_pca")

    assert sorted(barycenters) == ["A", "B"]
    assert np.allclose(barycenters["A"], np.array([0.1, 0.05, 1.05]))
    assert np.allclose(barycenters["B"], np.array([2.1, 2.05, 0.05]))


def test_barycenters_requires_categorical_obs(mini_adata):
    mini_adata.obs["plain"] = ["a", "a", "b", "b"]

    with pytest.raises(AttributeError, match="has no attribute 'cat'"):
        bt.sct.tl.barycenters(mini_adata, obs="plain")


def test_calculate_logfoldchanges_returns_non_empty_dataframe(mini_adata):
    result = bt.sct.tl.calculate_logfoldchanges(
        mini_adata,
        groupby="cluster",
    )

    assert set(result.columns) == {
        "group",
        "names",
        "logfoldchanges",
    }

    assert not result.empty
    assert set(result["group"]).issubset({"A", "B"})
    assert np.isfinite(result["logfoldchanges"]).all()


def test_calculate_logfoldchanges_filtering(mini_adata):
    filtered = bt.sct.tl.calculate_logfoldchanges(
        mini_adata,
        groupby="cluster",
        filter_logfoldchanges=lambda values: values > 0,
    )

    assert not filtered.empty
    assert (filtered["logfoldchanges"] > 0).all()


def test_calculate_logfoldchanges_rejects_invalid_filter(mini_adata):
    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'filter_logfoldchanges'",
    ):
        bt.sct.tl.calculate_logfoldchanges(
            mini_adata,
            groupby="cluster",
            filter_logfoldchanges="not callable",
        )


def test_hypergeometric_test_returns_probability(mini_adata):
    pvalue = bt.sct.tl.hypergeometric_test(
        mini_adata,
        signature=["g1", "g2"],
        markers=["g1", "g3"],
    )

    assert 0.0 <= pvalue <= 1.0


def test_smirnov_tests_stores_results(mini_adata):
    bt.sct.tl.smirnov_tests(
        mini_adata,
        groupby="cluster",
        groups=["A"],
        key_added="ks",
    )

    result = mini_adata.uns["ks"]["results"]
    assert result["group"].unique().tolist() == ["A"]
    assert {"statistics", "locations", "signs", "pvals", "pvals_adj"}.issubset(
        result.columns
    )


def test_smirnov_tests_copy_returns_dataframe(mini_adata):
    result = bt.sct.tl.smirnov_tests(
        mini_adata,
        groupby="cluster",
        groups=["A"],
        reference="B",
        corr_method="bonferroni",
        copy=True,
    )

    assert isinstance(result, pd.DataFrame)
    assert "pvals_adj" in result
    assert "smirnov_tests" not in mini_adata.uns
