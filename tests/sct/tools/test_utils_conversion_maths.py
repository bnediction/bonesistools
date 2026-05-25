#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _utils


def _toy_logfoldchange_adata():
    return ad.AnnData(
        X=np.array(
            [
                [4.0, 2.0, 1.0],
                [4.0, 2.0, 1.0],
                [1.0, 4.0, 1.0],
                [1.0, 4.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B"])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


def _toy_smirnov_adata():
    return ad.AnnData(
        X=np.array(
            [
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B"])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["separated", "identical"]),
    )


def test_choose_mtx_representation_selects_x_layer_and_raw(mini_adata):
    mini_adata.raw = mini_adata

    x_copy = bt.sct.tl.choose_mtx_representation(mini_adata)

    assert np.array_equal(x_copy, mini_adata.X)
    assert x_copy is not mini_adata.X
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
    full_representation = bt.sct.tl.choose_representation(mini_adata, use_rep=None)
    rep = bt.sct.tl.choose_representation(
        mini_adata,
        use_rep="X_pca",
        n_components=2,
    )

    assert full_representation is mini_adata.obsm["X_pca"]
    assert rep.shape == (4, 2)
    assert np.array_equal(rep, mini_adata.obsm["X_pca"][:, :2])


def test_choose_representation_reports_missing_key(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found in scdata.obsm"):
        bt.sct.tl.choose_representation(mini_adata, use_rep="missing")


def test_get_distances_and_connectivities_from_default_neighbors(mini_adata):
    assert _utils._get_distances(mini_adata) is mini_adata.obsp["distances"]
    assert _utils._get_connectivities(mini_adata) is mini_adata.obsp["connectivities"]


def test_get_distances_and_connectivities_from_explicit_keys(mini_adata):
    mini_adata.obsp["custom_distances"] = mini_adata.obsp["distances"].copy()
    mini_adata.obsp["custom_connectivities"] = mini_adata.obsp["connectivities"].copy()
    mini_adata.uns["custom_neighbors"] = {
        "distances_key": "custom_distances",
        "connectivities_key": "custom_connectivities",
    }

    assert (
        _utils._get_distances(
            mini_adata,
            obsp="custom_distances",
        )
        is mini_adata.obsp["custom_distances"]
    )
    assert (
        _utils._get_connectivities(
            mini_adata,
            obsp="custom_connectivities",
        )
        is mini_adata.obsp["custom_connectivities"]
    )
    assert (
        _utils._get_distances(
            mini_adata,
            neighbors_key="custom_neighbors",
        )
        is mini_adata.obsp["custom_distances"]
    )
    assert (
        _utils._get_connectivities(
            mini_adata,
            neighbors_key="custom_neighbors",
        )
        is mini_adata.obsp["custom_connectivities"]
    )


@pytest.mark.parametrize(
    "getter,obsp_key,neighbors_key,error_message",
    [
        (
            _utils._get_distances,
            "distances",
            "neighbors",
            "'obsp' and 'neighbors_key' cannot be both specified",
        ),
        (
            _utils._get_connectivities,
            "connectivities",
            "neighbors",
            "'obsp' and 'neighbors_key' cannot be both specified",
        ),
    ],
)
def test_get_distances_and_connectivities_reject_ambiguous_keys(
    mini_adata,
    getter,
    obsp_key,
    neighbors_key,
    error_message,
):
    with pytest.raises(ValueError, match=error_message):
        getter(mini_adata, obsp=obsp_key, neighbors_key=neighbors_key)


@pytest.mark.parametrize(
    "getter,error_message",
    [
        (_utils._get_distances, "distances not found"),
        (_utils._get_connectivities, "connectivities not found"),
    ],
)
def test_get_distances_and_connectivities_report_missing_defaults(
    mini_adata,
    getter,
    error_message,
):
    mini_adata.uns.clear()

    with pytest.raises(KeyError, match=error_message):
        getter(mini_adata)


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


def test_anndata_to_dataframe_handles_sparse_x_sparse_layer_and_obs_string(mini_adata):
    sparse_adata = mini_adata.copy()
    sparse_adata.X = csr_matrix(sparse_adata.X)

    x_df = bt.sct.tl.anndata_to_dataframe(sparse_adata)

    assert np.allclose(x_df, mini_adata.X)

    sparse_adata.layers["sparse_log_counts"] = csr_matrix(
        np.log1p(mini_adata.layers["counts"])
    )
    layer_df = bt.sct.tl.anndata_to_dataframe(
        sparse_adata,
        obs="cluster",
        layer="sparse_log_counts",
        is_log=True,
    )

    assert layer_df.columns.tolist() == ["g1", "g2", "g3", "cluster"]
    assert np.allclose(
        layer_df.loc[:, ["g1", "g2", "g3"]],
        mini_adata.layers["counts"],
    )
    assert layer_df["cluster"].tolist() == ["A", "A", "B", "B"]


def test_pairwise_distances_returns_or_stores_matrix(
    mini_adata,
    expected_mini_pca2_distances,
):
    distances = bt.sct.tl.pairwise_distances(
        mini_adata,
        use_rep="X_pca",
        n_components=2,
    )

    assert distances.shape == (4, 4)
    assert np.allclose(distances, expected_mini_pca2_distances)

    result = bt.sct.tl.pairwise_distances(
        mini_adata,
        use_rep="X_pca",
        key_added="custom_distances",
    )
    assert result is None
    assert "custom_distances" in mini_adata.obsp
    assert mini_adata.uns["custom_distances"]["use_rep"] == "X_pca"


def test_barycenters_returns_cluster_means(
    mini_adata,
    expected_mini_cluster_barycenters,
):
    barycenters = bt.sct.tl.barycenters(mini_adata, obs="cluster", use_rep="X_pca")

    assert sorted(barycenters) == ["A", "B"]
    assert np.allclose(barycenters["A"], expected_mini_cluster_barycenters["A"])
    assert np.allclose(barycenters["B"], expected_mini_cluster_barycenters["B"])


def test_barycenters_requires_categorical_obs(mini_adata):
    mini_adata.obs["plain"] = ["a", "a", "b", "b"]

    with pytest.raises(AttributeError, match="has no attribute 'cat'"):
        bt.sct.tl.barycenters(mini_adata, obs="plain")


def test_calculate_logfoldchanges_recovers_expected_deseq2_style_contrasts():
    adata = _toy_logfoldchange_adata()

    result = bt.sct.tl.calculate_logfoldchanges(
        adata,
        groupby="cluster",
    )
    observed = result.set_index(["group", "names"])["logfoldchanges"].to_dict()

    # Simple two-group DESeq2-style contrast with size factors equal to 1:
    # log2(mean normalized counts in group / mean normalized counts in reference).
    assert observed == pytest.approx(
        {
            ("A", "g1"): 2.0,
            ("A", "g2"): -1.0,
            ("A", "g3"): 0.0,
            ("B", "g1"): -2.0,
            ("B", "g2"): 1.0,
            ("B", "g3"): 0.0,
        }
    )


def test_calculate_logfoldchanges_filtering_keeps_expected_positive_ratios():
    adata = _toy_logfoldchange_adata()

    filtered = bt.sct.tl.calculate_logfoldchanges(
        adata,
        groupby="cluster",
        filter_logfoldchanges=lambda values: values > 0,
    )
    observed = filtered.set_index(["group", "names"])["logfoldchanges"].to_dict()

    assert observed == pytest.approx(
        {
            ("A", "g1"): 2.0,
            ("B", "g2"): 1.0,
        }
    )


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


def test_hypergeometric_test_returns_expected_probability():
    adata = ad.AnnData(
        X=np.ones((1, 10)),
        var=pd.DataFrame(index=[f"g{i}" for i in range(1, 11)]),
    )

    pvalue = bt.sct.tl.hypergeometric_test(
        adata,
        signature=["g1", "g2", "g3", "g4"],
        markers=["g1", "g2", "g5"],
    )

    # P(X >= 2) for X ~ Hypergeometric(M=10, n=4, N=3):
    # (C(4, 2) C(6, 1) + C(4, 3) C(6, 0)) / C(10, 3) = 1/3.
    assert pvalue == pytest.approx(1 / 3)


def test_smirnov_tests_recovers_expected_ks_statistics():
    adata = _toy_smirnov_adata()

    bt.sct.tl.smirnov_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        reference="B",
        corr_method="bonferroni",
        key_added="ks",
    )

    result = adata.uns["ks"]["results"].set_index("names")
    assert result["group"].unique().tolist() == ["A"]
    assert result.loc["separated", "statistics"] == pytest.approx(1.0)
    assert result.loc["separated", "pvals"] == pytest.approx(1 / 3)
    assert result.loc["separated", "pvals_adj"] == pytest.approx(2 / 3)
    assert result.loc["identical", "statistics"] == pytest.approx(0.0)
    assert result.loc["identical", "pvals"] == pytest.approx(1.0)
    assert result.loc["identical", "pvals_adj"] == pytest.approx(1.0)


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
