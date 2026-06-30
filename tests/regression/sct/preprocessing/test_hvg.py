#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import bonesistools as bt
from bonesistools.sctools.preprocessing import _hvg


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


def test_hvg_dense_matrix_adds_expected_outputs():
    adata = _hvg_adata()

    result = bt.sct.pp.hvg(adata, n_features=5)
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
        },
    }


def test_hvg_sparse_matrix_matches_dense_scores_and_selection():
    dense = _hvg_adata(sparse_input=False)
    sparse_adata = _hvg_adata(sparse_input=True)

    bt.sct.pp.hvg(dense, n_features=5)
    bt.sct.pp.hvg(sparse_adata, n_features=5)

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


def test_hvg_ranking_follows_descending_positive_scores():
    adata = _hvg_adata()

    bt.sct.pp.hvg(adata, n_features=5)
    var = cast(pd.DataFrame, adata.var)
    selected = cast(pd.DataFrame, var.loc[var["highly_variable"]])
    ordered = selected.sort_values("highly_variable_rank")

    assert ordered.index[0] == var["highly_variable_score"].idxmax()
    assert ordered["highly_variable_rank"].tolist() == [1.0, 2.0, 3.0, 4.0, 5.0]
    assert ordered["highly_variable_score"].is_monotonic_decreasing


def test_hvg_non_selected_features_have_nan_rank_and_non_positive_scores_are_excluded():
    adata = _hvg_adata()

    bt.sct.pp.hvg(adata, n_features=5)
    var = cast(pd.DataFrame, adata.var)
    non_selected = cast(pd.DataFrame, var.loc[~var["highly_variable"]])

    assert np.isnan(non_selected["highly_variable_rank"].to_numpy()).all()
    assert (var.loc[var["highly_variable"], "highly_variable_score"] > 0).all()
    assert not (var["highly_variable"] & (var["highly_variable_score"] <= 0)).any()


def test_hvg_warns_when_fewer_features_have_finite_scores_than_requested():
    adata = ad.AnnData(
        X=np.ones((1, 5), dtype=float),
        var=pd.DataFrame(index=["g{}".format(index) for index in range(5)]),
    )

    with pytest.warns(RuntimeWarning):
        bt.sct.pp.hvg(adata, n_features=3)

    var = cast(pd.DataFrame, adata.var)
    assert not var["highly_variable"].any()
    assert np.isnan(var["highly_variable_rank"].to_numpy()).all()
    assert np.isnan(var["highly_variable_score"].to_numpy()).all()
    assert adata.uns["highly_variable"]["n_selected"] == 0


def test_hvg_copy_returns_modified_copy_and_preserves_original():
    adata = _hvg_adata()

    copied = bt.sct.pp.hvg(adata, n_features=4, copy=True)

    assert copied is not None
    assert "highly_variable" not in adata.var
    assert "highly_variable" in copied.var
    assert int(cast(pd.DataFrame, copied.var)["highly_variable"].sum()) == 4


def test_hvg_expression_layer_and_custom_key():
    adata = _hvg_adata()

    bt.sct.pp.hvg(
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
        bt.sct.pp.hvg(mini_adata, method=cast(Any, "seurat_v3"))

    with pytest.raises(ValueError):
        bt.sct.pp.hvg(mini_adata, n_features=0)

    with pytest.raises(ValueError):
        bt.sct.pp.hvg(mini_adata, span=0)

    with pytest.raises(ValueError):
        bt.sct.pp.hvg(mini_adata, span=1.5)

    with pytest.raises(TypeError):
        bt.sct.pp.hvg(mini_adata, span=cast(Any, "0.3"))

    with pytest.raises(TypeError):
        bt.sct.pp.hvg(mini_adata, key_added=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.hvg(mini_adata, copy=cast(Any, "yes"))


def test_hvg_loess_reports_missing_skmisc(mini_adata, monkeypatch):
    original_import_module = _hvg.importlib.import_module

    def fake_import_module(name):
        if name == "skmisc.loess":
            raise ImportError("missing skmisc")
        return original_import_module(name)

    monkeypatch.setattr(_hvg.importlib, "import_module", fake_import_module)

    with pytest.raises(ImportError, match="scikit-misc"):
        bt.sct.pp.hvg(mini_adata, method="loess")


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

    bt.sct.pp.hvg(adata, n_features=30, method="loess")
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
