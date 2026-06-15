#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, Dict, Optional, Tuple, cast

import numpy as np
import pytest
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE, SpectralEmbedding

import bonesistools as bt


def test_spectral_stores_graph_embedding_and_metadata(mini_adata):
    expected = SpectralEmbedding(
        n_components=2,
        affinity="precomputed",
        random_state=0,
        n_jobs=1,
    ).fit_transform(mini_adata.obsp["connectivities"])

    result = bt.sct.tl.spectral(
        mini_adata,
        n_components=2,
        seed=0,
        n_jobs=1,
    )

    assert result is None
    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"] == {
        "method": "spectral",
        "neighbors_key": "neighbors",
        "connectivities_key": "connectivities",
        "n_components": 2,
        "n_neighbors": 3,
        "metric": None,
        "seed": 0,
        "eigen_solver": None,
        "n_jobs": 1,
    }


def test_spectral_copy_does_not_mutate_input(mini_adata):
    copied = bt.sct.tl.spectral(
        mini_adata,
        n_components=2,
        key_added="X_custom",
        seed=0,
        copy=True,
    )

    assert "X_custom" not in mini_adata.obsm
    assert copied is not None
    assert copied.obsm["X_custom"].shape == (mini_adata.n_obs, 2)


def test_spectral_wrapper_stores_spectral_embedding(mini_adata):
    bt.sct.tl.spectral(
        mini_adata,
        n_components=2,
        key_added="X_spectral",
        seed=0,
    )

    assert mini_adata.obsm["X_spectral"].shape == (mini_adata.n_obs, 2)
    assert mini_adata.uns["X_spectral"]["method"] == "spectral"


def test_spectral_accepts_random_state_object(mini_adata):
    expected = SpectralEmbedding(
        n_components=2,
        affinity="precomputed",
        random_state=np.random.RandomState(7),
        n_jobs=1,
    ).fit_transform(mini_adata.obsp["connectivities"])

    bt.sct.tl.spectral(
        mini_adata,
        n_components=2,
        seed=np.random.RandomState(7),
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "RandomState"


def test_spectral_accepts_numpy_random_module(mini_adata):
    np.random.seed(2)
    expected = SpectralEmbedding(
        n_components=2,
        affinity="precomputed",
        random_state=np.random.mtrand._rand,
        n_jobs=1,
    ).fit_transform(mini_adata.obsp["connectivities"])

    np.random.seed(2)
    bt.sct.tl.spectral(
        mini_adata,
        n_components=2,
        seed=np.random,
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "np.random"


def test_spectral_validates_neighbors_and_seed(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found"):
        bt.sct.tl.spectral(
            mini_adata,
            neighbors_key="missing",
        )

    with pytest.raises(ValueError, match="invalid argument value for 'seed'"):
        bt.sct.tl.spectral(mini_adata, seed=cast(Any, "bad"))


def test_tsne_embedding_stores_coordinates_and_metadata(mini_adata):
    result = bt.sct.tl.tsne(
        mini_adata,
        n_components=2,
        perplexity=1.0,
        max_iter=250,
        seed=0,
        n_jobs=1,
    )

    assert result is None
    assert mini_adata.obsm["X_tsne"].shape == (mini_adata.n_obs, 2)
    assert mini_adata.uns["X_tsne"] == {
        "method": "tsne",
        "representation": "X_pca",
        "n_pcs": None,
        "n_components": 2,
        "seed": 0,
        "perplexity": 1.0,
        "metric": "euclidean",
        "early_exaggeration": 12.0,
        "learning_rate": 1000.0,
        "max_iter": 250,
        "n_jobs": 1,
    }


def test_tsne_embedding_matches_sklearn(mini_adata):
    expected = TSNE(
        n_components=2,
        perplexity=1.0,
        early_exaggeration=12.0,
        learning_rate=1000.0,
        max_iter=250,
        metric="euclidean",
        random_state=np.random.RandomState(0),
        n_jobs=1,
        method="barnes_hut",
    ).fit_transform(mini_adata.obsm["X_pca"])

    bt.sct.tl.tsne(
        mini_adata,
        n_components=2,
        perplexity=1.0,
        max_iter=250,
        seed=0,
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_tsne"], expected)


def test_umap_embedding_uses_neighbors_graph_and_stores_metadata(
    mini_adata,
    monkeypatch,
):
    class FakeUMAP:
        parameters: Dict[str, Any] = {}
        input_matrix: Optional[np.ndarray] = None
        graph_shape: Optional[Tuple[int, int]] = None

    def fake_find_ab_params(spread, min_dist):
        FakeUMAP.parameters["find_ab_params"] = (spread, min_dist)
        return 1.7, 0.9

    def fake_simplicial_set_embedding(**kwargs):
        FakeUMAP.parameters.update(kwargs)
        FakeUMAP.input_matrix = kwargs["data"].copy()
        FakeUMAP.graph_shape = kwargs["graph"].shape
        return np.full((kwargs["data"].shape[0], kwargs["n_components"]), 2.0), {}

    module = ModuleType("umap")
    setattr(module, "__path__", [])
    umap_module = ModuleType("umap.umap_")
    setattr(umap_module, "find_ab_params", fake_find_ab_params)
    setattr(umap_module, "simplicial_set_embedding", fake_simplicial_set_embedding)
    monkeypatch.setitem(sys.modules, "umap", module)
    monkeypatch.setitem(sys.modules, "umap.umap_", umap_module)

    result = bt.sct.tl.umap(
        mini_adata,
        n_components=2,
        min_dist=0.2,
        spread=1.5,
        max_iter=17,
        alpha=0.75,
        gamma=1.25,
        negative_sample_rate=7,
        seed=0,
        n_jobs=2,
    )

    assert result is None
    input_matrix = FakeUMAP.input_matrix
    assert input_matrix is not None
    assert np.array_equal(input_matrix, mini_adata.obsm["X_pca"][:, :2])
    assert FakeUMAP.graph_shape == (mini_adata.n_obs, mini_adata.n_obs)
    assert FakeUMAP.parameters["n_components"] == 2
    assert FakeUMAP.parameters["initial_alpha"] == 0.75
    assert FakeUMAP.parameters["a"] == 1.7
    assert FakeUMAP.parameters["b"] == 0.9
    assert FakeUMAP.parameters["gamma"] == 1.25
    assert FakeUMAP.parameters["negative_sample_rate"] == 7
    assert FakeUMAP.parameters["n_epochs"] == 17
    assert FakeUMAP.parameters["init"] == "spectral"
    assert FakeUMAP.parameters["metric"] == "euclidean"
    assert np.allclose(mini_adata.obsm["X_umap"], 2.0)
    assert mini_adata.uns["X_umap"] == {
        "method": "umap",
        "neighbors_key": "neighbors",
        "connectivities_key": "connectivities",
        "n_components": 2,
        "seed": 0,
        "n_neighbors": 3,
        "metric": "euclidean",
        "min_dist": 0.2,
        "spread": 1.5,
        "max_iter": 17,
        "alpha": 0.75,
        "gamma": 1.25,
        "negative_sample_rate": 7,
        "init_pos": "spectral",
        "a": 1.7,
        "b": 0.9,
        "n_jobs": 2,
    }


def test_umap_requires_neighbors_graph_by_default(mini_adata):
    del mini_adata.uns["neighbors"]

    with pytest.raises(KeyError, match="bt.sct.tl.neighbors"):
        bt.sct.tl.umap(mini_adata)


def test_tsne_validates_perplexity(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'perplexity'"):
        bt.sct.tl.tsne(mini_adata, perplexity=mini_adata.n_obs)


def test_pca_stores_scores_loadings_and_metadata(mini_adata):
    expected = PCA(
        n_components=2,
        svd_solver="full",
        random_state=0,
    ).fit(mini_adata.X)

    result = bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        svd_solver="full",
        seed=0,
    )

    assert result is None
    assert np.allclose(
        np.abs(mini_adata.obsm["X_pca"]),
        np.abs(expected.transform(mini_adata.X)),
    )
    assert np.allclose(
        np.abs(mini_adata.varm["PCs"]),
        np.abs(expected.components_.T),
    )
    assert np.allclose(mini_adata.uns["pca"]["variance"], expected.explained_variance_)
    assert np.allclose(
        mini_adata.uns["pca"]["variance_ratio"],
        expected.explained_variance_ratio_,
    )
    assert mini_adata.uns["pca"]["params"] == {
        "zero_center": True,
        "var_subset": None,
        "layer": None,
        "use_raw": False,
        "svd_solver": "full",
        "key_added": "X_pca",
        "seed": 0,
        "n_jobs": 1,
    }


def test_pca_copy_layer_and_highly_variable_genes(mini_adata):
    mini_adata.var["highly_variable"] = [True, False, True]

    copied = bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        layer="counts",
        var_subset="highly_variable",
        key_added="X_custom",
        svd_solver="full",
        seed=0,
        copy=True,
    )

    assert "X_custom" not in mini_adata.obsm
    assert copied is not None
    assert copied.obsm["X_custom"].shape == (mini_adata.n_obs, 2)
    assert copied.varm["PCs"].shape == (mini_adata.n_vars, 2)
    assert np.allclose(copied.varm["PCs"][1, :], 0.0)
    assert copied.uns["pca"]["params"]["layer"] == "counts"
    assert copied.uns["pca"]["params"]["var_subset"] == "highly_variable"


def test_pca_accepts_variable_name_subset(mini_adata):
    bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        var_subset=["g1", "g3"],
        svd_solver="full",
        seed=0,
    )

    assert mini_adata.obsm["X_pca"].shape == (mini_adata.n_obs, 2)
    assert np.allclose(mini_adata.varm["PCs"][1, :], 0.0)
    assert mini_adata.uns["pca"]["params"]["var_subset"] == ("g1", "g3")


def test_pca_without_zero_center_preserves_sparse_input(mini_adata):
    sparse_expression_mtx = csr_matrix(mini_adata.X)
    mini_adata.X = sparse_expression_mtx
    expected = TruncatedSVD(
        n_components=2,
        algorithm="randomized",
        random_state=np.random.RandomState(10),
    )
    expected_scores = expected.fit_transform(sparse_expression_mtx)

    bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        zero_center=False,
        seed=10,
    )

    assert np.allclose(mini_adata.obsm["X_pca"], expected_scores)
    assert np.allclose(mini_adata.varm["PCs"], expected.components_.T)
    assert np.allclose(mini_adata.uns["pca"]["variance"], expected.explained_variance_)
    assert np.allclose(
        mini_adata.uns["pca"]["variance_ratio"],
        expected.explained_variance_ratio_,
    )
    assert mini_adata.uns["pca"]["params"]["svd_solver"] == "randomized"


def test_pca_with_sparse_input_densifies_for_dense_solver(mini_adata):
    sparse_expression_mtx = csr_matrix(mini_adata.X.astype(np.float32))
    mini_adata.X = sparse_expression_mtx
    expected = PCA(
        n_components=2,
        svd_solver="randomized",
        random_state=np.random.RandomState(10),
    )
    expected_scores = expected.fit_transform(sparse_expression_mtx.toarray())

    bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        zero_center=True,
        svd_solver="randomized",
        seed=10,
    )

    assert np.allclose(np.abs(mini_adata.obsm["X_pca"]), np.abs(expected_scores))
    assert np.allclose(
        np.abs(mini_adata.varm["PCs"]),
        np.abs(expected.components_.T),
    )
    assert np.allclose(mini_adata.uns["pca"]["variance"], expected.explained_variance_)
    assert np.allclose(
        mini_adata.uns["pca"]["variance_ratio"],
        expected.explained_variance_ratio_,
    )
    assert mini_adata.uns["pca"]["params"]["svd_solver"] == "randomized"


@pytest.mark.parametrize("zero_center", [True, False])
def test_pca_matches_scanpy_with_sparse_input_and_seed(mini_adata, zero_center):
    sc = pytest.importorskip("scanpy")

    sparse_expression_mtx = csr_matrix(mini_adata.X.astype(np.float32))
    mini_adata.X = sparse_expression_mtx
    mini_adata.var["highly_variable"] = [True, True, True]
    scanpy_adata = mini_adata.copy()
    bonesis_adata = mini_adata.copy()

    sc.tl.pca(
        scanpy_adata,
        n_comps=2,
        zero_center=zero_center,
        mask_var="highly_variable",
        random_state=10,
        copy=False,
    )
    bt.sct.tl.pca(
        bonesis_adata,
        n_components=2,
        zero_center=zero_center,
        var_subset="highly_variable",
        seed=10,
        copy=False,
    )

    assert np.allclose(
        np.abs(bonesis_adata.obsm["X_pca"]),
        np.abs(scanpy_adata.obsm["X_pca"]),
    )
    assert np.allclose(
        np.abs(bonesis_adata.varm["PCs"]),
        np.abs(scanpy_adata.varm["PCs"]),
    )
    assert np.allclose(
        bonesis_adata.uns["pca"]["variance"],
        scanpy_adata.uns["pca"]["variance"],
    )
    assert np.allclose(
        bonesis_adata.uns["pca"]["variance_ratio"],
        scanpy_adata.uns["pca"]["variance_ratio"],
    )


def test_pca_validates_arguments(mini_adata):
    with pytest.raises(KeyError, match="highly_variable"):
        bt.sct.tl.pca(mini_adata, n_components=2, var_subset="highly_variable")

    mini_adata.var["not_boolean"] = ["yes", "no", "yes"]
    with pytest.raises(TypeError):
        bt.sct.tl.pca(mini_adata, n_components=2, var_subset="not_boolean")

    with pytest.raises(KeyError):
        bt.sct.tl.pca(mini_adata, n_components=2, var_subset=["missing"])

    with pytest.raises(ValueError, match="invalid argument value for 'svd_solver'"):
        bt.sct.tl.pca(mini_adata, n_components=2, svd_solver=cast(Any, "bad"))

    with pytest.raises(ValueError, match="zero_center=False"):
        bt.sct.tl.pca(mini_adata, n_components=2, zero_center=False, svd_solver="full")

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.tl.pca(mini_adata, n_components=mini_adata.n_vars + 1)

    with pytest.raises(TypeError, match="unsupported argument type for 'key_added'"):
        bt.sct.tl.pca(mini_adata, n_components=2, key_added=cast(Any, 1))
