#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, Dict, Optional, Tuple, cast

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE, spectral_embedding

import bonesistools as bt


def _spectral_reference(
    graph,
    random_state,
    *,
    eigen_solver="arpack",
    eigen_tolerance=0.0,
):

    if eigen_solver == "lobpcg":
        graph = graph.astype(np.float64)
    return spectral_embedding(
        graph,
        n_components=2,
        random_state=random_state,
        eigen_solver=eigen_solver,
        eigen_tol=eigen_tolerance,
    )


def test_spectral_stores_graph_embedding_and_metadata(mini_adata):
    expected = _spectral_reference(
        mini_adata.obsp["connectivities"],
        random_state=0,
    )

    result = bt.omics.tl.spectral(
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
        "eigen_solver": "arpack",
        "eigen_tolerance": "auto",
        "n_jobs": 1,
    }


def test_spectral_copy_does_not_mutate_input(mini_adata):
    copied = bt.omics.tl.spectral(
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
    bt.omics.tl.spectral(
        mini_adata,
        n_components=2,
        key_added="X_spectral",
        seed=0,
    )

    assert mini_adata.obsm["X_spectral"].shape == (mini_adata.n_obs, 2)
    assert mini_adata.uns["X_spectral"]["method"] == "spectral"


def test_spectral_accepts_random_state_object(mini_adata):
    expected = _spectral_reference(
        mini_adata.obsp["connectivities"],
        random_state=np.random.RandomState(7),
    )

    bt.omics.tl.spectral(
        mini_adata,
        n_components=2,
        seed=np.random.RandomState(7),
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "RandomState"


def test_spectral_accepts_numpy_random_module(mini_adata):
    np.random.seed(2)
    expected = _spectral_reference(
        mini_adata.obsp["connectivities"],
        random_state=np.random.mtrand._rand,
    )

    np.random.seed(2)
    bt.omics.tl.spectral(
        mini_adata,
        n_components=2,
        seed=np.random,
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "np.random"


def test_spectral_validates_neighbors_and_seed(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found"):
        bt.omics.tl.spectral(
            mini_adata,
            neighbors_key="missing",
        )

    with pytest.raises(ValueError, match="invalid argument value for 'seed'"):
        bt.omics.tl.spectral(mini_adata, seed=cast(Any, "bad"))

    with pytest.raises(TypeError, match="unsupported argument type for 'eigen_solver'"):
        bt.omics.tl.spectral(mini_adata, eigen_solver=cast(Any, None))

    with pytest.raises(
        ValueError,
        match="invalid argument value for 'eigen_tolerance'",
    ):
        bt.omics.tl.spectral(mini_adata, eigen_tolerance=cast(Any, "invalid"))

    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'eigen_tolerance'",
    ):
        bt.omics.tl.spectral(mini_adata, eigen_tolerance=cast(Any, None))

    with pytest.raises(
        ValueError,
        match="invalid argument value for 'eigen_tolerance'",
    ):
        bt.omics.tl.spectral(mini_adata, eigen_tolerance=-1.0)

    with pytest.raises(ValueError, match="expected finite value"):
        bt.omics.tl.spectral(mini_adata, eigen_tolerance=np.nan)


def test_spectral_resolves_solver_tolerances(monkeypatch):
    observation_number = 20
    graph = csr_matrix(
        (observation_number, observation_number),
        dtype=np.float32,
    )
    adata = AnnData(np.zeros((observation_number, 1), dtype=np.float32))
    adata.obsp["connectivities"] = graph
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}

    parameters: Dict[str, Any] = {}

    def fake_spectral_embedding(
        matrix,
        *,
        n_components,
        eigen_solver,
        random_state,
        eigen_tol,
    ):
        parameters.update(
            {
                "n_components": n_components,
                "random_state": random_state,
                "eigen_solver": eigen_solver,
                "eigen_tol": eigen_tol,
                "dtype": matrix.dtype,
            }
        )
        return np.zeros((matrix.shape[0], n_components))

    monkeypatch.setattr(
        "sklearn.manifold.spectral_embedding",
        fake_spectral_embedding,
    )

    bt.omics.tl.spectral(adata, seed=0)

    assert parameters["eigen_solver"] == "arpack"
    assert parameters["eigen_tol"] == 0.0
    assert parameters["dtype"] == np.dtype(np.float32)
    assert adata.uns["X_se"]["eigen_solver"] == "arpack"
    assert adata.uns["X_se"]["eigen_tolerance"] == "auto"

    bt.omics.tl.spectral(
        adata,
        eigen_solver="lobpcg",
        eigen_tolerance="auto",
        seed=0,
    )

    assert parameters["eigen_solver"] == "lobpcg"
    assert parameters["eigen_tol"] is None
    assert parameters["dtype"] == np.dtype(np.float64)

    bt.omics.tl.spectral(
        adata,
        eigen_solver="lobpcg",
        eigen_tolerance=1e-10,
        seed=0,
    )

    assert parameters["eigen_tol"] == 1e-10
    assert adata.uns["X_se"]["eigen_tolerance"] == 1e-10


def test_tsne_embedding_stores_coordinates_and_metadata(mini_adata):
    result = bt.omics.tl.tsne(
        mini_adata,
        n_components=2,
        perplexity=1.0,
        n_iter=250,
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
        "n_iter": 250,
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

    bt.omics.tl.tsne(
        mini_adata,
        n_components=2,
        perplexity=1.0,
        n_iter=250,
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

    result = bt.omics.tl.umap(
        mini_adata,
        n_components=2,
        min_dist=0.2,
        spread=1.5,
        n_iter=17,
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
    assert FakeUMAP.parameters["init"] == "random"
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
        "n_iter": 17,
        "alpha": 0.75,
        "gamma": 1.25,
        "negative_sample_rate": 7,
        "init_pos": "random",
        "a": 1.7,
        "b": 0.9,
        "n_jobs": 2,
    }


def test_umap_requires_neighbors_graph_by_default(mini_adata):
    del mini_adata.uns["neighbors"]

    with pytest.raises(KeyError, match="bt.omics.tl.neighbors"):
        bt.omics.tl.umap(mini_adata)


def test_tsne_validates_perplexity(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'perplexity'"):
        bt.omics.tl.tsne(mini_adata, perplexity=mini_adata.n_obs)


def test_pca_stores_scores_loadings_and_metadata(mini_adata):
    expected = PCA(
        n_components=2,
        svd_solver="full",
        random_state=0,
    ).fit(mini_adata.X)

    result = bt.omics.tl.pca(
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
    loadings = cast(np.ndarray, mini_adata.varm["PCs"])
    assert np.allclose(
        np.abs(loadings),
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

    copied = bt.omics.tl.pca(
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
    loadings = cast(np.ndarray, copied.varm["PCs"])
    assert loadings.shape == (mini_adata.n_vars, 2)
    assert np.allclose(loadings[1, :], 0.0)
    assert copied.uns["pca"]["params"]["layer"] == "counts"
    assert copied.uns["pca"]["params"]["var_subset"] == "highly_variable"


def test_pca_accepts_variable_name_subset(mini_adata):
    bt.omics.tl.pca(
        mini_adata,
        n_components=2,
        var_subset=["g1", "g3"],
        svd_solver="full",
        seed=0,
    )

    assert mini_adata.obsm["X_pca"].shape == (mini_adata.n_obs, 2)
    loadings = cast(np.ndarray, mini_adata.varm["PCs"])
    assert np.allclose(loadings[1, :], 0.0)
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

    bt.omics.tl.pca(
        mini_adata,
        n_components=2,
        zero_center=False,
        seed=10,
    )

    assert np.allclose(mini_adata.obsm["X_pca"], expected_scores)
    loadings = cast(np.ndarray, mini_adata.varm["PCs"])
    assert np.allclose(loadings, expected.components_.T)
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

    bt.omics.tl.pca(
        mini_adata,
        n_components=2,
        zero_center=True,
        svd_solver="randomized",
        seed=10,
    )

    assert np.allclose(np.abs(mini_adata.obsm["X_pca"]), np.abs(expected_scores))
    loadings = cast(np.ndarray, mini_adata.varm["PCs"])
    assert np.allclose(
        np.abs(loadings),
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
    bt.omics.tl.pca(
        bonesis_adata,
        n_components=2,
        zero_center=zero_center,
        var_subset="highly_variable",
        seed=10,
        copy=False,
    )

    bonesis_scores = cast(np.ndarray, bonesis_adata.obsm["X_pca"])
    scanpy_scores = cast(np.ndarray, scanpy_adata.obsm["X_pca"])
    assert np.allclose(
        np.abs(bonesis_scores),
        np.abs(scanpy_scores),
    )
    bonesis_loadings = cast(np.ndarray, bonesis_adata.varm["PCs"])
    scanpy_loadings = cast(np.ndarray, scanpy_adata.varm["PCs"])
    assert np.allclose(
        np.abs(bonesis_loadings),
        np.abs(scanpy_loadings),
    )
    assert np.allclose(
        bonesis_adata.uns["pca"]["variance"],
        scanpy_adata.uns["pca"]["variance"],
    )
    assert np.allclose(
        bonesis_adata.uns["pca"]["variance_ratio"],
        scanpy_adata.uns["pca"]["variance_ratio"],
    )


def test_pca_preserves_dense_var_subset_layout():
    rng = np.random.RandomState(42)
    X = rng.normal(size=(200, 80)).astype(np.float32)
    obs = pd.DataFrame(index=["c" + str(i) for i in range(X.shape[0])])
    var = pd.DataFrame(index=["g" + str(i) for i in range(X.shape[1])])
    mask = np.arange(X.shape[1]) % 3 == 0
    var["highly_variable"] = mask

    bonesis_adata = AnnData(X.copy(), obs=obs.copy(), var=var.copy())
    expected = PCA(
        n_components=15,
        svd_solver="arpack",
        random_state=np.random.RandomState(10),
    )
    expected_scores = expected.fit_transform(X[:, mask])
    expected_loadings = np.zeros((X.shape[1], 15), dtype=expected.components_.dtype)
    expected_loadings[mask, :] = expected.components_.T

    bt.omics.tl.pca(
        bonesis_adata,
        n_components=15,
        zero_center=True,
        var_subset="highly_variable",
        seed=10,
        copy=False,
    )

    scores = cast(np.ndarray, bonesis_adata.obsm["X_pca"])
    assert np.array_equal(
        np.abs(scores),
        np.abs(expected_scores),
    )
    loadings = cast(np.ndarray, bonesis_adata.varm["PCs"])
    assert np.array_equal(
        np.abs(loadings),
        np.abs(expected_loadings),
    )
    assert np.array_equal(
        bonesis_adata.uns["pca"]["variance"],
        expected.explained_variance_,
    )
    assert np.array_equal(
        bonesis_adata.uns["pca"]["variance_ratio"],
        expected.explained_variance_ratio_,
    )


def test_pca_validates_arguments(mini_adata):
    with pytest.raises(KeyError, match="highly_variable"):
        bt.omics.tl.pca(mini_adata, n_components=2, var_subset="highly_variable")

    mini_adata.var["not_boolean"] = ["yes", "no", "yes"]
    with pytest.raises(TypeError):
        bt.omics.tl.pca(mini_adata, n_components=2, var_subset="not_boolean")

    with pytest.raises(KeyError):
        bt.omics.tl.pca(mini_adata, n_components=2, var_subset=["missing"])

    with pytest.raises(ValueError, match="invalid argument value for 'svd_solver'"):
        bt.omics.tl.pca(mini_adata, n_components=2, svd_solver=cast(Any, "bad"))

    with pytest.raises(ValueError, match="zero_center=False"):
        bt.omics.tl.pca(
            mini_adata, n_components=2, zero_center=False, svd_solver="full"
        )

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.omics.tl.pca(mini_adata, n_components=mini_adata.n_vars + 1)

    with pytest.raises(TypeError, match="unsupported argument type for 'key_added'"):
        bt.omics.tl.pca(mini_adata, n_components=2, key_added=cast(Any, 1))
