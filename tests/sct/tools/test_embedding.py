#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, Dict, Optional, cast

import numpy as np
import pytest
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, SpectralEmbedding

import bonesistools as bt


def test_embedding_stores_spectral_embedding_and_metadata(mini_adata):
    expected = SpectralEmbedding(
        n_components=2,
        n_neighbors=3,
        random_state=0,
        n_jobs=1,
    ).fit_transform(mini_adata.obsm["X_pca"][:, :2])

    result = bt.sct.tl.embedding(
        mini_adata,
        method="spectral",
        use_rep="X_pca",
        n_rep_components=2,
        n_components=2,
        n_neighbors=3,
        seed=0,
        n_jobs=1,
    )

    assert result is None
    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"] == {
        "method": "spectral",
        "use_rep": "X_pca",
        "n_rep_components": 2,
        "n_components": 2,
        "n_neighbors": 3,
        "seed": 0,
        "eigen_solver": None,
        "n_jobs": 1,
    }


def test_embedding_copy_does_not_mutate_input(mini_adata):
    copied = bt.sct.tl.embedding(
        mini_adata,
        method="spectral",
        n_components=2,
        n_neighbors=3,
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
        n_neighbors=3,
        key_added="X_spectral",
        seed=0,
    )

    assert mini_adata.obsm["X_spectral"].shape == (mini_adata.n_obs, 2)
    assert mini_adata.uns["X_spectral"]["method"] == "spectral"


def test_embedding_accepts_random_state_object(mini_adata):
    expected = SpectralEmbedding(
        n_components=2,
        n_neighbors=3,
        random_state=np.random.RandomState(7),
        n_jobs=1,
    ).fit_transform(mini_adata.obsm["X_pca"])

    bt.sct.tl.embedding(
        mini_adata,
        method="spectral",
        n_components=2,
        n_neighbors=3,
        seed=np.random.RandomState(7),
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "RandomState"


def test_embedding_accepts_numpy_random_module(mini_adata):
    np.random.seed(2)
    expected = SpectralEmbedding(
        n_components=2,
        n_neighbors=3,
        random_state=np.random.mtrand._rand,
        n_jobs=1,
    ).fit_transform(mini_adata.obsm["X_pca"])

    np.random.seed(2)
    bt.sct.tl.embedding(
        mini_adata,
        method="spectral",
        n_components=2,
        n_neighbors=3,
        seed=np.random,
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_se"], expected)
    assert mini_adata.uns["X_se"]["seed"] == "np.random"


def test_embedding_validates_method_representation_and_neighbors(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'method'"):
        bt.sct.tl.embedding(mini_adata, method=cast(Any, "bad"))

    with pytest.raises(ValueError, match="invalid argument value for 'method'"):
        bt.sct.tl.embedding(mini_adata, method=cast(Any, "pca"))

    with pytest.raises(KeyError, match="key 'missing' not found"):
        bt.sct.tl.embedding(mini_adata, use_rep="missing", n_neighbors=3)

    with pytest.raises(ValueError, match="smaller than number of observations"):
        bt.sct.tl.embedding(mini_adata, n_neighbors=mini_adata.n_obs)

    with pytest.raises(ValueError, match="invalid argument value for 'seed'"):
        bt.sct.tl.embedding(mini_adata, n_neighbors=3, seed=cast(Any, "bad"))


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
        "use_rep": "X_pca",
        "n_rep_components": None,
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

    bt.sct.tl.embedding(
        mini_adata,
        method="tsne",
        n_components=2,
        perplexity=1.0,
        max_iter=250,
        seed=0,
        n_jobs=1,
    )

    assert np.allclose(mini_adata.obsm["X_tsne"], expected)


def test_umap_embedding_uses_umap_learn_and_stores_metadata(mini_adata, monkeypatch):
    class FakeUMAP:
        parameters: Dict[str, Any] = {}
        input_matrix: Optional[np.ndarray] = None

        def __init__(self, **kwargs):
            FakeUMAP.parameters = kwargs

        def fit_transform(self, X):
            FakeUMAP.input_matrix = X.copy()
            return np.full((X.shape[0], self.parameters["n_components"]), 2.0)

    module = ModuleType("umap")
    setattr(module, "UMAP", FakeUMAP)
    monkeypatch.setitem(sys.modules, "umap", module)

    result = bt.sct.tl.umap(
        mini_adata,
        n_rep_components=2,
        n_components=2,
        n_neighbors=3,
        min_dist=0.2,
        spread=1.5,
        max_iter=17,
        alpha=0.75,
        gamma=1.25,
        negative_sample_rate=7,
        metric="cosine",
        seed=0,
        n_jobs=2,
    )

    assert result is None
    input_matrix = FakeUMAP.input_matrix
    assert input_matrix is not None
    assert np.array_equal(input_matrix, mini_adata.obsm["X_pca"][:, :2])
    assert FakeUMAP.parameters["n_neighbors"] == 3
    assert FakeUMAP.parameters["metric"] == "cosine"
    assert FakeUMAP.parameters["min_dist"] == 0.2
    assert FakeUMAP.parameters["spread"] == 1.5
    assert FakeUMAP.parameters["n_epochs"] == 17
    assert FakeUMAP.parameters["learning_rate"] == 0.75
    assert FakeUMAP.parameters["repulsion_strength"] == 1.25
    assert FakeUMAP.parameters["negative_sample_rate"] == 7
    assert np.allclose(mini_adata.obsm["X_umap"], 2.0)
    assert mini_adata.uns["X_umap"] == {
        "method": "umap",
        "use_rep": "X_pca",
        "n_rep_components": 2,
        "n_components": 2,
        "seed": 0,
        "n_neighbors": 3,
        "metric": "cosine",
        "min_dist": 0.2,
        "spread": 1.5,
        "max_iter": 17,
        "alpha": 0.75,
        "gamma": 1.25,
        "negative_sample_rate": 7,
        "init_pos": "spectral",
        "a": None,
        "b": None,
        "n_jobs": 2,
    }


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
        "use_highly_variable": None,
        "layer": None,
        "use_raw": False,
        "svd_solver": "full",
        "key_added": "X_pca",
        "seed": 0,
    }


def test_pca_copy_layer_and_highly_variable_genes(mini_adata):
    mini_adata.var["highly_variable"] = [True, False, True]

    copied = bt.sct.tl.pca(
        mini_adata,
        n_components=2,
        layer="counts",
        use_highly_variable=True,
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
    assert copied.uns["pca"]["params"]["use_highly_variable"] is True


def test_pca_validates_arguments(mini_adata):
    with pytest.raises(KeyError, match="highly_variable"):
        bt.sct.tl.pca(mini_adata, n_components=2, use_highly_variable=True)

    with pytest.raises(ValueError, match="invalid argument value for 'svd_solver'"):
        bt.sct.tl.pca(mini_adata, n_components=2, svd_solver=cast(Any, "bad"))

    with pytest.raises(ValueError, match="zero_center=False"):
        bt.sct.tl.pca(mini_adata, n_components=2, zero_center=False, svd_solver="full")

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.tl.pca(mini_adata, n_components=mini_adata.n_vars + 1)

    with pytest.raises(TypeError, match="unsupported argument type for 'key_added'"):
        bt.sct.tl.pca(mini_adata, n_components=2, key_added=cast(Any, 1))
