#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans

import bonesistools as bt


def _set_two_component_graph(adata, key):
    adata.obsp[key] = csr_matrix(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )


def test_kmeans_matches_sklearn_and_stores_metadata(mini_adata):
    result = bt.omics.tl.kmeans(
        mini_adata,
        n_clusters=2,
        representation="X_pca",
        n_pcs=2,
        n_init=3,
        algorithm="lloyd",
        key_added="kmeans",
        seed=10,
    )

    expected = cast(
        np.ndarray,
        cast(Any, KMeans)(
            n_clusters=2,
            init="k-means++",
            n_init=3,
            max_iter=300,
            tol=1e-4,
            random_state=np.random.RandomState(10),
            algorithm="lloyd",
            copy_x=True,
        ).fit_predict(mini_adata.obsm["X_pca"][:, :2]),
    )

    assert result is None
    clusters = mini_adata.obs["kmeans"]
    assert isinstance(clusters.dtype, pd.CategoricalDtype)
    assert clusters.astype(str).tolist() == [str(cluster) for cluster in expected]
    assert mini_adata.uns["kmeans"] == {
        "params": {
            "method": "kmeans",
            "n_clusters": 2,
            "representation": "X_pca",
            "n_pcs": 2,
            "n_init": 3,
            "seed": 10,
            "algorithm": "lloyd",
        }
    }


def test_kmeans_copy_and_custom_key(mini_adata):
    copied = bt.omics.tl.kmeans(
        mini_adata,
        n_clusters=2,
        representation="X_pca",
        n_pcs=1,
        n_init=2,
        key_added="clusters",
        seed=42,
        copy=True,
        init="random",
        max_iter=50,
        tol=0.0,
        algorithm="lloyd",
    )

    expected = cast(
        np.ndarray,
        cast(Any, KMeans)(
            n_clusters=2,
            init="random",
            n_init=2,
            max_iter=50,
            tol=0.0,
            random_state=np.random.RandomState(42),
            algorithm="lloyd",
            copy_x=True,
        ).fit_predict(mini_adata.obsm["X_pca"][:, :1]),
    )

    assert "clusters" not in mini_adata.obs
    assert copied is not None
    assert copied.obs["clusters"].astype(str).tolist() == [
        str(cluster) for cluster in expected
    ]
    assert copied.uns["clusters"]["params"]["n_pcs"] == 1
    assert copied.uns["clusters"]["params"]["init"] == "random"
    assert copied.uns["clusters"]["params"]["max_iter"] == 50
    assert copied.uns["clusters"]["params"]["tol"] == 0.0


def test_kmeans_seed_is_strictly_reproducible(mini_adata):
    first = mini_adata.copy()
    second = mini_adata.copy()

    bt.omics.tl.kmeans(
        first,
        n_clusters=2,
        representation="X_pca",
        n_init=5,
        seed=123,
        init="random",
    )
    bt.omics.tl.kmeans(
        second,
        n_clusters=2,
        representation="X_pca",
        n_init=5,
        seed=123,
        init="random",
    )

    assert (
        first.obs["kmeans"].astype(str).tolist()
        == second.obs["kmeans"].astype(str).tolist()
    )


def test_kmeans_validates_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.omics.tl.kmeans(mini_adata, n_clusters=0)

    with pytest.raises(ValueError):
        bt.omics.tl.kmeans(mini_adata, n_clusters=mini_adata.n_obs + 1)

    with pytest.raises(ValueError):
        bt.omics.tl.kmeans(mini_adata, n_pcs=0)

    with pytest.raises(ValueError):
        bt.omics.tl.kmeans(mini_adata, n_init=0)

    with pytest.raises(TypeError):
        bt.omics.tl.kmeans(mini_adata, key_added=cast(Any, object()))


def test_louvain_stores_clusters_and_metadata(mini_adata):
    _set_two_component_graph(mini_adata, "connectivities")

    result = bt.omics.tl.louvain(
        mini_adata,
        resolution=1.0,
        neighbors_key="neighbors",
        key_added="louvain",
        seed=10,
    )

    assert result is None
    clusters = mini_adata.obs["louvain"]
    assert isinstance(clusters.dtype, pd.CategoricalDtype)
    assert clusters.loc["c1"] == clusters.loc["c2"]
    assert clusters.loc["c3"] == clusters.loc["c4"]
    assert clusters.loc["c1"] != clusters.loc["c3"]
    assert mini_adata.uns["louvain"] == {
        "params": {
            "method": "louvain",
            "resolution": 1.0,
            "neighbors_key": "neighbors",
            "obsp": None,
            "weighted": True,
            "seed": 10,
        }
    }


def test_louvain_obsp_copy_and_unweighted(mini_adata):
    _set_two_component_graph(mini_adata, "custom_connectivities")

    copied = bt.omics.tl.louvain(
        mini_adata,
        resolution=1.0,
        neighbors_key=None,
        obsp="custom_connectivities",
        weighted=False,
        key_added="clusters",
        seed=10,
        copy=True,
    )

    assert "clusters" not in mini_adata.obs
    assert copied is not None
    assert copied.obs["clusters"].loc["c1"] == copied.obs["clusters"].loc["c2"]
    assert copied.obs["clusters"].loc["c3"] == copied.obs["clusters"].loc["c4"]
    assert copied.obs["clusters"].loc["c1"] != copied.obs["clusters"].loc["c3"]
    assert copied.uns["clusters"]["params"]["obsp"] == "custom_connectivities"
    assert copied.uns["clusters"]["params"]["weighted"] is False


def test_louvain_seed_is_strictly_reproducible():
    graph = csr_matrix(
        [
            [0, 1, 1, 0, 0, 0, 0, 1],
            [1, 0, 1, 1, 0, 0, 0, 0],
            [1, 1, 0, 1, 0, 0, 1, 0],
            [0, 1, 1, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 1],
            [0, 0, 1, 0, 1, 1, 0, 1],
            [1, 0, 0, 0, 0, 1, 1, 0],
        ],
        dtype=float,
    )
    first = ad.AnnData(
        X=np.zeros((8, 1)),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(8)]),
    )
    second = first.copy()
    first.obsp["connectivities"] = graph
    second.obsp["connectivities"] = graph.copy()
    first.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
    }
    second.uns["neighbors"] = first.uns["neighbors"].copy()

    bt.omics.tl.louvain(first, seed=123)
    bt.omics.tl.louvain(second, seed=123)

    assert (
        first.obs["louvain"].astype(str).tolist()
        == second.obs["louvain"].astype(str).tolist()
    )


def test_louvain_validates_graph_source_and_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.omics.tl.louvain(
            mini_adata, neighbors_key="neighbors", obsp="connectivities"
        )

    with pytest.raises(KeyError):
        bt.omics.tl.louvain(mini_adata, neighbors_key="missing")

    with pytest.raises(ValueError):
        bt.omics.tl.louvain(mini_adata, resolution=0)

    with pytest.raises(TypeError):
        bt.omics.tl.louvain(mini_adata, weighted=cast(Any, "bad"))

    with pytest.raises(TypeError):
        bt.omics.tl.louvain(mini_adata, key_added=cast(Any, object()))


def test_leiden_stores_clusters_and_metadata(mini_adata):
    _set_two_component_graph(mini_adata, "connectivities")

    result = bt.omics.tl.leiden(
        mini_adata,
        resolution=1.0,
        neighbors_key="neighbors",
        key_added="leiden",
        seed=10,
    )

    assert result is None
    clusters = mini_adata.obs["leiden"]
    assert isinstance(clusters.dtype, pd.CategoricalDtype)
    assert clusters.loc["c1"] == clusters.loc["c2"]
    assert clusters.loc["c3"] == clusters.loc["c4"]
    assert clusters.loc["c1"] != clusters.loc["c3"]
    assert mini_adata.uns["leiden"] == {
        "params": {
            "method": "leiden",
            "resolution": 1.0,
            "neighbors_key": "neighbors",
            "obsp": None,
            "directed": True,
            "weighted": True,
            "n_iterations": "auto",
            "seed": 10,
        }
    }


def test_leiden_obsp_copy_and_unweighted(mini_adata):
    _set_two_component_graph(mini_adata, "custom_connectivities")

    copied = bt.omics.tl.leiden(
        mini_adata,
        resolution=1.0,
        neighbors_key=None,
        obsp="custom_connectivities",
        weighted=False,
        directed=False,
        n_iterations=2,
        key_added="clusters",
        seed=10,
        copy=True,
    )

    assert "clusters" not in mini_adata.obs
    assert copied is not None
    assert copied.obs["clusters"].loc["c1"] == copied.obs["clusters"].loc["c2"]
    assert copied.obs["clusters"].loc["c3"] == copied.obs["clusters"].loc["c4"]
    assert copied.obs["clusters"].loc["c1"] != copied.obs["clusters"].loc["c3"]
    assert copied.uns["clusters"]["params"]["obsp"] == "custom_connectivities"
    assert copied.uns["clusters"]["params"]["directed"] is False
    assert copied.uns["clusters"]["params"]["weighted"] is False
    assert copied.uns["clusters"]["params"]["n_iterations"] == 2


def test_leiden_validates_graph_source_and_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.omics.tl.leiden(mini_adata, neighbors_key="neighbors", obsp="connectivities")

    with pytest.raises(KeyError):
        bt.omics.tl.leiden(mini_adata, neighbors_key="missing")

    with pytest.raises(ValueError):
        bt.omics.tl.leiden(mini_adata, resolution=0)

    with pytest.raises(ValueError):
        bt.omics.tl.leiden(mini_adata, n_iterations=0)

    with pytest.raises(ValueError):
        bt.omics.tl.leiden(mini_adata, n_iterations=-1)

    with pytest.raises(ValueError):
        bt.omics.tl.leiden(mini_adata, n_iterations=cast(Any, "bad"))

    with pytest.raises(TypeError):
        bt.omics.tl.leiden(mini_adata, directed=cast(Any, "bad"))
