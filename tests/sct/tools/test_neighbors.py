#!/usr/bin/env python

from typing import Any, cast

import networkx as nx
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph

import bonesistools as bt
from bonesistools.sctools.tools import _neighbors


def test_neighbors_stores_distances_connectivities_and_metadata(mini_adata):
    mini_adata.uns.clear()
    mini_adata.obsp.clear()
    close_distance = np.sqrt(0.05)
    medium_distance = np.sqrt(6.85)
    far_distance = np.sqrt(8.0)
    expected_distances = np.array(
        [
            [0.0, close_distance, far_distance, 0.0],
            [close_distance, 0.0, medium_distance, 0.0],
            [0.0, medium_distance, 0.0, close_distance],
            [0.0, far_distance, close_distance, 0.0],
        ]
    )
    expected_connectivities = np.array(
        [
            [0.0, 1.0, 0.5849658, 0.0],
            [1.0, 0.0, 0.8277489, 0.5849658],
            [0.5849658, 0.8277489, 0.0, 1.0],
            [0.0, 0.5849658, 1.0, 0.0],
        ],
    )

    result = bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        n_pcs=2,
        metric="euclidean",
        n_jobs=1,
    )

    assert result is None
    assert mini_adata.uns["neighbors"]["distances_key"] == "distances"
    assert mini_adata.uns["neighbors"]["connectivities_key"] == "connectivities"
    assert mini_adata.uns["neighbors"]["params"] == {
        "n_neighbors": 3,
        "n_pcs": 2,
        "representation": "X_pca",
        "backend": "exact",
        "metric": "euclidean",
        "connectivity_method": "fuzzy",
        "seed": 0,
    }
    assert np.allclose(
        mini_adata.obsp["distances"].toarray(),
        expected_distances,
    )
    assert np.allclose(
        mini_adata.obsp["connectivities"].toarray(),
        expected_connectivities,
    )
    assert np.all(np.asarray(mini_adata.obsp["connectivities"].sum(axis=0)) > 0)
    assert np.all(np.asarray(mini_adata.obsp["connectivities"].sum(axis=1)) > 0)


def test_neighbors_binary_connectivities_match_sklearn_spectral_graph(mini_adata):
    representation = mini_adata.obsm["X_pca"][:, :2]
    expected_connectivities = kneighbors_graph(
        representation,
        n_neighbors=3,
        mode="connectivity",
        metric="euclidean",
        include_self=True,
        n_jobs=1,
    )
    expected_connectivities = 0.5 * (
        expected_connectivities + expected_connectivities.T
    )

    bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        n_pcs=2,
        connectivity_method="binary",
        n_jobs=1,
    )

    assert mini_adata.uns["neighbors"]["params"]["connectivity_method"] == "binary"
    np.testing.assert_array_equal(
        mini_adata.obsp["connectivities"].toarray(),
        expected_connectivities.toarray(),
    )


def test_neighbors_custom_keys_and_copy(mini_adata):
    copied = bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        key_added="custom_neighbors",
        distances_key="custom_distances",
        connectivities_key="custom_connectivities",
        copy=True,
    )

    assert "custom_neighbors" not in mini_adata.uns
    assert "custom_distances" not in mini_adata.obsp
    assert "custom_connectivities" not in mini_adata.obsp
    assert copied is not None
    assert copied.uns["custom_neighbors"]["distances_key"] == "custom_distances"
    assert (
        copied.uns["custom_neighbors"]["connectivities_key"] == "custom_connectivities"
    )
    assert copied.obsp["custom_distances"].shape == (copied.n_obs, copied.n_obs)
    assert copied.obsp["custom_connectivities"].shape == (copied.n_obs, copied.n_obs)

    derived = bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        key_added="derived",
        copy=True,
    )
    assert derived is not None
    assert derived.uns["derived"]["distances_key"] == "derived_distances"
    assert derived.uns["derived"]["connectivities_key"] == "derived_connectivities"


def test_neighbors_uses_pynndescent_backend_with_sparse_transformer_output(
    mini_adata,
    monkeypatch,
):
    class FakePyNNDescentTransformer(object):
        calls = []

        def __init__(
            self,
            n_neighbors,
            metric,
            random_state,
            n_jobs,
        ):
            self.calls.append(
                {
                    "n_neighbors": n_neighbors,
                    "metric": metric,
                    "random_state": random_state,
                    "n_jobs": n_jobs,
                }
            )

        def fit_transform(self, matrix):
            assert matrix.shape == (mini_adata.n_obs, 2)
            indptr = np.arange(
                0,
                (mini_adata.n_obs + 1) * 3,
                3,
                dtype=np.int32,
            )
            indices = np.array(
                [
                    0,
                    1,
                    2,
                    1,
                    0,
                    2,
                    2,
                    3,
                    1,
                    3,
                    2,
                    0,
                ],
                dtype=np.int32,
            )
            data = np.array(
                [
                    0.0,
                    1.0,
                    2.0,
                    0.0,
                    1.0,
                    3.0,
                    0.0,
                    1.0,
                    3.0,
                    0.0,
                    1.0,
                    2.0,
                ],
                dtype=np.float32,
            )
            return csr_matrix(
                (data, indices, indptr),
                shape=(mini_adata.n_obs, mini_adata.n_obs),
            )

    monkeypatch.setattr(
        _neighbors,
        "_get_pynndescent_transformer",
        lambda: FakePyNNDescentTransformer,
    )

    bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        n_pcs=2,
        backend="pynndescent",
        metric="euclidean",
        seed=10,
        n_jobs=2,
    )

    expected_indices = np.array(
        [
            [0, 1, 2],
            [1, 0, 2],
            [2, 3, 1],
            [3, 2, 0],
        ],
        dtype=np.int32,
    )
    expected_distances = np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 3.0],
            [0.0, 1.0, 3.0],
            [0.0, 1.0, 2.0],
        ],
        dtype=np.float32,
    )
    expected_distance_mtx = _neighbors._sparse_distances_from_knn(
        expected_indices,
        expected_distances,
        mini_adata.n_obs,
    )
    expected_connectivities = _neighbors._umap_connectivities_from_knn(
        expected_indices,
        expected_distances,
        mini_adata.n_obs,
    )

    assert mini_adata.uns["neighbors"]["params"]["backend"] == "pynndescent"
    assert mini_adata.uns["neighbors"]["params"]["seed"] == 10
    assert FakePyNNDescentTransformer.calls[0]["n_neighbors"] == 3
    assert FakePyNNDescentTransformer.calls[0]["metric"] == "euclidean"
    assert FakePyNNDescentTransformer.calls[0]["n_jobs"] == 2
    assert np.array_equal(
        mini_adata.obsp["distances"].indptr,
        expected_distance_mtx.indptr,
    )
    assert np.array_equal(
        mini_adata.obsp["distances"].indices,
        expected_distance_mtx.indices,
    )
    assert np.allclose(mini_adata.obsp["distances"].data, expected_distance_mtx.data)
    assert np.array_equal(
        mini_adata.obsp["connectivities"].indptr,
        expected_connectivities.indptr,
    )
    assert np.array_equal(
        mini_adata.obsp["connectivities"].indices,
        expected_connectivities.indices,
    )
    assert np.allclose(
        mini_adata.obsp["connectivities"].data,
        expected_connectivities.data,
    )


def test_neighbors_explicit_pynndescent_requires_dependency(
    mini_adata,
    monkeypatch,
):
    def missing_pynndescent():
        raise ImportError(
            "backend='pynndescent' requires the optional dependency 'pynndescent'."
        )

    monkeypatch.setattr(
        _neighbors,
        "_get_pynndescent_transformer",
        missing_pynndescent,
    )

    with pytest.raises(ImportError):
        bt.sct.tl.neighbors(
            mini_adata,
            n_neighbors=3,
            representation="X_pca",
            n_pcs=2,
            backend="pynndescent",
        )


def test_neighbors_can_feed_shared_neighbors(mini_adata):
    mini_adata.uns.clear()
    mini_adata.obsp.clear()
    close_distance = np.sqrt(0.05)
    far_distance = np.sqrt(8.0)
    farther_distance = np.sqrt(9.25)
    expected_snn_distances = np.array(
        [
            [0.0, close_distance, far_distance, farther_distance],
            [close_distance, 0.0, 0.0, far_distance],
            [far_distance, 0.0, 0.0, close_distance],
            [farther_distance, far_distance, close_distance, 0.0],
        ]
    )
    expected_snn_connectivities = (
        np.array(
            [
                [0.0, 1.0, 1.0, 2.0],
                [1.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 1.0],
                [2.0, 1.0, 1.0, 0.0],
            ]
        )
        / 3.0
    )

    bt.sct.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        n_pcs=2,
    )
    bt.sct.tl.shared_neighbors(
        mini_adata,
        knn_key="neighbors",
        snn_key="snn_from_neighbors",
        prune_snn=0,
    )

    assert mini_adata.uns["snn_from_neighbors"]["params"]["knn_base"] == (
        "scdata.uns['neighbors']"
    )
    assert np.allclose(
        mini_adata.obsp["snn_from_neighbors_distances"].toarray(),
        expected_snn_distances,
    )
    assert np.allclose(
        mini_adata.obsp["snn_from_neighbors_connectivities"].toarray(),
        expected_snn_connectivities,
    )


def test_neighbors_validates_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.sct.tl.neighbors(mini_adata, n_neighbors=1)

    with pytest.raises(ValueError):
        bt.sct.tl.neighbors(mini_adata, n_neighbors=mini_adata.n_obs + 1)

    with pytest.raises(KeyError):
        bt.sct.tl.neighbors(mini_adata, n_neighbors=3, representation="missing")

    with pytest.raises(ValueError):
        bt.sct.tl.neighbors(mini_adata, n_neighbors=3, metric=cast(Any, "bad"))

    with pytest.raises(ValueError):
        bt.sct.tl.neighbors(mini_adata, n_neighbors=3, backend=cast(Any, "bad"))

    with pytest.raises(ValueError):
        bt.sct.tl.neighbors(
            mini_adata,
            n_neighbors=3,
            backend="pynndescent",
            connectivity_method="binary",
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'n_jobs'"):
        bt.sct.tl.neighbors(
            mini_adata,
            n_neighbors=3,
            n_jobs=cast(Any, "2"),
        )


def test_sparse_distance_neighbors_are_reordered_and_completed():
    sparse_distances = csr_matrix(
        (
            np.array([0.5, 0.0, 0.8, 0.4, 0.7, 0.0, 0.6], dtype=np.float32),
            np.array([1, 0, 2, 2, 0, 2, 1], dtype=np.int32),
            np.array([0, 3, 5, 7], dtype=np.int32),
        ),
        shape=(3, 3),
    )

    indices, distances = _neighbors._knn_arrays_from_sparse_distances(
        sparse_distances,
        n_neighbors=2,
    )

    np.testing.assert_array_equal(
        indices,
        np.array([[0, 1], [1, 2], [2, 1]], dtype=np.int32),
    )
    np.testing.assert_allclose(
        distances,
        np.array([[0.0, 0.5], [0.0, 0.4], [0.0, 0.6]], dtype=np.float32),
    )


def test_sparse_distance_neighbors_reject_short_rows():
    sparse_distances = csr_matrix(
        (
            np.array([0.0, 0.2, 0.0], dtype=np.float32),
            np.array([0, 1, 1], dtype=np.int32),
            np.array([0, 2, 3], dtype=np.int32),
        ),
        shape=(2, 2),
    )

    with pytest.raises(ValueError, match="fewer neighbors than expected"):
        _neighbors._knn_arrays_from_sparse_distances(
            sparse_distances,
            n_neighbors=2,
        )


def test_knn_graph_can_use_observation_names(mini_adata):
    graph = bt.sct.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        representation="X_pca",
        create_using=nx.DiGraph,
        index_or_name="name",
    )

    assert set(graph.nodes) == set(mini_adata.obs_names)
    assert graph.number_of_edges() == mini_adata.n_obs
    assert {
        (source, target): data["distance"]
        for source, target, data in graph.edges(data=True)
    } == pytest.approx(
        {
            ("c1", "c2"): np.sqrt(0.06),
            ("c2", "c1"): np.sqrt(0.06),
            ("c3", "c4"): np.sqrt(0.06),
            ("c4", "c3"): np.sqrt(0.06),
        }
    )

    index_graph = bt.sct.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        representation="X_pca",
        create_using=nx.DiGraph,
        index_or_name="index",
    )

    assert set(index_graph.nodes) == set(range(mini_adata.n_obs))


def test_knn_graph_forwards_metric_kwargs(mini_adata):
    graph = bt.sct.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        representation="X_pca",
        n_components=2,
        metric="minkowski",
        p=1,
        create_using=nx.DiGraph,
        index_or_name="name",
    )

    assert {
        (source, target): data["distance"]
        for source, target, data in graph.edges(data=True)
    } == pytest.approx(
        {
            ("c1", "c2"): 0.3,
            ("c2", "c1"): 0.3,
            ("c3", "c4"): 0.3,
            ("c4", "c3"): 0.3,
        }
    )


def test_kneighbors_graph_deprecated_alias(mini_adata):
    with pytest.warns(FutureWarning, match="kneighbors_graph"):
        graph = bt.sct.tl.kneighbors_graph(
            mini_adata,
            n_neighbors=1,
            representation="X_pca",
            create_using=nx.DiGraph,
            index_or_name="name",
        )

    assert set(graph.nodes) == set(mini_adata.obs_names)


def test_knn_graph_rejects_invalid_node_label_mode(mini_adata):
    with pytest.raises(TypeError, match="missing required argument: 'n_neighbors'"):
        bt.sct.tl.knn_graph(mini_adata)

    with pytest.raises(ValueError, match="invalid argument value for 'index_or_name'"):
        bt.sct.tl.knn_graph(
            mini_adata,
            n_neighbors=1,
            representation="X_pca",
            index_or_name=cast(Any, "bad"),
        )


def test_shared_neighbors_stores_distances_connectivities_and_metadata(
    mini_adata,
    expected_mini_pca2_distances,
    expected_mini_snn_connectivities,
):
    result = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=0,
        snn_key="snn",
        copy=True,
    )

    assert "snn" not in mini_adata.uns
    assert result.uns["snn"]["distances_key"] == "snn_distances"
    assert result.uns["snn"]["connectivities_key"] == "snn_connectivities"

    expected_distances = np.zeros_like(expected_mini_pca2_distances)
    expected_distances[0, 3] = expected_mini_pca2_distances[0, 3]
    expected_distances[1, 2] = expected_mini_pca2_distances[1, 2]
    expected_distances[2, 1] = expected_mini_pca2_distances[2, 1]
    expected_distances[3, 0] = expected_mini_pca2_distances[3, 0]

    assert np.allclose(
        result.obsp["snn_distances"].toarray(),
        expected_distances,
    )
    assert np.allclose(
        result.obsp["snn_connectivities"].toarray(),
        expected_mini_snn_connectivities / 3,
    )


def test_shared_neighbors_validates_source_graph_and_pruning(mini_adata):
    with pytest.raises(KeyError, match="neighborhood graph not found"):
        bt.sct.tl.shared_neighbors(mini_adata, knn_key="missing")

    with pytest.raises(ValueError, match="invalid argument value for 'prune_snn'"):
        bt.sct.tl.shared_neighbors(mini_adata, prune_snn=-1)
