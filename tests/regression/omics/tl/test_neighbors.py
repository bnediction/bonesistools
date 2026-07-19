#!/usr/bin/env python

import builtins
from typing import Any, cast

import anndata as ad
import networkx as nx
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph

import bonesistools as bt
from bonesistools.omics.tools import _neighbors


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

    result = bt.omics.tl.neighbors(
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
    expected_connectivities = cast(
        csr_matrix,
        kneighbors_graph(
            representation,
            n_neighbors=3,
            mode="connectivity",
            metric="euclidean",
            include_self=True,
            n_jobs=1,
        ),
    )
    expected_connectivities = cast(
        csr_matrix,
        0.5 * (expected_connectivities + expected_connectivities.T),
    )

    bt.omics.tl.neighbors(
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
    copied = bt.omics.tl.neighbors(
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

    derived = bt.omics.tl.neighbors(
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

    bt.omics.tl.neighbors(
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
        bt.omics.tl.neighbors(
            mini_adata,
            n_neighbors=3,
            representation="X_pca",
            n_pcs=2,
            backend="pynndescent",
        )


def test_neighbors_can_feed_shared_neighbors(mini_adata):
    mini_adata.uns.clear()
    mini_adata.obsp.clear()
    expected_snn_distances = np.array(
        [
            [0.0, 2 / 3, 2 / 3, 0.0],
            [2 / 3, 0.0, 0.0, 2 / 3],
            [2 / 3, 0.0, 0.0, 2 / 3],
            [0.0, 2 / 3, 2 / 3, 0.0],
        ]
    )
    expected_snn_connectivities = np.array(
        [
            [0.0, 1 / 3, 1 / 3, 1.0],
            [1 / 3, 0.0, 0.0, 1 / 3],
            [1 / 3, 0.0, 0.0, 1 / 3],
            [1.0, 1 / 3, 1 / 3, 0.0],
        ]
    )

    bt.omics.tl.neighbors(
        mini_adata,
        n_neighbors=3,
        representation="X_pca",
        n_pcs=2,
    )
    bt.omics.tl.shared_neighbors(
        mini_adata,
        neighbors_key="neighbors",
        key_added="snn_from_neighbors",
        prune=0,
    )

    assert mini_adata.uns["snn_from_neighbors"]["params"]["knn_base"] == (
        "scdata.uns['neighbors']"
    )
    assert mini_adata.uns["snn_from_neighbors"]["params"]["metric"] == "jaccard"
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
        bt.omics.tl.neighbors(mini_adata, n_neighbors=1)

    with pytest.raises(ValueError):
        bt.omics.tl.neighbors(mini_adata, n_neighbors=mini_adata.n_obs + 1)

    with pytest.raises(KeyError):
        bt.omics.tl.neighbors(mini_adata, n_neighbors=3, representation="missing")

    with pytest.raises(ValueError):
        bt.omics.tl.neighbors(mini_adata, n_neighbors=3, metric=cast(Any, "bad"))

    with pytest.raises(ValueError):
        bt.omics.tl.neighbors(mini_adata, n_neighbors=3, backend=cast(Any, "bad"))

    with pytest.raises(ValueError):
        bt.omics.tl.neighbors(
            mini_adata,
            n_neighbors=3,
            backend="pynndescent",
            connectivity_method="binary",
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'n_jobs'"):
        bt.omics.tl.neighbors(
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


def test_sparse_distance_neighbors_move_self_neighbor_before_equal_distances():
    sparse_distances = csr_matrix(
        (
            np.array([0.0, 0.0, 0.2, 0.0, 0.2], dtype=np.float32),
            np.array([1, 0, 0, 1, 0], dtype=np.int32),
            np.array([0, 3, 5], dtype=np.int32),
        ),
        shape=(2, 2),
    )

    indices, distances = _neighbors._knn_arrays_from_sparse_distances(
        sparse_distances,
        n_neighbors=2,
    )

    np.testing.assert_array_equal(
        indices,
        np.array([[0, 1], [1, 0]], dtype=np.int32),
    )
    np.testing.assert_allclose(
        distances,
        np.array([[0.0, 0.0], [0.0, 0.2]], dtype=np.float32),
    )


def test_pynndescent_backend_reports_missing_optional_dependency(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pynndescent":
            raise ImportError("missing pynndescent")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ImportError, match="backend='pynndescent' requires"):
        _neighbors._get_pynndescent_transformer()


def test_knn_arrays_are_sorted_canonically():
    indices = np.array(
        [
            [3, 0, 1, 2],
            [2, 0, 1, 3],
        ],
        dtype=np.int32,
    )
    distances = np.array(
        [
            [0.5, 1e-7, 0.5, 0.2],
            [0.7, 0.1, 1e-7, 0.1],
        ],
        dtype=np.float32,
    )

    sorted_indices, sorted_distances = _neighbors._sort_knn_arrays(
        knn_indices=indices,
        knn_distances=distances,
    )

    np.testing.assert_array_equal(
        sorted_indices,
        np.array(
            [
                [0, 2, 1, 3],
                [1, 0, 3, 2],
            ],
            dtype=np.int32,
        ),
    )
    np.testing.assert_allclose(
        sorted_distances,
        np.array(
            [
                [0.0, 0.2, 0.5, 0.5],
                [0.0, 0.1, 0.1, 0.7],
            ],
            dtype=np.float32,
        ),
    )


def test_umap_connectivities_ignore_numerical_self_distances():
    indices = np.array(
        [
            [0, 1, 2],
            [1, 0, 2],
            [2, 1, 0],
        ],
        dtype=np.int32,
    )
    distances = np.array(
        [
            [0.0, 0.2, 0.5],
            [0.0, 0.3, 0.6],
            [0.0, 0.4, 0.7],
        ],
        dtype=np.float32,
    )
    noisy_distances = distances.copy()
    noisy_distances[:, 0] = np.array([1e-7, 2e-7, 3e-7], dtype=np.float32)

    expected = _neighbors._umap_connectivities_from_knn(indices, distances, n_obs=3)
    actual = _neighbors._umap_connectivities_from_knn(
        indices,
        noisy_distances,
        n_obs=3,
    )

    np.testing.assert_array_equal(actual.toarray(), expected.toarray())


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
    graph = bt.omics.tl.knn_graph(
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

    index_graph = bt.omics.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        representation="X_pca",
        create_using=nx.DiGraph,
        index_or_name="index",
    )

    assert set(index_graph.nodes) == set(range(mini_adata.n_obs))


def test_knn_graph_forwards_metric_kwargs(mini_adata):
    graph = bt.omics.tl.knn_graph(
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


def test_knn_graph_supports_legacy_networkx_sparse_matrix_constructor(
    mini_adata,
    monkeypatch,
):
    from_scipy_sparse_array = _neighbors.nx.from_scipy_sparse_array
    calls = []

    def from_scipy_sparse_matrix(matrix, create_using, edge_attribute):
        calls.append(edge_attribute)
        return from_scipy_sparse_array(
            matrix,
            create_using=create_using,
            edge_attribute=edge_attribute,
        )

    monkeypatch.delattr(_neighbors.nx, "from_scipy_sparse_array")
    monkeypatch.setattr(
        _neighbors.nx,
        "from_scipy_sparse_matrix",
        from_scipy_sparse_matrix,
        raising=False,
    )

    graph = bt.omics.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        representation="X_pca",
        create_using=nx.DiGraph,
        index_or_name="index",
    )

    assert calls == ["distance"]
    assert graph.number_of_edges() == mini_adata.n_obs


def test_kneighbors_graph_deprecated_alias(mini_adata):
    with pytest.warns(FutureWarning, match="kneighbors_graph"):
        graph = getattr(bt.omics.tl, "kneighbors_graph")(
            mini_adata,
            n_neighbors=1,
            representation="X_pca",
            create_using=nx.DiGraph,
            index_or_name="name",
        )

    assert set(graph.nodes) == set(mini_adata.obs_names)


def test_knn_graph_rejects_invalid_node_label_mode(mini_adata):
    with pytest.raises(TypeError, match="missing required argument: 'n_neighbors'"):
        bt.omics.tl.knn_graph(mini_adata)

    with pytest.raises(ValueError, match="invalid argument value for 'index_or_name'"):
        bt.omics.tl.knn_graph(
            mini_adata,
            n_neighbors=1,
            representation="X_pca",
            index_or_name=cast(Any, "bad"),
        )


def test_shared_neighbors_stores_distances_connectivities_and_metadata(
    mini_adata,
    expected_mini_snn_connectivities,
):
    result = bt.omics.tl.shared_neighbors(
        mini_adata,
        prune=0,
        key_added="snn",
        copy=True,
    )

    assert "snn" not in mini_adata.uns
    assert result.uns["snn"]["distances_key"] == "snn_distances"
    assert result.uns["snn"]["connectivities_key"] == "snn_connectivities"
    assert result.uns["snn"]["params"]["prune"] == 0
    assert "prune_snn" not in result.uns["snn"]["params"]

    distances = result.obsp["snn_distances"]
    assert distances.nnz == 4
    assert np.all(distances.data == 0.0)
    assert np.allclose(
        result.obsp["snn_connectivities"].toarray(),
        expected_mini_snn_connectivities / 2,
    )


@pytest.mark.parametrize(
    ("metric", "expected"),
    (
        ("jaccard", 2 / 3),
        ("overlap", 1.0),
        ("binary-cosine", 2 / np.sqrt(6)),
        ("weighted-jaccard", 4 / 9),
        ("ranked-jaccard", 9 / 11),
    ),
)
def test_shared_neighbors_metrics_are_normalized(metric, expected):
    adata = ad.AnnData(X=np.zeros((5, 1), dtype=np.float32))
    adata.obsp["distances"] = csr_matrix(
        [
            [0.0, 0.0, 0.1, 0.4, 0.0],
            [0.0, 0.0, 0.2, 0.3, 0.5],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    adata.obsp["connectivities"] = csr_matrix(
        [
            [0.0, 0.0, 0.9, 0.4, 0.0],
            [0.0, 0.0, 0.6, 0.2, 0.5],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {"n_neighbors": 4},
    }

    bt.omics.tl.shared_neighbors(
        adata,
        prune=0,
        metric=metric,
    )

    connectivities = cast(csr_matrix, adata.obsp["shared_neighbors_connectivities"])
    distances = cast(csr_matrix, adata.obsp["shared_neighbors_distances"])
    assert connectivities.nnz == 2
    assert distances.nnz == 2
    assert connectivities[0, 1] == pytest.approx(expected)
    assert connectivities[1, 0] == pytest.approx(expected)
    assert distances[0, 1] == pytest.approx(1 - expected)
    assert np.all((connectivities.data >= 0) & (connectivities.data <= 1))
    assert adata.uns["shared_neighbors"]["params"]["metric"] == metric


def test_shared_neighbors_metrics_match_dense_definitions():
    neighborhoods = (
        (1, 2, 4),
        (2, 3),
        (0, 3, 4),
        (0, 1, 4),
        (1,),
    )
    n_obs = len(neighborhoods)
    source_distances = np.zeros((n_obs, n_obs), dtype=np.float32)
    source_weights = np.zeros((n_obs, n_obs), dtype=np.float32)
    ranks = []
    for row, neighbors in enumerate(neighborhoods):
        row_ranks = {}
        for rank, neighbor in enumerate(neighbors, start=1):
            source_distances[row, neighbor] = rank / 10
            source_weights[row, neighbor] = (1 + ((2 * row + neighbor) % 7)) / 10
            row_ranks[neighbor] = rank
        ranks.append(row_ranks)

    adata = ad.AnnData(X=np.zeros((n_obs, 1), dtype=np.float32))
    adata.obsp["distances"] = csr_matrix(source_distances)
    adata.obsp["connectivities"] = csr_matrix(source_weights)
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {"n_neighbors": 4},
    }

    maximum_ranked_score = sum(1 / (2 * rank) for rank in range(1, 4))
    for metric in (
        "jaccard",
        "overlap",
        "binary-cosine",
        "weighted-jaccard",
        "ranked-jaccard",
    ):
        expected = np.zeros((n_obs, n_obs), dtype=np.float32)
        for left in range(n_obs):
            left_neighbors = set(neighborhoods[left])
            for right in range(left + 1, n_obs):
                right_neighbors = set(neighborhoods[right])
                shared = left_neighbors & right_neighbors
                if not shared:
                    continue
                if metric == "jaccard":
                    value = len(shared) / len(left_neighbors | right_neighbors)
                elif metric == "overlap":
                    value = len(shared) / min(
                        len(left_neighbors),
                        len(right_neighbors),
                    )
                elif metric == "binary-cosine":
                    value = len(shared) / np.sqrt(
                        len(left_neighbors) * len(right_neighbors)
                    )
                elif metric == "weighted-jaccard":
                    intersection = sum(
                        min(
                            source_weights[left, neighbor],
                            source_weights[right, neighbor],
                        )
                        for neighbor in shared
                    )
                    union = sum(
                        max(
                            source_weights[left, neighbor],
                            source_weights[right, neighbor],
                        )
                        for neighbor in left_neighbors | right_neighbors
                    )
                    value = intersection / union
                else:
                    value = (
                        sum(
                            1 / (ranks[left][neighbor] + ranks[right][neighbor])
                            for neighbor in shared
                        )
                        / maximum_ranked_score
                    )
                expected[left, right] = value
                expected[right, left] = value

        bt.omics.tl.shared_neighbors(
            adata,
            prune=0,
            metric=cast(Any, metric),
        )

        connectivities = cast(
            csr_matrix,
            adata.obsp["shared_neighbors_connectivities"],
        )
        distances = cast(csr_matrix, adata.obsp["shared_neighbors_distances"])
        expected_distances = np.where(expected > 0, 1 - expected, 0)
        np.testing.assert_allclose(connectivities.toarray(), expected, rtol=1e-6)
        np.testing.assert_allclose(
            distances.toarray(),
            expected_distances,
            rtol=1e-6,
        )
        assert distances.nnz == np.count_nonzero(expected)


def test_shared_neighbors_supports_more_than_127_neighbors():
    n_obs = 132
    source_distances = csr_matrix(
        np.ones((n_obs, n_obs), dtype=np.float32)
        - np.eye(n_obs, dtype=np.float32)
    )
    adata = ad.AnnData(X=np.zeros((n_obs, 1), dtype=np.float32))
    adata.obsp["distances"] = source_distances
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "params": {"n_neighbors": n_obs},
    }

    bt.omics.tl.shared_neighbors(adata, prune=0)

    connectivities = cast(
        csr_matrix,
        adata.obsp["shared_neighbors_connectivities"],
    )
    assert connectivities[0, 1] == pytest.approx((n_obs - 2) / n_obs)


def test_shared_neighbors_deprecated_key_arguments_warn_and_delegate(mini_adata):
    deprecated_shared_neighbors = cast(Any, bt.omics.tl.shared_neighbors)

    with pytest.warns(FutureWarning) as warning_records:
        result = deprecated_shared_neighbors(
            mini_adata,
            knn_key="neighbors",
            snn_key="legacy_snn",
            prune=0,
            copy=True,
        )

    messages = tuple(str(record.message) for record in warning_records)
    assert any("`knn_key` is deprecated" in message for message in messages)
    assert any("`snn_key` is deprecated" in message for message in messages)
    assert "legacy_snn" in result.uns


def test_shared_neighbors_deprecated_and_current_keys_are_mutually_exclusive(
    mini_adata,
):
    deprecated_shared_neighbors = cast(Any, bt.omics.tl.shared_neighbors)

    with pytest.warns(FutureWarning, match="`knn_key` is deprecated"):
        with pytest.raises(TypeError, match="either 'knn_key' or 'neighbors_key'"):
            deprecated_shared_neighbors(
                mini_adata,
                neighbors_key="neighbors",
                knn_key="neighbors",
            )


def test_shared_neighbors_does_not_accept_deprecated_prune_name(mini_adata):
    with pytest.raises(TypeError, match="unexpected keyword argument 'prune_snn'"):
        cast(Any, bt.omics.tl.shared_neighbors)(mini_adata, prune_snn=0)


def test_shared_neighbors_options_are_keyword_only(mini_adata):
    with pytest.raises(TypeError):
        cast(Any, bt.omics.tl.shared_neighbors)(mini_adata, "neighbors")


def test_shared_neighbors_validates_source_graph_and_pruning(mini_adata):
    with pytest.raises(KeyError, match="neighborhood graph not found"):
        bt.omics.tl.shared_neighbors(mini_adata, neighbors_key="missing")

    with pytest.raises(ValueError, match="invalid argument value for 'prune'"):
        bt.omics.tl.shared_neighbors(mini_adata, prune=-1)

    with pytest.raises(ValueError, match="expected integer"):
        bt.omics.tl.shared_neighbors(mini_adata, prune=1.5)

    with pytest.raises(ValueError, match="invalid argument value for 'metric'"):
        bt.omics.tl.shared_neighbors(mini_adata, metric=cast(Any, "euclidean"))

    del mini_adata.obsp["connectivities"]
    with pytest.raises(KeyError, match="requires a source connectivity graph"):
        bt.omics.tl.shared_neighbors(mini_adata, metric="weighted-jaccard")
