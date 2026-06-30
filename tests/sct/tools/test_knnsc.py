#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import pytest

import bonesistools as bt
from tests.sct.toy_data import make_nestorowa_hvg_adata

ADATA = make_nestorowa_hvg_adata()


def _set_manual_shortest_path_lengths(estimator):
    estimator._obs = pd.Series(
        pd.Categorical(["A", "A", "B", "B"]),
        index=["c1", "c2", "c3", "c4"],
        name="cluster",
    )
    estimator._shortest_path_lengths_df = pd.DataFrame(
        {
            "A": [0.1, 0.8, 2.0, 1.5],
            "B": [2.0, 1.5, 0.2, 0.9],
        },
        index=estimator.obs.index,
    )
    estimator._cluster_counts = estimator.obs.value_counts()
    estimator._min_cluster_size = 1


def _grid_cloud(center, radius=0.2, steps=5):
    x_center, y_center = center
    return np.array(
        [
            (x, y)
            for x in np.linspace(x_center - radius, x_center + radius, steps)
            for y in np.linspace(y_center - radius, y_center + radius, steps)
        ]
    )


def _make_two_cluster_adata(cloud_a, cloud_b, names_a, names_b):
    coordinates = np.vstack([cloud_a, cloud_b])
    obs_names = names_a + names_b
    adata = ad.AnnData(
        X=np.zeros((coordinates.shape[0], 1)),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * len(cloud_a) + ["B"] * len(cloud_b))},
            index=obs_names,
        ),
    )
    adata.obsm["X_toy"] = coordinates
    return adata, coordinates, obs_names


def test_knn_graph_returns_graph_with_named_nodes():

    adata = ADATA.copy()

    adata.obsm["X_test"] = np.asarray(cast(Any, adata.X))[:, :5].copy()

    graph = bt.sct.tl.knn_graph(
        adata,
        representation="X_test",
        n_components=None,
        n_neighbors=5,
        index_or_name="name",
    )

    assert isinstance(graph, nx.DiGraph)
    assert graph.number_of_nodes() == adata.n_obs
    assert graph.number_of_edges() > 0
    assert set(graph.nodes) == set(adata.obs_names)


def test_knnsc_fits_and_keeps_deprecated_api_compatible():

    adata = ADATA[:200, :].copy()

    adata.obsm["X_test"] = np.asarray(cast(Any, adata.X))[:, :5].copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    with pytest.warns(FutureWarning, match="__init__"):
        estimator = bt.sct.tl.KNNSC(
            representation="X_test",
            n_components=5,
            n_neighbors=5,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)

        with pytest.warns(FutureWarning, match="obs"):
            estimator.fit(adata, obs="cluster", n_jobs=1)

    with pytest.warns(FutureWarning, match="compute_shortest_path_lengths"):
        estimator.compute_shortest_path_lengths(n_jobs=1)

    with pytest.warns(FutureWarning, match="shortest_path_lengths"):
        estimator.shortest_path_lengths(n_jobs=1)

    with pytest.warns(FutureWarning, match="knnbs"):
        subclusters = estimator.knnbs(size=5, key="knnbs")

    assert hasattr(estimator, "knn_graph")
    assert hasattr(estimator, "shortest_path_lengths_df")

    assert isinstance(estimator.knn_graph, nx.Graph)
    assert isinstance(estimator.shortest_path_lengths_df, pd.DataFrame)
    assert isinstance(subclusters, pd.Series)

    assert subclusters.name == "knnbs"
    assert subclusters.index.equals(adata.obs.index)
    assert subclusters.notna().sum() > 0
    assert set(subclusters.dropna().cat.categories) == set(
        adata.obs["cluster"].cat.categories
    )


def test_knnsc_fit_reconnects_disconnected_neighbor_graph():
    cloud_a = _grid_cloud(center=(-5.0, -5.0))
    cloud_b = _grid_cloud(center=(5.0, 5.0))
    adata, _, _ = _make_two_cluster_adata(
        cloud_a=cloud_a,
        cloud_b=cloud_b,
        names_a=[f"a{i}" for i in range(len(cloud_a))],
        names_b=[f"b{i}" for i in range(len(cloud_b))],
    )

    estimator = bt.sct.tl.KNNSC()

    with pytest.warns(UserWarning, match="not weakly connected"):
        estimator.fit(
            adata,
            cluster_key="cluster",
            representation="X_toy",
            n_components=None,
            n_neighbors=20,
            n_jobs=1,
        )

    cell_nodes = set(adata.obs_names)
    cell_clusters = adata.obs["cluster"].to_dict()
    bridge_distances = [
        data["distance"]
        for source, target, data in estimator.knn_graph.edges(data=True)
        if (
            source in cell_nodes
            and target in cell_nodes
            and cell_clusters[source] != cell_clusters[target]
        )
    ]

    assert nx.is_connected(estimator.knn_graph)
    assert len(bridge_distances) == 1
    assert min(bridge_distances) > 13.0


def test_knnsc_fit_uses_minimum_cluster_size_for_barycenter_sources():
    coordinates = np.array(
        [
            [0.0, 0.0],
            [0.1, 0.0],
            [0.2, 0.0],
            [1.0, 0.0],
            [1.1, 0.0],
            [2.0, 0.0],
        ]
    )
    adata = ad.AnnData(
        X=np.zeros((coordinates.shape[0], 1)),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "A", "B", "B", "C"],
                    categories=["A", "B", "C"],
                )
            },
            index=[f"cell_{i}" for i in range(coordinates.shape[0])],
        ),
    )
    adata.obsm["X_toy"] = coordinates
    estimator = bt.sct.tl.KNNSC()

    estimator.fit(
        adata,
        cluster_key="cluster",
        representation="X_toy",
        n_components=None,
        n_neighbors=2,
        min_cluster_size=2,
        n_jobs=1,
    )

    assert set(adata.obs_names) <= set(estimator.knn_graph.nodes)
    assert "A" in estimator.knn_graph.nodes
    assert "B" in estimator.knn_graph.nodes
    assert "C" not in estimator.knn_graph.nodes
    assert list(estimator.obs.cat.categories) == ["A", "B"]
    assert pd.isna(estimator.obs.loc["cell_5"])
    assert list(estimator.shortest_path_lengths_df.columns) == ["A", "B"]


def test_knnsc_rejects_requested_clusters_below_minimum_size():
    coordinates = np.array(
        [
            [0.0, 0.0],
            [0.1, 0.0],
            [0.2, 0.0],
            [1.0, 0.0],
        ]
    )
    adata = ad.AnnData(
        X=np.zeros((coordinates.shape[0], 1)),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "A", "B"],
                    categories=["A", "B"],
                )
            },
            index=[f"cell_{i}" for i in range(coordinates.shape[0])],
        ),
    )
    adata.obsm["X_toy"] = coordinates
    estimator = bt.sct.tl.KNNSC()

    estimator.fit(
        adata,
        cluster_key="cluster",
        representation="X_toy",
        n_components=None,
        n_neighbors=2,
        min_cluster_size=2,
        n_jobs=1,
    )

    with pytest.raises(ValueError, match="B.*size=1.*min_cluster_size=2"):
        estimator.predict(peripheral_clusters=["B"])


def test_knnsc_select_central_cells_recovers_nearest_single_cluster_cells():
    central_points = np.array(
        [
            [0.0, 0.0],
            [0.05, 0.0],
            [-0.05, 0.0],
            [0.0, 0.05],
            [0.0, -0.05],
        ]
    )
    remote_points = np.array(
        [
            [8.0, 0.0],
            [-8.0, 0.0],
            [0.0, 8.0],
            [0.0, -8.0],
            [6.0, 6.0],
            [-6.0, -6.0],
            [6.0, -6.0],
            [-6.0, 6.0],
            [10.0, 0.0],
            [-10.0, 0.0],
        ]
    )
    coordinates = np.vstack([central_points, remote_points])
    obs_names = [f"central_{i}" for i in range(len(central_points))] + [
        f"remote_{i}" for i in range(len(remote_points))
    ]
    adata = ad.AnnData(
        X=np.zeros((coordinates.shape[0], 1)),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * len(coordinates))},
            index=obs_names,
        ),
    )
    adata.obsm["X_toy"] = coordinates
    estimator = bt.sct.tl.KNNSC()

    estimator.fit(
        adata,
        cluster_key="cluster",
        representation="X_toy",
        n_components=None,
        n_neighbors=10,
        n_jobs=1,
    )

    central = estimator.predict(
        subcluster_size=5,
        central_clusters=["A"],
    )

    assert central.name == "knnsc"
    assert set(central[central == "A"].index) == {
        f"central_{i}" for i in range(len(central_points))
    }


def test_knnsc_select_peripheral_cells_recovers_remote_cells():
    core_a = _grid_cloud(center=(-5.0, -5.0), radius=0.1, steps=4)
    remote_a = _grid_cloud(center=(-6.0, -6.0), radius=0.04, steps=3)[[0, 2, 4, 6, 8]]
    core_b = _grid_cloud(center=(5.0, 5.0), radius=0.1, steps=4)
    remote_b = _grid_cloud(center=(6.0, 6.0), radius=0.04, steps=3)[[0, 2, 4, 6, 8]]
    cloud_a = np.vstack([core_a, remote_a])
    cloud_b = np.vstack([core_b, remote_b])
    names_a = [f"a_core_{i}" for i in range(len(core_a))] + [
        f"a_remote_{i}" for i in range(len(remote_a))
    ]
    names_b = [f"b_core_{i}" for i in range(len(core_b))] + [
        f"b_remote_{i}" for i in range(len(remote_b))
    ]
    adata, _, _ = _make_two_cluster_adata(
        cloud_a=cloud_a,
        cloud_b=cloud_b,
        names_a=names_a,
        names_b=names_b,
    )
    estimator = bt.sct.tl.KNNSC()

    with pytest.warns(UserWarning, match="not weakly connected"):
        estimator.fit(
            adata,
            cluster_key="cluster",
            representation="X_toy",
            n_components=None,
            n_neighbors=16,
            n_jobs=1,
        )

    peripheral = estimator.select_peripheral_cells(subcluster_size=5)

    assert peripheral.name == "knnsc"
    assert set(peripheral[peripheral == "A"].index) == {
        f"a_remote_{i}" for i in range(5)
    }
    assert set(peripheral[peripheral == "B"].index) == {
        f"b_remote_{i}" for i in range(5)
    }


def test_knnsc_validates_init_and_fit_arguments_and_repr(mini_adata):
    estimator = bt.sct.tl.KNNSC()

    assert repr(estimator) == "KNNSC()"

    with pytest.raises(ValueError, match="invalid argument value for 'n_neighbors'"):
        estimator.fit(mini_adata, cluster_key="cluster", n_neighbors=0)

    with pytest.raises(ValueError, match="expected integer"):
        estimator.fit(mini_adata, cluster_key="cluster", n_neighbors=cast(Any, 1.5))

    with pytest.raises(TypeError, match="unsupported argument type for 'n_neighbors'"):
        estimator.fit(mini_adata, cluster_key="cluster", n_neighbors=cast(Any, "5"))

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            n_components=0,
            n_neighbors=5,
        )

    with pytest.raises(ValueError, match="expected integer"):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            n_components=cast(Any, 1.2),
            n_neighbors=5,
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'n_components'"):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            n_components=cast(Any, "2"),
            n_neighbors=5,
        )

    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'representation'",
    ):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            representation=cast(Any, object()),
            n_neighbors=5,
        )

    with pytest.raises(ValueError, match="invalid argument value for 'metric'"):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            representation=None,
            n_components=None,
            n_neighbors=5,
            metric=cast(Any, "not-a-metric"),
        )

    with pytest.raises(TypeError, match="missing required argument: 'n_neighbors'"):
        estimator.fit(mini_adata, cluster_key="cluster")

    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'metric_kwargs'",
    ):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            n_neighbors=5,
            metric_kwargs=cast(Any, "bad"),
        )

    non_categorical = mini_adata.copy()
    non_categorical.obs["cluster"] = non_categorical.obs["cluster"].astype(str)

    with pytest.raises(AttributeError, match="has no attribute 'cat'"):
        estimator.fit(
            non_categorical,
            cluster_key="cluster",
            n_neighbors=3,
        )

    with pytest.raises(ValueError, match="no non-empty clusters"):
        estimator.fit(
            mini_adata,
            cluster_key="cluster",
            n_neighbors=3,
            min_cluster_size=mini_adata.n_obs + 1,
        )

    with pytest.warns(FutureWarning, match="__init__"):
        deprecated_estimator = bt.sct.tl.KNNSC(
            n_components=cast(Any, 2.0),
            n_neighbors=cast(Any, 5.0),
            metric="cosine",
        )

    assert "KNNSC(representation=X_pca, n_components=2, n_neighbors=5" in repr(
        deprecated_estimator
    )

    with pytest.warns(FutureWarning, match="metric_kwds"):
        assert deprecated_estimator.metric_kwds == {}
    with pytest.warns(FutureWarning, match="metric_kwds"):
        deprecated_estimator.metric_kwds = {"p": 1}
    assert deprecated_estimator.metric_kwargs == {"p": 1}

    with pytest.warns(FutureWarning, match="obs"):
        with pytest.raises(TypeError, match="received both 'cluster_key'"):
            estimator.fit(
                mini_adata,
                cluster_key="cluster",
                n_neighbors=3,
                obs="cluster",
            )

    with pytest.raises(TypeError, match="missing required argument: 'cluster_key'"):
        estimator.fit(mini_adata, n_neighbors=3)


def test_knnsc_deprecated_graph_property_and_class_alias():
    estimator = bt.sct.tl.KNNSC()
    graph = nx.Graph()
    graph.add_edge("c1", "c2", distance=1.0)

    with pytest.warns(FutureWarning, match="kneighbors_graph"):
        estimator.kneighbors_graph = graph

    with pytest.warns(FutureWarning, match="kneighbors_graph"):
        returned = estimator.kneighbors_graph

    assert returned is graph

    with pytest.warns(FutureWarning) as warning_records:
        deprecated = bt.sct.tl.Knnbs(n_neighbors=3)

    assert isinstance(deprecated, bt.sct.tl.KNNSC)
    assert any("Knnbs" in str(record.message) for record in warning_records)


def test_knnsc_fit_stores_params_and_updates_repr(mini_adata):
    estimator = bt.sct.tl.KNNSC()

    assert repr(estimator) == "KNNSC()"
    assert not hasattr(estimator, "params_")
    with pytest.raises(AttributeError, match="has not been fitted"):
        estimator.params_

    estimator.fit(
        mini_adata,
        cluster_key="cluster",
        representation="X_pca",
        n_components=2,
        n_neighbors=3,
        metric="euclidean",
        metric_kwargs={},
        min_cluster_size=1,
        method="dijkstra",
    )

    assert estimator.params_ == {
        "cluster_key": "cluster",
        "representation": "X_pca",
        "n_components": 2,
        "n_neighbors": 3,
        "metric": "euclidean",
        "metric_kwargs": {},
        "min_cluster_size": 1,
        "method": "dijkstra",
    }
    assert (
        "KNNSC(cluster_key=cluster, representation=X_pca, "
        "n_components=2, n_neighbors=3"
    ) in repr(estimator)

    read_only_attributes = [
        ("knn_graph", nx.Graph()),
        ("shortest_path_lengths_df", pd.DataFrame()),
        ("cluster_counts", pd.Series(dtype=int)),
        ("obs", pd.Series(dtype="category")),
        ("params_", {}),
    ]
    for name, value in read_only_attributes:
        with pytest.raises(AttributeError):
            setattr(estimator, name, value)


def test_knnsc_deprecated_selection_modes_and_overlap_error():
    estimator = bt.sct.tl.KNNSC()
    _set_manual_shortest_path_lengths(estimator)

    with pytest.warns(FutureWarning, match="find_closest"):
        closest = estimator.find_closest_cells_to_self_barycenter(
            size=1,
            key="closest",
        )
    with pytest.warns(FutureWarning, match="find_furthest"):
        furthest = estimator.find_furthest_cells_to_other_barycenters(
            size=1,
            key="furthest",
        )
    with pytest.warns(FutureWarning, match="knnbs"):
        mixed = estimator.knnbs(
            size=1,
            subclusters_maximizing_distances=["A"],
            subclusters_minimizing_distances=["B"],
        )

    assert closest.name == "closest"
    assert furthest.name == "furthest"
    assert closest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert furthest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert mixed.dropna().to_dict() == {"c1": "A", "c3": "B"}

    with pytest.warns(FutureWarning, match="find_closest"):
        all_from_a = estimator.find_closest_cells_to_self_barycenter(
            size=5,
            clusters=["A"],
        )
    assert all_from_a.dropna().to_dict() == {"c1": "A", "c2": "A"}

    with pytest.warns(FutureWarning, match="knnbs"):
        with pytest.raises(RuntimeError, match="not disjoint"):
            estimator.knnbs(
                subclusters_maximizing_distances=["A"],
                subclusters_minimizing_distances=["A"],
            )


def test_knnsc_new_api_names_are_available():
    estimator = bt.sct.tl.KNNSC()
    _set_manual_shortest_path_lengths(estimator)

    closest = estimator.select_central_cells(subcluster_size=1, key="closest")
    closest_without_name = estimator.select_central_cells(subcluster_size=1, key=None)
    furthest = estimator.select_peripheral_cells(subcluster_size=1, key="furthest")
    furthest_all_from_a = estimator.select_peripheral_cells(
        subcluster_size=5,
        clusters=["A"],
    )
    mixed = estimator.predict(
        subcluster_size=1,
        peripheral_clusters=["A"],
        central_clusters=["B"],
    )
    peripheral_only = estimator.predict(
        subcluster_size=1,
        peripheral_clusters=["A"],
    )
    central_only = estimator.predict(
        subcluster_size=1,
        central_clusters=["B"],
    )

    assert closest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert closest_without_name.name is None
    assert closest_without_name.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert furthest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert furthest_all_from_a.dropna().to_dict() == {"c1": "A", "c2": "A"}
    assert mixed.name == "knnsc"
    assert mixed.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert peripheral_only.dropna().to_dict() == {"c1": "A"}
    assert central_only.dropna().to_dict() == {"c3": "B"}
    assert list(central_only.cat.categories) == ["B"]


def test_knnsc_reports_empty_and_missing_candidate_clusters():
    estimator = bt.sct.tl.KNNSC()
    estimator._obs = pd.Series(
        pd.Categorical(
            ["A", "A", "B", "B"],
            categories=["A", "B"],
        ),
        index=["c1", "c2", "c3", "c4"],
        name="cluster",
    )
    estimator._cluster_counts = pd.Series(
        [2, 2, 0],
        index=pd.CategoricalIndex(["A", "B", "C"]),
    )
    estimator._min_cluster_size = 1

    with pytest.raises(ValueError) as exc_info:
        estimator._validate_candidate_clusters(
            ["C", "missing"],
            argument="peripheral_clusters",
        )

    message = str(exc_info.value)
    assert "empty: C" in message
    assert "not found: missing" in message


def test_shared_neighbors_pruning_and_inplace_modes(
    mini_adata,
    expected_mini_snn_connectivities,
):
    result = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=None,
        normalize_connectivities=False,
        distances_key="raw_snn_distances",
        connectivities_key="raw_snn_connectivities",
        copy=True,
    )

    assert "raw_snn_distances" in result.obsp
    assert "raw_snn_connectivities" in result.obsp
    assert result.uns["shared_neighbors"]["params"]["metric"] == "euclidean"
    assert np.allclose(
        result.obsp["raw_snn_connectivities"].toarray(),
        expected_mini_snn_connectivities,
    )

    dense_adata = mini_adata.copy()
    dense_adata.obsp["distances"] = dense_adata.obsp["distances"].toarray()
    dense_result = bt.sct.tl.shared_neighbors(
        dense_adata,
        prune_snn=1,
        snn_key="dense_snn",
        copy=True,
    )

    assert np.allclose(
        dense_result.obsp["dense_snn_connectivities"].toarray(),
        expected_mini_snn_connectivities / 3,
    )

    returned = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=0,
        snn_key="snn_inplace",
        copy=False,
    )

    assert returned is None
    assert "snn_inplace" in mini_adata.uns

    with pytest.raises(ValueError, match="expected prune_snn < n_neighbors"):
        bt.sct.tl.shared_neighbors(mini_adata, prune_snn=2)
