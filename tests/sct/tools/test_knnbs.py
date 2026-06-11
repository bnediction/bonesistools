#!/usr/bin/env python

import warnings

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import pytest

import bonesistools as bt

ADATA = bt.sct.datasets.nestorowa()


def _set_manual_shortest_path_lengths(estimator):
    estimator.obs = pd.Series(
        pd.Categorical(["A", "A", "B", "B"]),
        index=["c1", "c2", "c3", "c4"],
        name="cluster",
    )
    estimator.shortest_path_lengths_df = pd.DataFrame(
        {
            "A": [0.1, 0.8, 2.0, 1.5],
            "B": [2.0, 1.5, 0.2, 0.9],
        },
        index=estimator.obs.index,
    )


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

    adata.obsm["X_test"] = adata.X[:, :5].copy()

    graph = bt.sct.tl.knn_graph(
        adata,
        n_neighbors=5,
        use_rep="X_test",
        index_or_name="name",
    )

    assert isinstance(graph, nx.DiGraph)
    assert graph.number_of_nodes() == adata.n_obs
    assert graph.number_of_edges() > 0
    assert set(graph.nodes) == set(adata.obs_names)


def test_knnbs_fits_and_returns_subclusters_with_deprecated_api():

    adata = ADATA[:200, :].copy()

    adata.obsm["X_test"] = adata.X[:, :5].copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    knnbs = bt.sct.tl.KNNSC(
        n_neighbors=5,
        use_rep="X_test",
        n_components=5,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)

        with pytest.warns(DeprecationWarning, match="obs"):
            knnbs.fit(adata, obs="cluster", n_jobs=1)

    with pytest.warns(DeprecationWarning, match="shortest_path_lengths"):
        knnbs.shortest_path_lengths(n_jobs=1)

    with pytest.warns(DeprecationWarning, match="knnbs"):
        subclusters = knnbs.knnbs(size=5, key="knnbs")

    assert hasattr(knnbs, "knn_graph")
    assert hasattr(knnbs, "shortest_path_lengths_df")

    assert isinstance(knnbs.knn_graph, nx.Graph)
    assert isinstance(knnbs.shortest_path_lengths_df, pd.DataFrame)
    assert isinstance(subclusters, pd.Series)

    assert subclusters.name == "knnbs"
    assert subclusters.index.equals(adata.obs.index)
    assert subclusters.notna().sum() > 0
    assert set(subclusters.dropna().cat.categories) == set(
        adata.obs["cluster"].cat.categories
    )


def test_knnbs_fit_reconnects_disconnected_neighbor_graph():
    cloud_a = _grid_cloud(center=(-5.0, -5.0))
    cloud_b = _grid_cloud(center=(5.0, 5.0))
    adata, _, _ = _make_two_cluster_adata(
        cloud_a=cloud_a,
        cloud_b=cloud_b,
        names_a=[f"a{i}" for i in range(len(cloud_a))],
        names_b=[f"b{i}" for i in range(len(cloud_b))],
    )

    estimator = bt.sct.tl.KNNSC(n_neighbors=20, use_rep="X_toy")

    with pytest.warns(UserWarning, match="not weakly connected"):
        estimator.fit(adata, cluster_key="cluster", n_jobs=1)

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


def test_knnbs_select_central_cells_recovers_nearest_single_cluster_cells():
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
    estimator = bt.sct.tl.KNNSC(n_neighbors=10, use_rep="X_toy")

    estimator.fit(adata, cluster_key="cluster", n_jobs=1)

    estimator.compute_shortest_path_lengths(n_jobs=1)
    central = estimator.select_central_cells(subcluster_size=5)

    assert set(central[central == "A"].index) == {
        f"central_{i}" for i in range(len(central_points))
    }


def test_knnbs_select_peripheral_cells_recovers_remote_cells():
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
    estimator = bt.sct.tl.KNNSC(n_neighbors=16, use_rep="X_toy")

    with pytest.warns(UserWarning, match="not weakly connected"):
        estimator.fit(adata, cluster_key="cluster", n_jobs=1)

    estimator.compute_shortest_path_lengths(n_jobs=1)
    peripheral = estimator.select_peripheral_cells(subcluster_size=5)

    assert set(peripheral[peripheral == "A"].index) == {
        f"a_remote_{i}" for i in range(5)
    }
    assert set(peripheral[peripheral == "B"].index) == {
        f"b_remote_{i}" for i in range(5)
    }


def test_knnbs_validates_init_and_fit_arguments_and_repr(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'n_neighbors'"):
        bt.sct.tl.KNNSC(n_neighbors=0)

    with pytest.raises(ValueError, match="expected integer"):
        bt.sct.tl.KNNSC(n_neighbors=1.5)

    with pytest.raises(TypeError, match="unsupported argument type for 'n_neighbors'"):
        bt.sct.tl.KNNSC(n_neighbors="5")

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.tl.KNNSC(n_neighbors=5, n_components=0)

    with pytest.raises(ValueError, match="expected integer"):
        bt.sct.tl.KNNSC(n_neighbors=5, n_components=1.2)

    with pytest.raises(TypeError, match="unsupported argument type for 'n_components'"):
        bt.sct.tl.KNNSC(n_neighbors=5, n_components="2")

    with pytest.raises(TypeError, match="unsupported argument type for 'use_rep'"):
        bt.sct.tl.KNNSC(n_neighbors=5, use_rep=None)

    with pytest.raises(ValueError, match="invalid argument value for 'metric'"):
        bt.sct.tl.KNNSC(n_neighbors=5, metric="not-a-metric")

    estimator = bt.sct.tl.KNNSC(
        n_neighbors=5.0,
        n_components=2.0,
        metric="cosine",
    )

    assert "KNNSC(n_neighbors=5" in repr(estimator)

    with pytest.warns(DeprecationWarning, match="metric_kwds"):
        assert estimator.metric_kwds == {}
    with pytest.warns(DeprecationWarning, match="metric_kwds"):
        estimator.metric_kwds = {"p": 1}
    assert estimator.metric_kwargs == {"p": 1}

    with pytest.warns(DeprecationWarning, match="obs"):
        with pytest.raises(TypeError, match="received both 'cluster_key'"):
            estimator.fit(mini_adata, cluster_key="cluster", obs="cluster")

    with pytest.raises(TypeError, match="missing required argument: 'cluster_key'"):
        estimator.fit(mini_adata)


def test_knnbs_deprecated_selection_modes_and_overlap_error():
    knnbs = bt.sct.tl.KNNSC(n_neighbors=1)
    _set_manual_shortest_path_lengths(knnbs)

    with pytest.warns(DeprecationWarning, match="find_closest"):
        closest = knnbs.find_closest_cells_to_self_barycenter(size=1, key="closest")
    with pytest.warns(DeprecationWarning, match="find_furthest"):
        furthest = knnbs.find_furthest_cells_to_other_barycenters(
            size=1,
            key="furthest",
        )
    with pytest.warns(DeprecationWarning, match="knnbs"):
        mixed = knnbs.knnbs(
            size=1,
            subclusters_maximizing_distances=["A"],
            subclusters_minimizing_distances=["B"],
        )

    assert closest.name == "closest"
    assert furthest.name == "furthest"
    assert closest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert furthest.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert mixed.dropna().to_dict() == {"c1": "A", "c3": "B"}

    with pytest.warns(DeprecationWarning, match="find_closest"):
        all_from_a = knnbs.find_closest_cells_to_self_barycenter(
            size=5,
            clusters=["A"],
        )
    assert all_from_a.dropna().to_dict() == {"c1": "A", "c2": "A"}

    with pytest.warns(DeprecationWarning, match="knnbs"):
        with pytest.raises(RuntimeError, match="not disjoint"):
            knnbs.knnbs(
                subclusters_maximizing_distances=["A"],
                subclusters_minimizing_distances=["A"],
            )


def test_knnbs_new_api_names_are_available():
    estimator = bt.sct.tl.KNNSC(
        n_neighbors=1,
    )
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
    assert mixed.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert peripheral_only.dropna().to_dict() == {"c1": "A"}
    assert central_only.dropna().to_dict() == {"c3": "B"}


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
