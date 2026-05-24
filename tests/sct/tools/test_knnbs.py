#!/usr/bin/env python

import warnings

import bonesistools as bt
import networkx as nx
import numpy as np
import pandas as pd
import pytest

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


def test_kneighbors_graph_returns_graph_with_named_nodes():

    adata = ADATA.copy()

    adata.obsm["X_test"] = adata.X[:, :5].copy()

    graph = bt.sct.tl.kneighbors_graph(
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

    knnbs = bt.sct.tl.Knnbs(
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

    assert hasattr(knnbs, "kneighbors_graph")
    assert hasattr(knnbs, "shortest_path_lengths_df")

    assert isinstance(knnbs.kneighbors_graph, nx.Graph)
    assert isinstance(knnbs.shortest_path_lengths_df, pd.DataFrame)
    assert isinstance(subclusters, pd.Series)

    assert subclusters.name == "knnbs"
    assert subclusters.index.equals(adata.obs.index)
    assert subclusters.notna().sum() > 0
    assert set(subclusters.dropna().cat.categories) == set(
        adata.obs["cluster"].cat.categories
    )


def test_knnbs_validates_init_arguments_and_repr():
    with pytest.raises(ValueError, match="invalid argument value for 'n_neighbors'"):
        bt.sct.tl.Knnbs(n_neighbors=0)

    with pytest.raises(ValueError, match="expected integer"):
        bt.sct.tl.Knnbs(n_neighbors=1.5)

    with pytest.raises(TypeError, match="unsupported argument type for 'n_neighbors'"):
        bt.sct.tl.Knnbs(n_neighbors="5")

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.tl.Knnbs(n_neighbors=5, n_components=0)

    with pytest.raises(ValueError, match="expected integer"):
        bt.sct.tl.Knnbs(n_neighbors=5, n_components=1.2)

    with pytest.raises(TypeError, match="unsupported argument type for 'n_components'"):
        bt.sct.tl.Knnbs(n_neighbors=5, n_components="2")

    with pytest.raises(TypeError, match="unsupported argument type for 'use_rep'"):
        bt.sct.tl.Knnbs(n_neighbors=5, use_rep=None)

    with pytest.raises(ValueError, match="invalid argument value for 'metric'"):
        bt.sct.tl.Knnbs(n_neighbors=5, metric="not-a-metric")

    estimator = bt.sct.tl.Knnbs(
        n_neighbors=5.0,
        n_components=2.0,
        metric="cosine",
    )

    assert "Knnbs(n_neighbors=5.0" in repr(estimator)

    with pytest.warns(DeprecationWarning, match="metric_kwds"):
        assert estimator.metric_kwds == {}
    with pytest.warns(DeprecationWarning, match="metric_kwds"):
        estimator.metric_kwds = {"p": 1}
    assert estimator.metric_kwargs == {"p": 1}

def test_knnbs_distance_selection_modes_and_overlap_error_with_deprecated_api():
    knnbs = bt.sct.tl.Knnbs(n_neighbors=1)
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
    estimator = bt.sct.tl.Knnbs(
        n_neighbors=1,
    )
    _set_manual_shortest_path_lengths(estimator)

    closest = estimator.select_central_cells(subcluster_size=1, key="closest")
    closest_without_name = estimator.select_central_cells(subcluster_size=1, key=None)
    furthest = estimator.select_peripheral_cells(subcluster_size=1, key="furthest")
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
    assert mixed.dropna().to_dict() == {"c1": "A", "c3": "B"}
    assert peripheral_only.dropna().to_dict() == {"c1": "A"}
    assert central_only.dropna().to_dict() == {"c3": "B"}


def test_shared_neighbors_pruning_and_inplace_modes(mini_adata):
    result = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=None,
        metric=None,
        normalize_connectivities=False,
        distances_key="raw_snn_distances",
        connectivities_key="raw_snn_connectivities",
        copy=True,
    )

    assert "raw_snn_distances" in result.obsp
    assert "raw_snn_connectivities" in result.obsp
    assert result.uns["shared_neighbors"]["params"]["metric"] == "euclidean"
    assert np.all(result.obsp["raw_snn_connectivities"].data >= 1)

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
