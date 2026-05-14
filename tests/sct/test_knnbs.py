import warnings

import bonesistools as bt
import networkx as nx
import pandas as pd

ADATA = bt.sct.datasets.nestorowa()


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


def test_knnbs_fits_and_returns_subclusters():

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

        knnbs.fit(adata, obs="cluster", n_jobs=1)

    knnbs.shortest_path_lengths(n_jobs=1)

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
