#!/usr/bin/env python

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pytest
from matplotlib.axes import Axes
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.plotting import _graph as graph_plotting


def _add_paga_edges(adata):
    adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.4, 0.0]])


def _add_trajectory_graph(adata, graph_key="epg"):
    graph = nx.Graph()
    graph.add_node("n1", pos=np.array([0.0, 0.0, 0.0]), label="N1")
    graph.add_node("n2", pos=np.array([2.0, 2.0, 1.0]), label="N2")
    graph.add_edge("n1", "n2")
    adata.uns[graph_key] = graph


@pytest.mark.parametrize("representation", ["X_paga_2d", "X_pca"])
def test_paga_reuses_axes_with_2d_and_3d_representations(mini_adata, representation):
    _add_paga_edges(mini_adata)
    mini_adata.obsm["X_paga_2d"] = mini_adata.obsm["X_pca"][:, :2].copy()

    fig, ax = plt.subplots()
    returned_ax = bt.sct.pl.paga(
        mini_adata,
        obs="cluster",
        representation=representation,
        edges="paga_edges",
        threshold=0.1,
        ax=ax,
        with_labels=True,
        node_color={"A": "red", "B": "blue"},
    )

    assert returned_ax is ax
    assert isinstance(returned_ax, Axes)
    plt.close(fig)


def test_paga_writes_outfile(mini_adata, tmp_path):
    _add_paga_edges(mini_adata)

    outfile = tmp_path / "paga.png"
    result = bt.sct.pl.paga(
        mini_adata,
        obs="cluster",
        representation="X_pca",
        edges="paga_edges",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()


def test_paga_graph_matches_scanpy_transition_orientation(
    mini_adata,
    expected_mini_cluster_barycenters,
):
    mini_adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])

    graph = graph_plotting._paga_graph(
        mini_adata,
        obs="cluster",
        representation="X_pca",
        edges="paga_edges",
        threshold=0.1,
    )

    assert graph.is_directed()
    assert set(graph.nodes) == {"A", "B"}
    assert list(graph.edges) == [("B", "A")]
    assert np.allclose(graph.nodes["A"]["pos"], expected_mini_cluster_barycenters["A"])
    assert np.allclose(graph.nodes["B"]["pos"], expected_mini_cluster_barycenters["B"])


def test_paga_graph_reads_nested_connectivities_as_undirected(mini_adata):
    bt.sct.tl.paga(mini_adata, groupby="cluster")

    graph = graph_plotting._paga_graph(
        mini_adata,
        obs="cluster",
        representation="X_pca",
        edges="connectivities",
        threshold=0.1,
    )

    assert not graph.is_directed()
    assert set(graph.edges) == {("A", "B")}


def test_paga_uses_nested_connectivities_after_tools_paga(mini_adata):
    bt.sct.tl.paga(mini_adata, groupby="cluster")

    fig, ax = plt.subplots()
    returned_ax = bt.sct.pl.paga(
        mini_adata,
        obs="cluster",
        representation="X_pca",
        threshold=0.1,
        ax=ax,
    )

    assert returned_ax is ax
    plt.close(fig)


def test_trajectory_uses_custom_graph_key_and_draws_labels(mini_adata):
    _add_trajectory_graph(mini_adata, graph_key="custom_key")

    fig, ax = bt.sct.pl.trajectory(
        mini_adata,
        obs="cluster",
        representation="X_pca",
        graph_key="custom_key",
        show_labels=True,
        graph={"linewidth": 5},
    )

    assert len(ax.lines) == 1
    assert ax.lines[0].get_linewidth() == 5
    assert {text.get_text() for text in ax.texts} == {"N1", "N2"}
    plt.close(fig)


def test_trajectory_deprecates_obsm(mini_adata):
    _add_trajectory_graph(mini_adata)

    with pytest.warns(FutureWarning, match="`obsm` is deprecated"):
        fig, ax = bt.sct.pl.trajectory(
            mini_adata,
            obs="cluster",
            obsm="X_pca",
        )

    assert isinstance(ax, Axes)
    plt.close(fig)


def test_trajectory_missing_graph_key_raises_clear_key_error(mini_adata):
    with pytest.raises(KeyError, match="missing_graph"):
        bt.sct.pl.trajectory(
            mini_adata,
            obs="cluster",
            representation="X_pca",
            graph_key="missing_graph",
        )


def test_paga_deprecates_obsm(mini_adata):
    _add_paga_edges(mini_adata)

    with pytest.warns(FutureWarning, match="`obsm` is deprecated"):
        ax = bt.sct.pl.paga(
            mini_adata,
            obs="cluster",
            obsm="X_pca",
            edges="paga_edges",
        )

    assert isinstance(ax, Axes)


def test_embedding_deprecated_graph_arguments_warn_and_route_overlay(mini_adata):
    _add_trajectory_graph(mini_adata)

    with pytest.warns(FutureWarning, match="add_graph"):
        fig, ax = bt.sct.pl.embedding(
            mini_adata,
            obs="cluster",
            representation="X_pca",
            add_graph=True,
        )

    assert len(ax.lines) == 1
    plt.close(fig)

    with pytest.warns(FutureWarning, match="add_labels_to_graph"):
        fig, ax = bt.sct.pl.embedding(
            mini_adata,
            obs="cluster",
            representation="X_pca",
            add_labels_to_graph=True,
        )

    assert {text.get_text() for text in ax.texts} == {"N1", "N2"}
    plt.close(fig)
