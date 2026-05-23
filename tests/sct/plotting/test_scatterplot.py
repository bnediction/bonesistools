#!/usr/bin/env python

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import bonesistools as bt
from bonesistools.sctools.plotting import _colors

bt.sct.pl.set_default_params()

ADATA = bt.sct.datasets.nestorowa()


def test_embedding_plot_with_test_representation():
    adata = ADATA.copy()

    adata.obsm["X_test"] = adata.X[:, :2].copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    fig, ax = bt.sct.pl.embedding_plot(
        adata,
        obs="cluster",
        use_rep="X_test",
        xlabel="X1",
        ylabel="X2",
        figwidth=6,
        s=2,
        alpha=1,
        add_legend=True,
        lgd_params={
            "title": "clusters",
            "ncol": 1,
            "markerscale": 5,
            "frameon": True,
            "edgecolor": bt.sct.pl.get_color("black"),
            "shadow": False,
        },
        n_components=2,
        background_visible=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    plt.close(fig)


def test_embedding_plot_continuous_3d_with_title_and_labels(mini_adata):
    mini_adata.obs["score_with_nan"] = mini_adata.obs["score"].astype(float)
    mini_adata.obs.loc["c3", "score_with_nan"] = np.nan

    fig, ax = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="score_with_nan",
        use_rep="X_pca",
        n_components=3,
        title={"label": "continuous"},
        xlabel="x",
        ylabel="y",
        zlabel="z",
        background_visible=False,
        s=8,
    )

    assert isinstance(fig, Figure)
    assert ax.name == "3d"
    assert ax.get_title() == "continuous"
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"
    assert ax.get_zlabel() == "z"
    plt.close(fig)


def test_embedding_plot_discrete_handles_nan_mapping_labels_and_graph(mini_adata):
    mini_adata.obs["label_with_nan"] = mini_adata.obs["cluster"].astype(object)
    mini_adata.obs.loc["c4", "label_with_nan"] = np.nan
    mini_adata.obs["label_with_nan"] = mini_adata.obs["label_with_nan"].astype(
        "category"
    )

    epg = nx.Graph()
    epg.add_node("n1", pos=np.array([0.0, 0.0]))
    epg.add_node("n2", pos=np.array([2.0, 2.0]))
    flat_tree = nx.Graph()
    flat_tree.add_node("n1", pos=np.array([0.0, 0.0]), label="N1")
    flat_tree.add_node("n2", pos=np.array([2.0, 2.0]), label="N2")
    flat_tree.add_edge("n1", "n2", nodes=["n1", "n2"])
    mini_adata.uns["epg"] = epg
    mini_adata.uns["flat_tree"] = flat_tree

    fig, ax = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="label_with_nan",
        use_rep="X_pca",
        colors={"A": [1.0, 0.0, 0.0], "B": [0.0, 0.0, 1.0]},
        add_legend=True,
        add_labels=True,
        add_graph=True,
        add_labels_to_graph=True,
        automatic_resize=True,
        text={"fontsize": 8},
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_legend() is not None
    assert {text.get_text() for text in ax.texts}.issuperset({"A", "B", "N1", "N2"})
    assert len(ax.lines) >= 1
    plt.close(fig)


def test_embedding_plot_writes_outfile_and_validates_inputs(mini_adata, tmp_path):
    outfile = tmp_path / "embedding.png"

    result = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.pl.embedding_plot(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            n_components=4,
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'title'"):
        bt.sct.pl.embedding_plot(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            title=object(),
        )

    with pytest.raises(KeyError):
        bt.sct.pl.embedding_plot(mini_adata, obs="missing", use_rep="X_pca")

    with pytest.raises(KeyError):
        bt.sct.pl.embedding_plot(mini_adata, obs="cluster", use_rep="missing")


def test_embedding_plot_discrete_3d_reuses_axes_and_customizes_ticks(mini_adata):
    mini_adata.obs["cluster_with_nan"] = mini_adata.obs["cluster"].astype(object)
    mini_adata.obs.loc["c2", "cluster_with_nan"] = np.nan
    mini_adata.obs["cluster_with_nan"] = mini_adata.obs["cluster_with_nan"].astype(
        "category"
    )

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    returned_fig, returned_ax = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="cluster_with_nan",
        use_rep="X_pca",
        n_components=3,
        ax=ax,
        colors={
            "A": np.array([1.0, 0.0, 0.0, 1.0]),
            "B": np.array([0.0, 0.0, 1.0, 1.0]),
        },
        add_legend=True,
        legend_params={"loc": "upper left"},
        tick_params={"labelsize": 6},
        title="discrete 3d",
        nan={"facecolors": "gray", "edgecolors": "none", "alpha": 0.2},
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.name == "3d"
    assert ax.get_title() == "discrete 3d"
    assert ax.get_legend() is not None
    plt.close(fig)


def test_embedding_plot_continuous_2d_colormap_axes_and_defaults(mini_adata):
    called = []
    mini_adata.obs["score_with_nan"] = mini_adata.obs["score"].astype(float)
    mini_adata.obs.loc["c1", "score_with_nan"] = np.nan

    fig, ax = plt.subplots()

    returned_fig, returned_ax = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="score_with_nan",
        use_rep="X_pca",
        colors=plt.get_cmap("viridis"),
        n_components=2,
        ax=ax,
        default_parameters=lambda: called.append(True),
        xlabel=None,
        ylabel=None,
        xtick_params={"labelsize": 5},
        ytick_params={"labelsize": 5},
        nan={"facecolor": "lightgray", "edgecolor": "black", "alpha": 0.4},
        colorbar_scale=0.5,
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert called == [True]
    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == ""
    assert len(fig.axes) == 2
    plt.close(fig)


def test_embedding_plot_discrete_3d_default_colors_legend_labels_and_graph(
    mini_adata,
    monkeypatch,
):
    categories = [f"cat{i}" for i in range(4)]
    monkeypatch.setattr(_colors, "QUALITATIVE_COLORS", _colors.QUALITATIVE_COLORS[:1])
    mini_adata.obs["many_categories"] = pd.Categorical(
        categories,
        categories=categories,
    )

    epg = nx.Graph()
    epg.add_node("n1", pos=np.array([0.0, 0.0, 0.0]))
    epg.add_node("n2", pos=np.array([2.0, 2.0, 1.0]))
    flat_tree = nx.Graph()
    flat_tree.add_node("n1", pos=np.array([0.0, 0.0, 0.0]), label="N1")
    flat_tree.add_node("n2", pos=np.array([2.0, 2.0, 1.0]), label="N2")
    flat_tree.add_edge("n1", "n2", nodes=["n1", "n2"])
    mini_adata.uns["epg"] = epg
    mini_adata.uns["flat_tree"] = flat_tree

    fig, ax = bt.sct.pl.embedding_plot(
        mini_adata,
        obs="many_categories",
        use_rep="X_pca",
        n_components=3,
        add_legend=True,
        add_labels=True,
        add_graph=True,
        add_labels_to_graph=True,
        background_visible=False,
    )

    assert ax.name == "3d"
    assert ax.get_legend() is not None
    assert {text.get_text() for text in ax.texts}.issuperset(
        {"cat0", "cat1", "cat2", "cat3", "N1", "N2"}
    )
    assert len(ax.lines) >= 1
    plt.close(fig)
