#!/usr/bin/env python

import networkx as nx
import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import bonesistools as bt

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
