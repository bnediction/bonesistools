#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.figure import Figure

import bonesistools as bt
from bonesistools.sctools.plotting import _scatterplot

bt.sct.pl.set_default_params(tex=False)

ADATA = bt.sct.datasets.nestorowa()


def test_embedding_plot_with_test_representation():
    adata = ADATA.copy()

    adata.obsm["X_test"] = adata.X[:, :2].copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    fig, ax = bt.sct.pl.embedding(
        adata,
        obs="cluster",
        use_rep="X_test",
        xlabel="X1",
        ylabel="X2",
        figwidth=6,
        s=2,
        alpha=1,
        show_legend=True,
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


def test_embedding_alias_with_test_representation(mini_adata):
    fig, ax = bt.sct.pl.embedding(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    plt.close(fig)


def test_embedding_plot_continuous_3d_with_title_and_labels(mini_adata):
    mini_adata.obs["score_with_nan"] = mini_adata.obs["score"].astype(float)
    mini_adata.obs.loc["c3", "score_with_nan"] = np.nan

    fig, ax = bt.sct.pl.embedding(
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


def test_embedding_plot_discrete_handles_nan_mapping_and_labels(mini_adata):
    mini_adata.obs["label_with_nan"] = mini_adata.obs["cluster"].astype(object)
    mini_adata.obs.loc["c4", "label_with_nan"] = np.nan
    mini_adata.obs["label_with_nan"] = mini_adata.obs["label_with_nan"].astype(
        "category"
    )

    fig, ax = bt.sct.pl.embedding(
        mini_adata,
        obs="label_with_nan",
        use_rep="X_pca",
        colors={"A": [1.0, 0.0, 0.0], "B": [0.0, 0.0, 1.0]},
        show_legend=True,
        show_labels=True,
        automatic_resize=True,
        text={"fontsize": 8},
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_legend() is not None
    assert {text.get_text() for text in ax.texts} == {"A", "B"}
    plt.close(fig)


def test_embedding_discrete_uses_colors_from_uns(mini_adata):
    mini_adata.uns["cluster_color"] = {"A": "red", "B": "blue"}

    fig, ax = bt.sct.pl.embedding(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
    )

    facecolors = [collection.get_facecolor()[0] for collection in ax.collections]

    assert len(facecolors) == 2
    np.testing.assert_allclose(facecolors[0], to_rgba("red", alpha=0.3))
    np.testing.assert_allclose(facecolors[1], to_rgba("blue", alpha=0.3))
    plt.close(fig)


def test_embedding_plot_writes_outfile_and_validates_inputs(mini_adata, tmp_path):
    outfile = tmp_path / "embedding.png"

    result = bt.sct.pl.embedding(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()

    with pytest.raises(ValueError, match="invalid argument value for 'n_components'"):
        bt.sct.pl.embedding(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            n_components=4,
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'title'"):
        bt.sct.pl.embedding(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            title=object(),
        )

    with pytest.raises(KeyError):
        bt.sct.pl.embedding(mini_adata, obs="missing", use_rep="X_pca")

    with pytest.raises(KeyError):
        bt.sct.pl.embedding(mini_adata, obs="cluster", use_rep="missing")


def test_embedding_plot_discrete_3d_reuses_axes_and_customizes_ticks(mini_adata):
    mini_adata.obs["cluster_with_nan"] = mini_adata.obs["cluster"].astype(object)
    mini_adata.obs.loc["c2", "cluster_with_nan"] = np.nan
    mini_adata.obs["cluster_with_nan"] = mini_adata.obs["cluster_with_nan"].astype(
        "category"
    )

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    returned_fig, returned_ax = bt.sct.pl.embedding(
        mini_adata,
        obs="cluster_with_nan",
        use_rep="X_pca",
        n_components=3,
        ax=ax,
        colors={
            "A": np.array([1.0, 0.0, 0.0, 1.0]),
            "B": np.array([0.0, 0.0, 1.0, 1.0]),
        },
        show_legend=True,
        legend_params={"loc": "upper left"},
        tick_params={"labelsize": 6},
        title="discrete 3d",
        nan={"facecolors": "gray", "edgecolors": "none", "alpha": 0.2},
        ztick_params={"labelsize": 5},
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

    returned_fig, returned_ax = bt.sct.pl.embedding(
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


def test_embedding_plot_continuous_accepts_colormap_name(mini_adata):
    mini_adata.obs["score"] = mini_adata.obs["score"].astype(float)

    fig, ax = bt.sct.pl.embedding(
        mini_adata,
        obs="score",
        use_rep="X_pca",
        colors="gnuplot",
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert len(fig.axes) == 2
    plt.close(fig)


def test_embedding_plot_discrete_3d_default_colors_legend_and_labels(
    mini_adata,
    monkeypatch,
):
    categories = [f"cat{i}" for i in range(4)]
    monkeypatch.setattr(
        _scatterplot,
        "QUALITATIVE_COLORS",
        _scatterplot.QUALITATIVE_COLORS[:1],
    )
    mini_adata.obs["many_categories"] = pd.Categorical(
        categories,
        categories=categories,
    )

    fig, ax = bt.sct.pl.embedding(
        mini_adata,
        obs="many_categories",
        use_rep="X_pca",
        n_components=3,
        show_legend=True,
        show_labels=True,
        background_visible=False,
        lgd_params={"title": "groups"},
    )

    assert ax.name == "3d"
    assert ax.get_legend() is not None
    assert {text.get_text() for text in ax.texts} == {
        "cat0",
        "cat1",
        "cat2",
        "cat3",
    }
    plt.close(fig)


def test_embedding_plot_discrete_listed_colormap_skips_unused_category(mini_adata):
    mini_adata.obs["cluster_with_unused"] = pd.Categorical(
        ["A", "A", "B", "B"],
        categories=["A", "B", "C"],
    )

    with pytest.warns(RuntimeWarning, match="Mean of empty slice"):
        fig, ax = bt.sct.pl.embedding(
            mini_adata,
            obs="cluster_with_unused",
            use_rep="X_pca",
            colors=ListedColormap(["red", "blue", "green"]),
            show_legend=True,
            show_labels=True,
            text={"verticalalignment": "bottom"},
        )

    handles, labels = ax.get_legend_handles_labels()
    assert len(handles) == 2
    assert labels == ["A", "B"]
    assert {text.get_text() for text in ax.texts} == {"A", "B"}
    plt.close(fig)


def test_embedding_plot_rejects_unsupported_object_dtype(mini_adata):
    mini_adata.obs["object_values"] = pd.Series(
        [(1,), (2,), (1,), (2,)],
        index=mini_adata.obs.index,
        dtype=object,
    )

    with pytest.raises(TypeError, match="unsupported dtype"):
        bt.sct.pl.embedding(
            mini_adata,
            obs="object_values",
            use_rep="X_pca",
        )
