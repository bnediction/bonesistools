#!/usr/bin/env python

import pytest

import bonesistools as bt

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import numpy as np

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


def test_kde_plot_returns_matplotlib_objects():

    adata = ADATA.copy()

    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    gene = adata.var_names[0]

    fig, ax = bt.sct.pl.kde_plot(
        adata,
        gene=gene,
        xlabel="count",
        ylabel="density",
        show_legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_ecdf_plot_returns_matplotlib_objects():

    adata = ADATA.copy()

    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    gene = adata.var_names[0]

    fig, ax = bt.sct.pl.ecdf_plot(
        adata,
        gene=gene,
        xlabel="count",
        ylabel="ECDF",
        show_legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_boxplot_without_hue():

    adata = ADATA.copy()

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    fig, ax, bps = bt.sct.pl.boxplot(
        adata,
        obs="n_counts",
        groupby="label",
        sort="ascending",
        showpoints=True,
        showlegend=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert isinstance(bps, dict)

    assert "boxes" in bps
    assert "medians" in bps
    assert len(bps["boxes"]) == adata.obs["label"].nunique()

    plt.close(fig)


def test_boxplot_with_hue():

    adata = ADATA.copy()

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    adata.obs["condition"] = "condition2"
    adata.obs.iloc[:1000, adata.obs.columns.get_loc("condition")] = "condition1"
    adata.obs["condition"] = adata.obs["condition"].astype("category")

    fig, _, bps = bt.sct.pl.boxplot(
        adata,
        obs="n_counts",
        groupby="label",
        hue="condition",
        showpoints=False,
        showlegend=True,
    )

    assert isinstance(bps, dict)
    assert set(bps.keys()) == {"condition1", "condition2"}

    for condition in adata.obs["condition"].cat.categories:
        assert "boxes" in bps[condition]
        assert len(bps[condition]["boxes"]) == adata.obs["label"].nunique()

    plt.close(fig)


def test_boxplot_invalid_sort():

    adata = ADATA.copy()

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    with pytest.raises(ValueError):
        bt.sct.pl.boxplot(
            adata,
            obs="n_counts",
            groupby="label",
            sort="invalid",
        )


def test_boxplot_outfile(tmp_path):

    mpl.rcParams["text.usetex"] = False

    adata = bt.sct.datasets.nestorowa()

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    outfile = tmp_path / "boxplot.png"

    result = bt.sct.pl.boxplot(
        adata,
        obs="n_counts",
        groupby="label",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()
