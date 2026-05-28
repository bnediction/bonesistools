#!/usr/bin/env python

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.plotting import _density

ADATA = bt.sct.datasets.nestorowa()


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
    plt.close(fig)


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
    plt.close(fig)


def test_kde_plot_with_obs_layer_mapping_ax_outfile_and_errors(mini_adata, tmp_path):
    fig, ax = plt.subplots()
    called = []
    mini_adata.layers["sparse_counts"] = csr_matrix(mini_adata.layers["counts"])

    returned_fig, returned_ax = bt.sct.pl.kde_plot(
        mini_adata,
        gene="g1",
        layer="sparse_counts",
        obs="cluster",
        colors={
            "A": np.array([1.0, 0.0, 0.0, 1.0]),
            "B": np.array([0.0, 0.0, 1.0, 1.0]),
        },
        not_all=True,
        clip=True,
        title="kde",
        default_parameters=lambda: called.append(True),
        ax=ax,
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.get_title() == "kde"
    assert ax.get_legend() is not None
    assert called == [True]
    plt.close(fig)

    outfile = tmp_path / "kde.png"
    assert bt.sct.pl.kde_plot(mini_adata, gene="g1", outfile=outfile) is None
    assert outfile.exists()

    with pytest.raises(ValueError, match="invalid argument values"):
        bt.sct.pl.kde_plot(mini_adata, gene="g1", not_all=True)

    with pytest.raises(KeyError, match="key 'missing' not found in adata.obs"):
        bt.sct.pl.kde_plot(mini_adata, gene="g1", obs="missing")


def test_kde_plot_default_generated_and_listed_colormap_colors(
    mini_adata,
    monkeypatch,
):
    fig, ax = bt.sct.pl.kde_plot(
        mini_adata,
        gene="g1",
        obs="cluster",
        not_all=True,
    )

    assert ax.get_legend() is not None
    plt.close(fig)

    monkeypatch.setattr(_density, "QUALITATIVE_COLORS", _density.QUALITATIVE_COLORS[:1])
    adata = ad.AnnData(
        X=np.array([[1.0], [2.0], [2.0], [3.0], [3.0], [4.0], [4.0], [5.0]]),
        obs=pd.DataFrame(
            {"many_clusters": pd.Categorical(["A", "A", "B", "B", "C", "C", "D", "D"])},
            index=[f"c{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    fig, ax = bt.sct.pl.kde_plot(
        adata,
        gene="g1",
        obs="many_clusters",
        title={"label": "kde generated colors"},
    )

    assert ax.get_title() == "kde generated colors"
    assert ax.get_legend() is not None
    plt.close(fig)

    mini_adata.obs["cluster"] = mini_adata.obs["cluster"].cat.remove_unused_categories()
    fig, ax = bt.sct.pl.kde_plot(
        mini_adata,
        gene="g1",
        obs="cluster",
        colors=ListedColormap(["red", "blue"]),
        not_all=True,
    )

    assert ax.get_legend() is not None
    plt.close(fig)


def test_ecdf_plot_with_obs_mapping_ax_outfile(mini_adata, tmp_path):
    fig, ax = plt.subplots()
    called = []
    mini_adata.layers["sparse_counts"] = csr_matrix(mini_adata.layers["counts"])

    returned_fig, returned_ax = bt.sct.pl.ecdf_plot(
        mini_adata,
        gene="g1",
        layer="sparse_counts",
        obs="cluster",
        xlabel=None,
        ylabel="ecdf",
        default_parameters=lambda: called.append(True),
        ax=ax,
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == "ecdf"
    assert ax.get_legend() is not None
    assert called == [True]
    plt.close(fig)

    outfile = tmp_path / "ecdf.png"
    assert bt.sct.pl.ecdf_plot(mini_adata, gene="g1", outfile=outfile) is None
    assert outfile.exists()


def test_ecdf_plot_uses_mapping_colors_for_groups(mini_adata):
    fig, ax = bt.sct.pl.ecdf_plot(
        mini_adata,
        gene="g1",
        obs="cluster",
        colors={"A": "red", "B": "blue"},
        show_legend=True,
    )

    assert ax.get_legend() is not None
    assert [line.get_color() for line in ax.lines] == [_density.gray, "red", "blue"]
    plt.close(fig)
