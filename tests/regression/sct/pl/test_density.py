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
from tests.regression.sct.toy_data import make_nestorowa_hvg_adata

ADATA = make_nestorowa_hvg_adata()


def test_density_returns_matplotlib_objects():
    adata = ADATA.copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")
    feature = adata.var_names[0]

    fig, ax = bt.sct.pl.density(
        adata,
        feature=feature,
        xlabel={"label": "count", "fontsize": 13},
        ylabel={"label": "density", "fontsize": 14},
        legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_xlabel() == "count"
    assert ax.get_ylabel() == "density"
    assert ax.xaxis.label.get_fontsize() == 13
    assert ax.yaxis.label.get_fontsize() == 14
    plt.close(fig)


def test_cdf_returns_matplotlib_objects():
    adata = ADATA.copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")
    feature = adata.var_names[0]

    fig, ax = bt.sct.pl.cdf(
        adata,
        feature=feature,
        xlabel="count",
        ylabel="ECDF",
        legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    plt.close(fig)


def test_density_with_obs_expression_mapping_ax_outfile_and_errors(
    mini_adata,
    tmp_path,
):
    fig, ax = plt.subplots()
    mini_adata.layers["sparse_counts"] = csr_matrix(mini_adata.layers["counts"])

    returned_fig, returned_ax = bt.sct.pl.density(
        mini_adata,
        feature="g1",
        expression="sparse_counts",
        obs="cluster",
        colors={
            "A": np.array([1.0, 0.0, 0.0, 1.0]),
            "B": np.array([0.0, 0.0, 1.0, 1.0]),
        },
        show_global=False,
        clip_outliers=True,
        title="kde",
        legend={"title": "groups"},
        ax=ax,
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.get_title() == "kde"
    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "groups"
    plt.close(fig)

    outfile = tmp_path / "kde.png"
    assert bt.sct.pl.density(mini_adata, feature="g1", outfile=outfile) is None
    assert outfile.exists()

    with pytest.raises(ValueError, match="invalid argument values"):
        bt.sct.pl.density(mini_adata, feature="g1", show_global=False)

    with pytest.raises(KeyError):
        bt.sct.pl.density(mini_adata, feature="g1", obs="missing")


def test_density_default_generated_and_listed_colormap_colors(
    mini_adata,
    monkeypatch,
):
    fig, ax = bt.sct.pl.density(
        mini_adata,
        feature="g1",
        obs="cluster",
        show_global=False,
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

    fig, ax = bt.sct.pl.density(
        adata,
        feature="g1",
        obs="many_clusters",
        title={"label": "kde generated colors"},
    )

    assert ax.get_title() == "kde generated colors"
    assert ax.get_legend() is not None
    plt.close(fig)

    mini_adata.obs["cluster"] = mini_adata.obs["cluster"].cat.remove_unused_categories()
    fig, ax = bt.sct.pl.density(
        mini_adata,
        feature="g1",
        obs="cluster",
        colors=ListedColormap(["red", "blue"]),
        show_global=False,
    )

    assert ax.get_legend() is not None
    plt.close(fig)


def test_cdf_with_obs_mapping_ax_outfile(mini_adata, tmp_path):
    fig, ax = plt.subplots()
    mini_adata.layers["sparse_counts"] = csr_matrix(mini_adata.layers["counts"])

    returned_fig, returned_ax = bt.sct.pl.cdf(
        mini_adata,
        feature="g1",
        expression="sparse_counts",
        obs="cluster",
        xlabel=None,
        ylabel={"label": "ecdf", "fontsize": 13},
        ax=ax,
    )

    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == "ecdf"
    assert ax.yaxis.label.get_fontsize() == 13
    assert ax.get_legend() is not None
    plt.close(fig)

    outfile = tmp_path / "ecdf.png"
    assert bt.sct.pl.cdf(mini_adata, feature="g1", outfile=outfile) is None
    assert outfile.exists()


def test_cdf_uses_mapping_colors_for_groups(mini_adata):
    fig, ax = bt.sct.pl.cdf(
        mini_adata,
        feature="g1",
        obs="cluster",
        colors={"A": "red", "B": "blue"},
        legend=True,
    )

    assert ax.get_legend() is not None
    assert [line.get_color() for line in ax.lines] == [_density.gray, "red", "blue"]
    plt.close(fig)


def test_cdf_accepts_explicit_legend_kwargs(mini_adata):
    fig, ax = bt.sct.pl.cdf(
        mini_adata,
        feature="g1",
        obs="cluster",
        legend={"title": "cluster"},
    )

    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "cluster"
    plt.close(fig)


@pytest.mark.parametrize("function", [bt.sct.pl.density, bt.sct.pl.cdf])
def test_density_functions_deprecate_gene_and_layer(mini_adata, function):
    mini_adata.layers["sparse_counts"] = csr_matrix(mini_adata.layers["counts"])

    with pytest.warns(FutureWarning, match="`gene` is deprecated"):
        fig, ax = function(
            mini_adata,
            gene="g1",
            obs="cluster",
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(FutureWarning, match="`layer` is deprecated"):
        fig, ax = function(
            mini_adata,
            feature="g1",
            layer="sparse_counts",
            obs="cluster",
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(FutureWarning, match="`gene` is deprecated"):
        with pytest.raises(TypeError, match="feature.*gene"):
            function(mini_adata, feature="g1", gene="g2")

    with pytest.warns(FutureWarning, match="`layer` is deprecated"):
        with pytest.raises(TypeError, match="expression.*layer"):
            function(
                mini_adata,
                feature="g1",
                expression="counts",
                layer="sparse_counts",
            )


def test_density_deprecates_clip(mini_adata):
    with pytest.warns(FutureWarning, match="`clip` is deprecated"):
        fig, ax = bt.sct.pl.density(
            mini_adata,
            feature="g1",
            obs="cluster",
            clip=True,
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(FutureWarning, match="`clip` is deprecated"):
        with pytest.raises(TypeError, match="clip_outliers.*clip"):
            bt.sct.pl.density(
                mini_adata,
                feature="g1",
                clip_outliers=True,
                clip=True,
            )


def test_density_deprecates_not_all(mini_adata):
    with pytest.warns(FutureWarning, match="`not_all` is deprecated"):
        fig, ax = bt.sct.pl.density(
            mini_adata,
            feature="g1",
            obs="cluster",
            not_all=True,
        )

    legend = ax.get_legend()
    assert legend is not None
    assert "all" not in [text.get_text() for text in legend.texts]
    plt.close(fig)

    with pytest.warns(FutureWarning, match="`not_all` is deprecated"):
        with pytest.raises(TypeError, match="show_global.*not_all"):
            bt.sct.pl.density(
                mini_adata,
                feature="g1",
                obs="cluster",
                show_global=False,
                not_all=True,
            )


@pytest.mark.parametrize("function", [bt.sct.pl.density, bt.sct.pl.cdf])
def test_density_functions_deprecate_show_legend_without_effect(
    mini_adata,
    function,
):
    with pytest.warns(DeprecationWarning, match="'show_legend' is deprecated"):
        fig, ax = function(
            mini_adata,
            feature="g1",
            obs="cluster",
            show_legend=False,
        )

    assert ax.get_legend() is not None
    plt.close(fig)
