#!/usr/bin/env python

from types import MappingProxyType
from typing import Any, Dict, Sequence, cast

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import bonesistools as bt
from bonesistools.sctools.plotting import _distribution
from tests.regression.sct.toy_data import make_nestorowa_hvg_adata

ADATA = make_nestorowa_hvg_adata()


def _n_counts(adata: ad.AnnData) -> np.ndarray:
    return np.asarray(cast(Any, adata.X).sum(axis=1)).flatten()


def _median_lines(bp: _distribution.BoxPlots) -> Sequence[Line2D]:
    return cast(Sequence[Line2D], bp["medians"])


def test_distribution_without_hue():
    adata = ADATA.copy()
    adata.obs["n_counts"] = _n_counts(adata)

    fig, ax, bps = bt.sct.pl.distribution(
        adata,
        obs="n_counts",
        groupby="label",
        sort="ascending",
        points=True,
        legend=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert isinstance(bps, dict)
    assert "boxes" in bps
    assert "medians" in bps
    assert len(bps["boxes"]) == adata.obs["label"].nunique()
    plt.close(fig)


def test_distribution_with_hue():
    adata = ADATA.copy()
    adata.obs["n_counts"] = _n_counts(adata)
    condition = np.full(adata.n_obs, "condition2", dtype=object)
    condition[:1000] = "condition1"
    adata.obs["condition"] = pd.Categorical(condition)

    fig, _, bps = bt.sct.pl.distribution(
        adata,
        obs="n_counts",
        groupby="label",
        hue="condition",
        points=False,
    )

    assert isinstance(bps, dict)
    assert set(bps.keys()) == {"condition1", "condition2"}
    for condition in adata.obs["condition"].cat.categories:
        assert "boxes" in bps[condition]
        assert len(bps[condition]["boxes"]) == adata.obs["label"].nunique()
    plt.close(fig)


def test_distribution_invalid_sort():
    adata = ADATA.copy()
    adata.obs["n_counts"] = _n_counts(adata)

    with pytest.raises(ValueError):
        bt.sct.pl.distribution(
            adata,
            obs="n_counts",
            groupby="label",
            sort=cast(Any, "invalid"),
        )


def test_distribution_sort_preserve_keeps_category_order(mini_adata):
    mini_adata.obs["cluster"] = pd.Categorical(
        mini_adata.obs["cluster"],
        categories=["B", "A"],
        ordered=True,
    )

    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        sort="preserve",
    )

    assert [label.get_text() for label in ax.get_xticklabels()] == ["B", "A"]
    plt.close(fig)


def test_distribution_outfile(tmp_path):
    mpl.rcParams["text.usetex"] = False

    adata = ADATA.copy()
    adata.obs["n_counts"] = _n_counts(adata)
    outfile = tmp_path / "distribution.png"

    result = bt.sct.pl.distribution(
        adata,
        obs="n_counts",
        groupby="label",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()


def test_distribution_without_groupby_and_validation_errors(mini_adata):
    fig, ax, bps = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        title={"label": "score"},
        points=True,
        median=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_title() == "score"
    assert len(bps["boxes"]) == 1
    assert all(median.get_linewidth() == 0 for median in _median_lines(bps))
    plt.close(fig)

    with pytest.raises(ValueError, match="invalid argument values"):
        bt.sct.pl.distribution(mini_adata, obs="score", hue="cluster")

    with pytest.raises(TypeError, match="unsupported argument type for 'title'"):
        bt.sct.pl.distribution(mini_adata, obs="score", title=cast(Any, object()))


def test_distribution_with_hue_custom_colors_and_hidden_medians(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, bps = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        colors={"ctrl": [1.0, 0.0, 0.0], "stim": [0.0, 0.0, 1.0]},
        points={"colors": {"ctrl": [1.0, 0.5, 0.5], "stim": [0.5, 0.5, 1.0]}},
        median=False,
    )

    assert isinstance(ax, Axes)
    assert set(bps) == {"ctrl", "stim"}
    assert ax.get_legend() is not None
    for bp in bps.values():
        assert all(median.get_linewidth() == 0 for median in _median_lines(bp))
    plt.close(fig)

    with pytest.raises(TypeError, match="box_colors"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            box_colors={"ctrl": "red", "stim": "blue"},
        )

    with pytest.raises(TypeError, match="point_colors"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            point_colors={"ctrl": "pink", "stim": "lightblue"},
        )


def test_distribution_accepts_explicit_legend_kwargs(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        legend={"title": "condition", "fontsize": 12},
    )

    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "condition"
    assert legend.get_texts()[0].get_fontsize() == 12
    plt.close(fig)


def test_distribution_accepts_empty_and_non_dict_legend_mapping(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        legend={},
    )

    assert ax.get_legend() is not None
    plt.close(fig)

    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        legend=MappingProxyType({"title": "condition"}),
    )

    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "condition"
    plt.close(fig)


def test_distribution_rejects_legend_fontsize(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    with pytest.raises(TypeError, match="legend=\\{'fontsize': \\.\\.\\.\\}"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            legend_fontsize=12,
        )


def test_distribution_hue_defaults_and_listed_colormaps(mini_adata, monkeypatch):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, bps = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        title="score by condition",
        points=True,
    )

    assert ax.get_title() == "score by condition"
    assert set(bps) == {"ctrl", "stim"}
    assert ax.get_legend() is not None
    plt.close(fig)

    monkeypatch.setattr(
        _distribution,
        "QUALITATIVE_COLORS",
        _distribution.QUALITATIVE_COLORS[:1],
    )
    hue_adata = ad.AnnData(
        X=np.ones((8, 1)),
        obs=pd.DataFrame(
            {
                "score": [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5],
                "group": pd.Categorical(["A", "B"] * 4),
                "hue4": pd.Categorical(
                    ["h1", "h1", "h2", "h2", "h3", "h3", "h4", "h4"]
                ),
            },
            index=[f"c{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    box_colors = ListedColormap(["red", "green", "blue", "purple"])
    point_colors = ListedColormap(["pink", "lightgreen", "lightblue", "plum"])

    fig, _, bps = bt.sct.pl.distribution(
        hue_adata,
        obs="score",
        groupby="group",
        hue="hue4",
        colors=box_colors,
        points={"colors": point_colors},
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)

    fig, _, bps = bt.sct.pl.distribution(
        hue_adata,
        obs="score",
        groupby="group",
        hue="hue4",
        points=False,
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)

    fig, _, bps = bt.sct.pl.distribution(
        hue_adata,
        obs="score",
        groupby="group",
        hue="hue4",
        points=True,
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)


def test_distribution_accepts_existing_axes(mini_adata):
    fig, ax = plt.subplots()

    result = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        ax=ax,
    )

    assert result is not None
    returned_fig, returned_ax, _ = result
    assert returned_fig is fig
    assert returned_ax is ax
    plt.close(fig)


def test_distribution_violin_without_hue(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        colors=["red", "blue"],
        showextrema=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert "bodies" in artists
    assert len(artists["bodies"]) == mini_adata.obs["cluster"].nunique()
    assert "cbars" not in artists
    body = artists["bodies"][0]
    np.testing.assert_allclose(body.get_facecolor()[0][:3], bt.sct.pl.rgba("red")[:3])
    assert body.get_facecolor()[0][3] == 1.0
    plt.close(fig)


def test_distribution_violin_without_hue_uses_default_palette(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        showextrema=False,
    )

    expected_colors = _distribution.QUALITATIVE_COLORS[: len(artists["bodies"])]
    for body, color in zip(artists["bodies"], expected_colors):
        np.testing.assert_allclose(
            body.get_facecolor()[0][:3],
            bt.sct.pl.rgba(color)[:3],
        )
        assert body.get_facecolor()[0][3] == 1.0

    assert "cmedians" in artists
    np.testing.assert_allclose(
        artists["cmedians"].get_colors(),
        [bt.sct.pl.rgba("C1")],
    )
    assert all(
        pattern is not None for _, pattern in artists["cmedians"].get_linestyle()
    )
    np.testing.assert_allclose(artists["cmedians"].get_linewidths(), [1.0])
    np.testing.assert_allclose(
        artists["cmedians"].get_segments(),
        [
            [[-0.25, 0.4], [0.25, 0.4]],
            [[0.55, 0.5], [1.05, 0.5]],
        ],
    )

    plt.close(fig)


def test_distribution_violin_accepts_alpha(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        alpha=0.5,
        showextrema=False,
    )

    for body in artists["bodies"]:
        assert body.get_facecolor()[0][3] == 0.5

    plt.close(fig)


def test_distribution_violin_mean_matches_boxplot_style(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        widths=0.5,
        mean=True,
        showextrema=False,
    )

    assert "cmeans" in artists
    np.testing.assert_allclose(
        artists["cmeans"].get_colors(),
        [bt.sct.pl.rgba("C2")],
    )
    assert all(
        pattern is not None for _, pattern in artists["cmeans"].get_linestyle()
    )
    np.testing.assert_allclose(artists["cmeans"].get_linewidths(), [1.0])
    np.testing.assert_allclose(
        artists["cmeans"].get_segments(),
        [
            [[-0.25, 0.4], [0.25, 0.4]],
            [[0.55, 0.5], [1.05, 0.5]],
        ],
    )

    plt.close(fig)


def test_distribution_box_summary_statistics_accept_artist_kwargs(mini_adata):
    fig, _, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        median={"color": "tab:red", "linewidth": 3.0},
        mean={
            "marker": "D",
            "markerfacecolor": "tab:green",
            "markeredgecolor": "tab:green",
        },
    )

    assert all(median.get_color() == "tab:red" for median in artists["medians"])
    assert all(median.get_linewidth() == 3.0 for median in artists["medians"])
    assert artists["means"]
    assert all(mean.get_marker() == "D" for mean in artists["means"])
    assert all(
        mean.get_markerfacecolor() == "tab:green" for mean in artists["means"]
    )
    plt.close(fig)


def test_distribution_violin_summary_statistics_accept_artist_kwargs(mini_adata):
    fig, _, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        median={"color": "tab:red", "linewidth": 2.0},
        mean={"color": "tab:green", "linewidth": 3.0},
        showextrema=False,
    )

    assert "cmedians" in artists
    assert "cmeans" in artists
    np.testing.assert_allclose(
        artists["cmedians"].get_colors(),
        [bt.sct.pl.rgba("tab:red")],
    )
    np.testing.assert_allclose(artists["cmedians"].get_linewidths(), [2.0])
    np.testing.assert_allclose(
        artists["cmeans"].get_colors(),
        [bt.sct.pl.rgba("tab:green")],
    )
    np.testing.assert_allclose(artists["cmeans"].get_linewidths(), [3.0])
    plt.close(fig)


def test_distribution_single_violin_sets_width_relative_xlim(mini_adata):
    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        kind="violin",
        widths=0.5,
        showextrema=False,
    )

    np.testing.assert_allclose(ax.get_xlim(), (-0.325, 0.325))

    plt.close(fig)


def test_distribution_points_use_alpha(mini_adata):
    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        points=True,
        alpha=0.4,
    )

    points = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert len(points) == 1
    assert points[0].get_alpha() == 0.4

    plt.close(fig)


def test_distribution_points_explicit_alpha_takes_precedence(mini_adata):
    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        points={"alpha": 0.2},
        alpha=0.8,
    )

    points = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert len(points) == 1
    assert points[0].get_alpha() == 0.2

    plt.close(fig)


def test_distribution_points_accept_non_dict_mapping(mini_adata):
    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        points=MappingProxyType({"alpha": 0.25, "s": 4}),
    )

    points = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert len(points) == 1
    assert points[0].get_alpha() == 0.25
    assert points[0].get_sizes()[0] == 4
    plt.close(fig)


def test_distribution_violin_default_clip_uses_observed_range(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        showextrema=False,
    )

    groups = mini_adata.obs["cluster"].cat.categories
    for body, group in zip(artists["bodies"], groups):
        observed = mini_adata.obs.loc[mini_adata.obs["cluster"] == group, "score"]
        vertices = body.get_paths()[0].vertices
        np.testing.assert_allclose(vertices[:, 1].min(), observed.min())
        np.testing.assert_allclose(vertices[:, 1].max(), observed.max())

    plt.close(fig)


def test_distribution_violin_clip_none_extends_kde_range(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        clip=None,
        showextrema=False,
    )

    groups = mini_adata.obs["cluster"].cat.categories
    for body, group in zip(artists["bodies"], groups):
        observed = mini_adata.obs.loc[mini_adata.obs["cluster"] == group, "score"]
        vertices = body.get_paths()[0].vertices
        assert vertices[:, 1].min() < observed.min()
        assert vertices[:, 1].max() > observed.max()

    plt.close(fig)


def test_distribution_violin_cut_zero_uses_observed_range(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        cut=0,
        showextrema=False,
    )

    groups = mini_adata.obs["cluster"].cat.categories
    for body, group in zip(artists["bodies"], groups):
        observed = mini_adata.obs.loc[mini_adata.obs["cluster"] == group, "score"]
        vertices = body.get_paths()[0].vertices
        np.testing.assert_allclose(vertices[:, 1].min(), observed.min())
        np.testing.assert_allclose(vertices[:, 1].max(), observed.max())

    plt.close(fig)


def test_distribution_violin_clip_limits_kde_range(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        clip=(0, None),
        showextrema=False,
    )

    for body in artists["bodies"]:
        vertices = body.get_paths()[0].vertices
        np.testing.assert_allclose(vertices[:, 1].min(), 0.0)

    plt.close(fig)


def test_distribution_box_accepts_backend_kwargs(mini_adata):
    fig, _, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        notch=True,
    )

    assert "boxes" in artists
    plt.close(fig)


def test_distribution_violin_with_points(mini_adata):
    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        points=True,
    )

    assert "bodies" in artists
    assert any(isinstance(collection, PathCollection) for collection in ax.collections)
    plt.close(fig)


def test_distribution_violin_with_hue_legend_and_colors(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, artists = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        kind="violin",
        colors={"ctrl": "red", "stim": "blue"},
        legend={"title": "condition"},
    )

    assert set(artists) == {"ctrl", "stim"}
    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "condition"
    ctrl_body = cast(_distribution.ViolinPlots, artists["ctrl"])["bodies"][0]
    stim_body = cast(_distribution.ViolinPlots, artists["stim"])["bodies"][0]
    np.testing.assert_allclose(
        ctrl_body.get_facecolor()[0][:3],
        bt.sct.pl.rgba("red")[:3],
    )
    np.testing.assert_allclose(
        stim_body.get_facecolor()[0][:3],
        bt.sct.pl.rgba("blue")[:3],
    )
    plt.close(fig)


def test_distribution_violin_outfile(tmp_path, mini_adata):
    mpl.rcParams["text.usetex"] = False
    outfile = tmp_path / "violin.png"

    result = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        kind="violin",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()


def test_distribution_explicit_figure_arguments(mini_adata):
    fig, ax, _ = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        figwidth=4.0,
        figheight=3.0,
        xlabel="cluster",
        ylabel="score",
    )

    assert fig.get_figwidth() == 4.0
    assert fig.get_figheight() == 3.0
    assert ax.get_xlabel() == "cluster"
    assert ax.get_ylabel() == "score"
    plt.close(fig)


def test_distribution_validates_kind_and_backend_arguments(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'kind'"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            kind=cast(Any, "density"),
        )

    with pytest.raises(TypeError, match="boxplot"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            boxplot={"notch": True},
        )

    with pytest.raises(TypeError, match="violin"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            kind="violin",
            violin={"showextrema": False},
        )

    with pytest.raises(ValueError, match="unsupported keyword argument"):
        bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            kind="violin",
            notch=True,
        )


def test_distribution_position_and_point_helper_validation():
    with pytest.raises(ValueError, match="groups' and 'hues'"):
        _distribution.__get_box_positions(
            widths=0.5,
            hues=(2, 0.1),
        )

    with pytest.raises(ValueError, match="2-length tuple"):
        _distribution.__get_box_positions(
            widths=0.5,
            groups=cast(Any, (1, 2, 3)),
        )

    with pytest.raises(ValueError, match="expected None or 2-length tuple"):
        _distribution.__get_box_positions(
            widths=0.5,
            groups=cast(Any, "bad"),
        )

    with pytest.raises(ValueError, match="2-length tuple"):
        _distribution.__get_box_positions(
            widths=0.5,
            groups=(1, 0.2),
            hues=cast(Any, (1, 2, 3)),
        )

    with pytest.raises(ValueError, match="expected None or 2-length tuple"):
        _distribution.__get_box_positions(
            widths=0.5,
            groups=(1, 0.2),
            hues=cast(Any, "bad"),
        )


@pytest.mark.parametrize("deprecated", ["showlegend", "show_legend"])
@pytest.mark.parametrize("value", [False, True])
def test_distribution_deprecates_show_legend_without_effect(
    mini_adata,
    deprecated,
    value,
):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")
    kwargs = {deprecated: value}

    with pytest.warns(DeprecationWarning, match=f"'{deprecated}' is deprecated"):
        fig, ax, _ = bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            **kwargs,
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(DeprecationWarning, match=f"'{deprecated}' is deprecated"):
        fig, ax, _ = bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            legend=False,
            **kwargs,
        )

    assert ax.get_legend() is None
    plt.close(fig)


@pytest.mark.parametrize(
    ("deprecated", "value"),
    [
        ("show_median", False),
        ("show_medians", False),
        ("showmedians", False),
        ("show_mean", True),
        ("show_means", True),
        ("showmeans", True),
        ("showcaps", False),
        ("showbox", False),
        ("showfliers", False),
        ("showpoints", True),
    ],
)
def test_distribution_deprecates_legacy_show_kwargs(mini_adata, deprecated, value):
    with pytest.warns(FutureWarning, match=f"`{deprecated}` is deprecated"):
        fig, _, _ = bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            **{deprecated: value},
        )

    plt.close(fig)


@pytest.mark.parametrize("deprecated", ["lgd_params", "legend_params"])
def test_distribution_deprecates_legacy_legend_kwargs(mini_adata, deprecated):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")
    kwargs: Dict[str, Any] = {deprecated: {"loc": "upper left"}}

    with pytest.warns(FutureWarning, match=f"`{deprecated}` is deprecated"):
        result = bt.sct.pl.distribution(
            mini_adata,
            obs="score",
            groupby="cluster",
            hue="condition",
            **kwargs,
        )

    assert result is not None
    fig, ax, _ = result
    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(FutureWarning, match=f"`{deprecated}` is deprecated"):
        with pytest.raises(TypeError, match=f"{deprecated}.*legend"):
            bt.sct.pl.distribution(
                mini_adata,
                obs="score",
                groupby="cluster",
                hue="condition",
                legend={"loc": "upper left"},
                **kwargs,
            )
