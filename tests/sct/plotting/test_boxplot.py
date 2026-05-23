#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import pytest

import bonesistools as bt
from bonesistools.sctools.plotting import _boxplot

ADATA = bt.sct.datasets.nestorowa()


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

    adata = ADATA.copy()
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


def test_boxplot_without_groupby_and_validation_errors(mini_adata):
    fig, ax, bps = bt.sct.pl.boxplot(
        mini_adata,
        obs="score",
        title={"label": "score"},
        showpoints=True,
        showmedians=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_title() == "score"
    assert len(bps["boxes"]) == 1
    assert all(median.get_linewidth() == 0 for median in bps["medians"])
    plt.close(fig)

    with pytest.raises(ValueError, match="invalid argument values"):
        bt.sct.pl.boxplot(mini_adata, obs="score", hue="cluster")

    with pytest.raises(TypeError, match="unsupported argument type for 'title'"):
        bt.sct.pl.boxplot(mini_adata, obs="score", title=object())


def test_boxplot_with_hue_custom_colors_and_hidden_medians(mini_adata):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, bps = bt.sct.pl.boxplot(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        box_colors={"ctrl": [1.0, 0.0, 0.0], "stim": [0.0, 0.0, 1.0]},
        point_colors={"ctrl": [1.0, 0.5, 0.5], "stim": [0.5, 0.5, 1.0]},
        showpoints=True,
        showmedians=False,
        showlegend=True,
    )

    assert isinstance(ax, Axes)
    assert set(bps) == {"ctrl", "stim"}
    assert ax.get_legend() is not None
    for bp in bps.values():
        assert all(median.get_linewidth() == 0 for median in bp["medians"])
    plt.close(fig)


def test_boxplot_position_and_point_helper_validation():
    with pytest.raises(ValueError, match="groups' and 'hues'"):
        _boxplot.__get_box_positions(widths=0.5, hues=(2, 0.1))

    with pytest.raises(ValueError, match="2-length tuple"):
        _boxplot.__get_box_positions(widths=0.5, groups=(1, 2, 3))

    with pytest.raises(ValueError, match="expected None or 2-length tuple"):
        _boxplot.__get_box_positions(widths=0.5, groups="bad")

    with pytest.raises(ValueError, match="2-length tuple"):
        _boxplot.__get_box_positions(widths=0.5, groups=(1, 0.2), hues=(1, 2, 3))

    with pytest.raises(ValueError, match="expected None or 2-length tuple"):
        _boxplot.__get_box_positions(widths=0.5, groups=(1, 0.2), hues="bad")
