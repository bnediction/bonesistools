#!/usr/bin/env python

from typing import Any, Sequence, cast

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import bonesistools as bt
from bonesistools.sctools.plotting import _distribution

ADATA = bt.sct.datasets.nestorowa()


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
        showpoints=False,
        showlegend=True,
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
        showpoints=True,
        showmedians=False,
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
        assert all(median.get_linewidth() == 0 for median in _median_lines(bp))
    plt.close(fig)


def test_distribution_hue_defaults_and_listed_colormaps(mini_adata, monkeypatch):
    mini_adata.obs["condition"] = ["ctrl", "stim", "ctrl", "stim"]
    mini_adata.obs["condition"] = mini_adata.obs["condition"].astype("category")

    fig, ax, bps = bt.sct.pl.distribution(
        mini_adata,
        obs="score",
        groupby="cluster",
        hue="condition",
        title="score by condition",
        showpoints=True,
        showlegend=True,
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
        box_colors=box_colors,
        point_colors=point_colors,
        showpoints=False,
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)

    fig, _, bps = bt.sct.pl.distribution(
        hue_adata,
        obs="score",
        groupby="group",
        hue="hue4",
        showpoints=False,
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)

    fig, _, bps = bt.sct.pl.distribution(
        hue_adata,
        obs="score",
        groupby="group",
        hue="hue4",
        showpoints=True,
    )

    assert set(bps) == {"h1", "h2", "h3", "h4"}
    plt.close(fig)


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
