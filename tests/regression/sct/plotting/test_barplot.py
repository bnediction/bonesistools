#!/usr/bin/env python

from typing import Any, Dict, cast

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.container import BarContainer
from matplotlib.figure import Figure

import bonesistools as bt
from bonesistools.sctools.plotting import _barplot


def _bar_container(ax: Axes, index: int) -> BarContainer:
    return cast(BarContainer, ax.containers[index])


def test_composition_plots_expected_proportions(mini_adata):
    mini_adata.obs["state"] = pd.Categorical(
        ["G0", "G0", "G1", "G1"],
        categories=["G0", "G1", "G2"],
    )
    mini_adata.obs["kind"] = pd.Categorical(
        ["alpha", "beta", "alpha", "alpha"],
        categories=["alpha", "beta"],
    )
    mini_adata.uns["kind_color"] = {"alpha": "black", "beta": "tab:blue"}
    mini_adata.uns["state_color"] = {
        "G0": "tab:green",
        "G1": "tab:orange",
        "G2": "tab:purple",
    }

    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="kind",
        groupby="state",
        group_order=["G0", "G1"],
        obs_order=["alpha", "beta"],
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert [patch.get_height() for patch in _bar_container(ax, 0).patches] == [
        0.5,
        1.0,
    ]
    assert [patch.get_height() for patch in _bar_container(ax, 1).patches] == [
        0.5,
        0.0,
    ]
    assert ax.get_ylim() == (0.0, 1.0)
    assert ax.get_legend() is not None
    assert ax.get_xticklabels()[0].get_color() == "tab:green"
    assert ax.get_xticklabels()[0].get_fontweight() == "bold"
    plt.close(fig)


def test_composition_can_plot_counts_without_legend(mini_adata):
    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        normalize=False,
        legend=False,
        colors={"b1": "red", "b2": "blue"},
    )

    assert [patch.get_height() for patch in _bar_container(ax, 0).patches] == [1, 1]
    assert [patch.get_height() for patch in _bar_container(ax, 1).patches] == [1, 1]
    assert ax.get_legend() is None
    plt.close(fig)


def test_composition_legend_uses_rcparams_fontsize(mini_adata):
    with mpl.rc_context({"legend.fontsize": 7}):
        fig, ax = bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
        )

    legend = ax.get_legend()
    assert legend is not None
    assert [text.get_fontsize() for text in legend.get_texts()] == [7, 7]
    plt.close(fig)


def test_composition_accepts_explicit_legend_kwargs(mini_adata):
    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        legend={"title": "batch"},
    )

    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "batch"
    plt.close(fig)


def test_composition_horizontal_orientation(mini_adata):
    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        orientation="horizontal",
        xlabel="proportion",
        ylabel="cluster",
        group_colors={"A": "tab:red", "B": "tab:blue"},
    )

    assert _bar_container(ax, 0).patches[0].get_width() == 0.5
    assert ax.get_xlabel() == "proportion"
    assert ax.get_ylabel() == "cluster"
    assert ax.get_xlim() == (0.0, 1.0)
    assert ax.get_yticklabels()[0].get_color() == "tab:red"
    plt.close(fig)


def test_composition_uses_embedding_like_color_fallback(mini_adata, monkeypatch):
    monkeypatch.setattr(_barplot, "QUALITATIVE_COLORS", _barplot.QUALITATIVE_COLORS[:1])
    generated = ["red", "green"]

    def fake_colormap(color_number):
        assert color_number == 2
        return generated

    monkeypatch.setattr(_barplot, "generate_colormap", fake_colormap)

    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
    )

    assert len(ax.containers) == 2
    plt.close(fig)


def test_composition_accepts_colormap_name(mini_adata):
    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        colors="tab20",
    )

    assert len(ax.containers) == 2
    plt.close(fig)


def test_composition_accepts_registered_colormap(mini_adata):
    fig, ax = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        colors=bt.sct.pl.get_colormap("earth"),
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert len(ax.containers) == 2
    plt.close(fig)


def test_composition_outfile_and_validation(mini_adata, tmp_path):
    mpl.rcParams["text.usetex"] = False
    outfile = tmp_path / "composition.png"

    result = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()

    with pytest.raises(KeyError, match="missing colors"):
        bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            colors={"b1": "red"},
        )

    with pytest.raises(ValueError, match="empty groups"):
        bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            group_order=["A", "B", "C"],
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'title'"):
        bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            title=cast(Any, object()),
        )

    with pytest.raises(ValueError, match="orientation"):
        bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            orientation=cast(Any, "diagonal"),
        )


@pytest.mark.parametrize("deprecated", ["showlegend", "show_legend"])
@pytest.mark.parametrize("value", [False, True])
def test_composition_deprecates_show_legend_without_effect(
    mini_adata,
    deprecated,
    value,
):
    kwargs = {deprecated: value}

    with pytest.warns(DeprecationWarning, match=f"'{deprecated}' is deprecated"):
        fig, ax = bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            **kwargs,
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(DeprecationWarning, match=f"'{deprecated}' is deprecated"):
        fig, ax = bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            legend=False,
            **kwargs,
        )

    assert ax.get_legend() is None
    plt.close(fig)


@pytest.mark.parametrize("deprecated", ["lgd_params", "legend_params"])
def test_composition_deprecates_legacy_legend_kwargs(mini_adata, deprecated):
    kwargs: Dict[str, Any] = {deprecated: {"loc": "upper left"}}

    with pytest.warns(FutureWarning, match=f"`{deprecated}` is deprecated"):
        fig, ax = bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            **kwargs,
        )

    assert ax.get_legend() is not None
    plt.close(fig)

    with pytest.warns(FutureWarning, match=f"`{deprecated}` is deprecated"):
        with pytest.raises(TypeError, match=f"{deprecated}.*legend"):
            bt.sct.pl.composition(
                mini_adata,
                obs="batch",
                groupby="cluster",
                legend={"loc": "upper left"},
                **kwargs,
            )
