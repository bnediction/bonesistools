#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import bonesistools as bt
from bonesistools.sctools.plotting import _barplot


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

    fig, ax, table = bt.sct.pl.composition(
        mini_adata,
        obs="kind",
        groupby="state",
        group_order=["G0", "G1"],
        obs_order=["alpha", "beta"],
    )

    expected = pd.DataFrame(
        {
            "alpha": [0.5, 1.0],
            "beta": [0.5, 0.0],
        },
        index=pd.Index(["G0", "G1"], name="state"),
        columns=pd.Index(["alpha", "beta"], name="kind"),
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    pd.testing.assert_frame_equal(table, expected)
    assert ax.get_ylim() == (0.0, 1.0)
    assert ax.get_legend() is not None
    assert ax.get_xticklabels()[0].get_color() == "tab:green"
    assert ax.get_xticklabels()[0].get_fontweight() == "bold"
    plt.close(fig)


def test_composition_can_plot_counts_without_legend(mini_adata):
    fig, ax, table = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        normalize=False,
        showlegend=False,
        colors={"b1": "red", "b2": "blue"},
    )

    expected = pd.DataFrame(
        {
            "b1": [1, 1],
            "b2": [1, 1],
        },
        index=pd.Index(["A", "B"], name="cluster"),
        columns=pd.Index(["b1", "b2"], name="batch"),
    )

    pd.testing.assert_frame_equal(table, expected)
    assert ax.get_legend() is None
    plt.close(fig)


def test_composition_horizontal_orientation(mini_adata):
    fig, ax, table = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        orientation="horizontal",
        xlabel="proportion",
        ylabel="cluster",
        group_colors={"A": "tab:red", "B": "tab:blue"},
    )

    assert table.loc["A", "b1"] == 0.5
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

    fig, ax, table = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
    )

    assert list(table.columns) == ["b1", "b2"]
    assert len(ax.containers) == 2
    plt.close(fig)


def test_composition_accepts_colormap_name(mini_adata):
    fig, ax, table = bt.sct.pl.composition(
        mini_adata,
        obs="batch",
        groupby="cluster",
        colors="tab20",
    )

    assert list(table.columns) == ["b1", "b2"]
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
            title=object(),
        )

    with pytest.raises(ValueError, match="orientation"):
        bt.sct.pl.composition(
            mini_adata,
            obs="batch",
            groupby="cluster",
            orientation="diagonal",
        )
