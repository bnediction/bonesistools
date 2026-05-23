#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt

import bonesistools as bt


def test_set_default_params_and_axis_update_matplotlib_state():
    bt.sct.pl.set_default_params()

    assert mpl.rcParams["text.usetex"] is True
    assert mpl.rcParams["axes.spines.top"] is False

    fig, ax = plt.subplots()
    ax.set_title("temporary")

    result = bt.sct.pl.set_default_axis(ax)

    assert result is None
    assert ax.get_title() == ""
    assert ax.xaxis.get_major_formatter().format_data(1.5) == "1.5"
    assert ax.yaxis.get_major_formatter().format_data(1.5) == "1.5"
    plt.close(fig)
