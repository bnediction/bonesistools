#!/usr/bin/env python

from typing import Optional

import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes._axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter

from ._colors import (
    CLASSIC_COLORS,
    black,
    blue,
    darkred,
    green,
    maroon,
    navy,
    orange,
    pink,
    purple,
    red,
    skyblue,
    teal,
    violet,
)


def set_default_params(tex: bool = True):
    """
    Set default parameters for matplotlib.

    Parameters
    ----------
    tex: bool (default: True)
        If True, enable LaTeX rendering for text. Set to False in
        environments where LaTeX is not installed.
    """

    mpl.rcParams.update(mpl.rcParamsDefault)

    font = {"weight": "normal", "size": 12}
    mpl.rc("font", **font)

    mpl.rcParams["text.usetex"] = tex
    mpl.rcParams["lines.linewidth"] = 1.5
    mpl.rcParams["legend.fontsize"] = 12
    mpl.rcParams["legend.markerscale"] = 3

    mpl.rc(
        "axes",
        **{
            "spines.top": False,
            "spines.bottom": True,
            "spines.left": True,
            "spines.right": False,
            "xmargin": 0,
            "ymargin": 0,
            "zmargin": 0,
            "labelsize": 14,
        },
    )

    mpl.rcParams["axes.prop_cycle"] = cycler.cycler(
        color=[
            blue,
            red,
            green,
            orange,
            purple,
            pink,
            skyblue,
            teal,
            violet,
            navy,
            darkred,
            maroon,
        ]
    )

    mpl.rc(
        "boxplot",
        **{
            "notch": False,
            "vertical": True,
            "whiskers": 2,
            "showmeans": False,
            "showcaps": True,
            "showbox": True,
            "showfliers": True,
            "flierprops.markerfacecolor": black,
            "flierprops.linewidth": 2.0,
            "whiskerprops.color": black,
            "whiskerprops.linewidth": 2.0,
            "whiskerprops.linestyle": ":",
            "capprops.color": black,
            "capprops.linewidth": 2.0,
            "medianprops.color": blue,
            "medianprops.linewidth": 2.0,
            "meanprops.color": blue,
            "meanprops.linewidth": 2.0,
            "meanprops.linestyle": "-",
        },
    )

    return None


cmap = ListedColormap(colors=CLASSIC_COLORS, name="default", N=None)

mpl.colormaps.register(cmap)


def set_default_axis(ax: Optional[Axes] = None):
    """
    Set default parameters for matplotlib axes.
    """

    if ax is None:
        ax = plt.gca()
    ax.set_title("")
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%g"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    return None
