#!/usr/bin/env python

import cycler
import matplotlib as mpl, matplotlib.pyplot as plt

from typing import Optional
from matplotlib.axes._axes import Axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap
from . import _colors

def set_default_params():

    mpl.rcParams.update(mpl.rcParamsDefault)

    font = {"family" : "normal",
            "weight" : "normal",
            "size"   : 12}
    mpl.rc("font", **font)

    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["lines.linewidth"] = 1.5

    mpl.rcParams.update({
        "axes.spines.top"    : False,
        "axes.spines.bottom" : True,
        "axes.spines.left"   : True,
        "axes.spines.right"  : False
    })

    __margin = 0
    mpl.rcParams["axes.xmargin"] = __margin
    mpl.rcParams["axes.ymargin"] = __margin
    mpl.rcParams["axes.zmargin"] = __margin
    mpl.rcParams["axes.labelsize"] = 14

    mpl.rcParams["axes.prop_cycle"] = cycler.cycler(color=[
        _colors.blue,
        _colors.red,
        _colors.green,
        _colors.orange,
        _colors.purple,
        _colors.skyblue,
        _colors.teal,
        _colors.pink,
        _colors.violet,
        _colors.darkblue
    ])

    return None

cmap = ListedColormap(
    colors  = _colors.COLORS,
    name    = "default",
    N       = None
)

mpl.colormaps.register(cmap)

def set_default_axis(
    ax: Optional[Axes] = None
):
    if ax is None:
        ax = plt.gca()
    ax.set_title("")
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%g"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    return None
