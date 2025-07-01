#!/usr/bin/env python

import cycler
import matplotlib as mpl, matplotlib.pyplot as plt

from typing import Optional
from matplotlib.axes._axes import Axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap
from . import _colors

def set_default_params():
    """
    Set default parameters for matplotlib.
    """

    mpl.rcParams.update(mpl.rcParamsDefault)

    font = {"family" : "normal",
            "weight" : "normal",
            "size"   : 12}
    mpl.rc("font", **font)

    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["lines.linewidth"] = 1.5

    mpl.rc(
        "axes",
        **{
            "spines.top"    : False,
            "spines.bottom" : True,
            "spines.left"   : True,
            "spines.right"  : False,
            "xmargin"       : 0,
            "ymargin"       : 0,
            "zmargin"       : 0,
            "labelsize"     : 14
        }
    )

    mpl.rcParams["axes.prop_cycle"] = cycler.cycler(color=[
        _colors.blue,
        _colors.red,
        _colors.green,
        _colors.orange,
        _colors.purple,
        _colors.pink,
        _colors.skyblue,
        _colors.teal,
        _colors.violet,
        _colors.navy,
        _colors.darkred,
        _colors.maroon
    ])

    mpl.rc(
        "boxplot",
        **{
            "notch"                         : False,
            "vertical"                      : True,
            "whiskers"                      : 2,
            "showmeans"                     : False,
            "showcaps"                      : True,
            "showbox"                       : True,
            "showfliers"                    : True,
            "flierprops.markerfacecolor"    : None,
            "flierprops.markerfacecolor"    : _colors.black,
            "flierprops.linewidth"          : 2.0,
            "whiskerprops.color"            : _colors.black,
            "whiskerprops.linewidth"        : 2.0,
            "whiskerprops.linestyle"        : ":",
            "capprops.color"                : _colors.black,
            "capprops.linewidth"            : 2.0,
            "medianprops.color"             : _colors.blue,
            "medianprops.linewidth"         : 2.0,
            "meanprops.color"               : _colors.blue,
            "meanprops.linewidth"           : 2.0,
            "meanprops.linestyle"           : "-"
        }
    )

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
