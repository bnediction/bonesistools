#!/usr/bin/env python

from collections.abc import Mapping
from typing import (
    Optional,
    Union,
    Sequence,
    Tuple,
    List,
    Literal,
    Any
)
from ._typing import RGB
from .._typing import (
    ScData,
    anndata_or_mudata_checker
)

from pathlib import Path

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.lines import Line2D
from matplotlib.colors import Colormap
from itertools import cycle
Colors = Union[Sequence[RGB], cycle, Colormap]

from ._colors import (
    generate_colormap,
    QUALITATIVE_COLORS,
    black,
    gray
)

BoxItem = Literal["whiskers", "caps", "boxes", "medians", "fliers", "means"]
BoxPlots = Mapping[BoxItem, List[Line2D]]

def __get_box_positions(
    widths: float,
    groups: Optional[Tuple[int, float]] = None,
    hues: Optional[Tuple[int, float]] = None,
) -> np.ndarray:

    if groups is None and hues is not None:
        raise ValueError(f"invalid argument value for 'groups' and 'hues': expected either both specified or only 'groupbs' specified")

    if not groups is None:
        if isinstance(groups, Tuple):
            if not len(groups) == 2:
                raise ValueError(f"invalid argument value for 'groups': expected 2-length tuple but received {len(groups)}-length tuple")
        else:
            raise ValueError(f"invalid argument value for 'groups': expected {None} or 2-length tuple but received {groups}")

    if not hues is None:
        if isinstance(hues, Tuple):
            if not len(hues) == 2:
                raise ValueError(f"invalid argument value for 'hues': expected 2-length tuple but received {len(hues)}-length tuple")
        else:
            raise ValueError(f"invalid argument value for 'hues': expected {None} or 2-length tuple but received {hues}")
    
    if groups is None and hues is None:
        return np.zeros(shape=(1,))
    elif groups is not None and hues is None:
        return np.array(np.arange(0, groups[0]) * (widths + groups[1]))
    else:
        within_widths = widths + hues[1]
        within_positions = np.array(np.arange(0, hues[0]) * within_widths)
        between_positions = within_widths + widths + groups[1]
        positions = np.tile(within_positions, (groups[0], 1))
        for i, row in enumerate(positions):
            row[...] = row + between_positions * i
        return positions

def __apply_box_colors(
    bp,
    color,
    items: BoxItem,
):
    for k in bp.keys():
        if k in items:
            plt.setp(bp.get(k), color=color)
    
    return None

@anndata_or_mudata_checker
def boxplot(
    scdata: ScData, # type: ignore
    obs: str,
    groupby: Optional[str] = None,
    hue: Optional[str] = None,
    notch: Optional[bool] = None,
    sym: Optional[str] = None,
    patch_artist: Optional[bool] = None,
    vert: Optional[bool] = None,
    title: Optional[Union[str, dict]] = None,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colormap] = None,
    apply_colors_to: Tuple[BoxItem, ...] = ["whiskers", "caps", "boxes"],
    showmedians: bool = True,
    showmeans: bool = False,
    showcaps: bool = True,
    showbox: bool = True,
    showfliers: Optional[bool] = None,
    showpoints: Optional[bool] = None,
    showlegend: bool = True,
    outfile: Optional[Path] = None,
    **kwargs: Mapping[str, Any]
) -> Tuple[Figure, Axes, BoxPlots]:

    fig = plt.figure()
    ax = fig.subplots()
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
    fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else 5)

    if title:
        if isinstance(title, str):
            fig.canvas.manager.set_window_title(title)
            ax.set_title(title)
        elif isinstance(title, dict):
            fig.canvas.manager.set_window_title(title["label"])
            ax.set_title(**title)
        else:
            raise TypeError(f"unsupported argument type for 'title': expected {str} or {dict}, but received {type(title)}")

    if not "boxplot" in kwargs:
        kwargs["boxplot"] = {}
    
    if showfliers is None:
        showfliers = True if showpoints is not True else False
    if showpoints is None:
        showpoints = True if showfliers is not True else False

    if groupby is None and hue is None:
        groups = None
        hues = None
        series = scdata.obs[obs]
    elif groupby is not None and hue is None:
        groups = scdata.obs[groupby].cat.categories
        hues = None
        series = scdata.obs.groupby(by=groupby)[obs]
    elif groupby is None and hue is not None:
        raise ValueError(f"invalid argument value for 'groupby' and 'hue': expected either both specified or only 'groupby' specified")
    else:
        groups = scdata.obs[groupby].cat.categories
        hues = scdata.obs[hue].cat.categories
        series = [scdata.obs[scdata.obs[hue] == cat].groupby(by=[groupby])[obs] for cat in hues]

    positions = __get_box_positions(
        widths=widths,
        groups=None if groups is None else (len(groups), groupby_spacing),
        hues=None if hues is None else (len(hues), hue_spacing)
    )

    bps = []

    if hue is None:
        bps.append(
            plt.boxplot(
                x=series.dropna() if groupby is None else [series.get_group(group).dropna() for group in groups],
                positions=positions,
                notch=notch,
                sym=sym,
                vert=vert,
                widths=widths,
                patch_artist=patch_artist,
                showmeans=showmeans,
                showcaps=showcaps,
                showbox=showbox,
                showfliers=showfliers,
                **kwargs["boxplot"]
            )
        )
    else:
        if colors is None:
            if len(QUALITATIVE_COLORS) >= len(hues):
                colors = QUALITATIVE_COLORS[0:len(hues)]
            else:
                colors = generate_colormap(color_number=len(hues))
        elif isinstance(colors, Mapping):
            colors = [colors[_] for _ in hues]
        if hasattr(colors, "colors"):
            colors = colors.colors

        if apply_colors_to:
            if not "medianprops" in kwargs["boxplot"]:
                kwargs["boxplot"]["medianprops"] = {}
            if not "color" in kwargs["boxplot"]["medianprops"]:
                kwargs["boxplot"]["medianprops"]["color"] = black

        series_iterator = iter(series)
        for pos in positions.transpose():
            hue_series = next(series_iterator)
            bps.append(
                plt.boxplot(
                    x=[hue_series.get_group(group).dropna() for group in hue_series.groups.keys()],
                    positions=pos,
                    notch=notch,
                    sym=sym,
                    vert=vert,
                    widths=widths,
                    patch_artist=patch_artist,
                    showmeans=showmeans,
                    showcaps=showcaps,
                    showbox=showbox,
                    showfliers=showfliers,
                    **kwargs["boxplot"]
                ))

        if apply_colors_to:
            for i, bp in enumerate(bps):
                __apply_box_colors(
                    bp=bp,
                    color=colors[i],
                    items=apply_colors_to
                )
                plt.plot([], c=colors[i], label=hues[i])
            if showlegend is True:
                plt.legend()

    if showmedians is False:
        for bp in bps:
            for median in bp["medians"]:
                median.set(linewidth=0)

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax, bps
