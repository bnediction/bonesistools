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
from pandas import Series
from pandas.core.groupby.generic import SeriesGroupBy
from ._typing import RGB
from .._typing import (
    ScData,
    anndata_or_mudata_checker
)

from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.lines import Line2D
from matplotlib.colors import Colormap
from itertools import cycle
from ._colors import (
    generate_colormap,
    QUALITATIVE_COLORS,
    gray,
    black,
    white
)

Colors = Union[Sequence[RGB], cycle, Colormap]
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
        return positions.transpose()

def __apply_box_colors(
    bp,
    color,
    items: BoxItem,
):
    for k in bp.keys():
        if k in items:
            plt.setp(bp.get(k), color=color)
    
    return None

def __add_points(
    series: Union[Series, SeriesGroupBy, Mapping[str, SeriesGroupBy]],
    positions: np.ndarray,
    scale: float,
    groups: Optional[Sequence] = None,
    hues: Optional[Sequence] = None,
    **kwargs: Mapping[str, Any]
):
    ax = plt.gca()
    if groups is None:
        pos = positions[0]
        y = series.dropna()
        x = np.random.normal(pos, scale=scale, size=len(y))
        ax.scatter(x, y, **kwargs)
    elif groups is not None and hues is None:
        for i, group in enumerate(groups):
            pos = positions[i]
            y = series.get_group(group).dropna()
            x = np.random.normal(pos, scale=scale, size=len(y))
            ax.scatter(x, y, **kwargs)
    elif groups is None and hues is not None:
        raise ValueError(f"invalid argument value for 'groups' and 'hues': expected either both specified or only 'groups' specified")
    else:
        for i, hue in enumerate(hues):
            for j, group in enumerate(groups):
                pos = positions[i, j]
                y = series[hue].get_group(group).dropna()
                x = np.random.normal(pos, scale=scale, size=len(y))
                ax.scatter(x, y, **kwargs[hue])
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
    sort: Optional[Literal["ascending", "descending"]] = None,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    box_colors: Optional[Colormap] = None,
    point_colors: Optional[Colormap] = None,
    boxitems_to_color: Tuple[BoxItem, ...] = ["whiskers", "caps", "boxes"],
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
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 6)
    fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else 5 if groupby is None else 8)

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

    if groupby is None:
        if hue is None:
            groups = None
            hues = None
            series = scdata.obs[obs]
        else:
            raise ValueError(f"invalid argument value for 'groupby' and 'hue': expected either both specified or only 'groupby' specified")
    else:
        if sort is None:
            groups = scdata.obs[groupby].cat.categories
        elif sort in ["ascending", "descending"]:
            series = scdata.obs.groupby(by=groupby)[obs]
            groups = series.median().sort_values(ascending=(sort=="ascending")).index
        else:
            raise ValueError(f"invalid argument value for 'sort': expected 'ascending' or 'descending' but received {sort}")
        if hue is None:
            hues = None
        else:
            hues = scdata.obs[hue].cat.categories
            series = {h: scdata.obs[scdata.obs[hue] == h].groupby(by=[groupby])[obs] for h in hues}

    positions = __get_box_positions(
        widths=widths,
        groups=None if groups is None else (len(groups), groupby_spacing),
        hues=None if hues is None else (len(hues), hue_spacing)
    )

    if hue is None:
        if point_colors is None:
            point_colors = gray
        bps = plt.boxplot(
            x=series.dropna() if groupby is None else [series.get_group(g).dropna() for g in groups],
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
        if showmedians is False:
            for median in bps["medians"]:
                median.set(linewidth=0)
    else:
        bps = {}

        if box_colors is None:
            if showpoints is True:
                box_colors = [black]*len(hues)
            else:
                if len(QUALITATIVE_COLORS) >= len(hues):
                    box_colors = QUALITATIVE_COLORS[0:len(hues)]
                else:
                    box_colors = generate_colormap(color_number=len(hues))
        if hasattr(box_colors, "colors"):
            box_colors = box_colors.colors
        if not isinstance(box_colors, dict):
            box_colors = {h: box_colors[i] for i, h in enumerate(hues)}
        
        if point_colors is None:
            if showpoints is True:
                if len(QUALITATIVE_COLORS) >= len(hues):
                    point_colors = QUALITATIVE_COLORS[0:len(hues)]
                else:
                    point_colors = generate_colormap(color_number=len(hues))
            else:
                point_colors = [white]*len(hues)
        if hasattr(point_colors, "colors"):
            point_colors = point_colors.colors
        if not isinstance(point_colors, dict):
            point_colors = {h: point_colors[i] for i, h in enumerate(hues)}

        if not "medianprops" in kwargs["boxplot"]:
            kwargs["boxplot"]["medianprops"] = {}
        if "color" not in kwargs["boxplot"]["medianprops"] and "medians" not in boxitems_to_color:
            kwargs["boxplot"]["medianprops"]["color"] = black

        positions_iterator = iter(positions)
        for h, values in series.items():
            bps[h] = plt.boxplot(
                x=[values.get_group(g).dropna() for g in groups],
                positions=next(positions_iterator),
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
        
        if showmedians is False:
            for bp in bps.values():
                for median in bp["medians"]:
                    median.set(linewidth=0)

        if boxitems_to_color:
            for h in hues:
                __apply_box_colors(
                    bp=bps[h],
                    color=box_colors[h],
                    items=boxitems_to_color
                )

        if showlegend is True:
            handles = []
            for h in hues:
                handles.append(mpatches.Patch(edgecolor=box_colors[h], facecolor=point_colors[h], label=h))
            plt.legend(handles=handles)
    
    if showpoints:
        if not "scatter" in kwargs:
            kwargs["scatter"] = {}
        if hue is None:
            kwargs["scatter"].update(
                {
                    "s": 1 if not "s" in kwargs["scatter"] else kwargs["scatter"]["s"],
                    "alpha": 0.7 if not "alpha" in kwargs["scatter"] else kwargs["scatter"]["alpha"],
                    "facecolors": point_colors if not "facecolors" in kwargs["scatter"] else kwargs["scatter"]["facecolors"],
                    "edgecolors": "none" if not "edgecolors" in kwargs["scatter"] else kwargs["edgecolors"]["s"]
                }
            )
        else:
            for h in hues:
                if h not in kwargs["scatter"]:
                    kwargs["scatter"][h] = {}
                kwargs["scatter"][h].update(
                    {
                        "s": 1 if not "s" in kwargs["scatter"][h] else kwargs["scatter"][h]["s"],
                        "alpha": 0.7 if not "alpha" in kwargs["scatter"][h] else kwargs["scatter"][h]["alpha"],
                        "facecolors": point_colors[h],
                        "edgecolors": "none" if not "edgecolors" in kwargs["scatter"][h] else kwargs["edgecolors"][h]["s"]
                    }
                )
        __add_points(
            series=series,
            positions=positions,
            scale=widths/8,
            groups=groups,
            hues=hues,
            **kwargs["scatter"]
        )

    ylim_min, ylim_max = ax.get_ylim()
    ylim_diff = ylim_max - ylim_min
    plt.ylim(ylim_min - ylim_diff/50, ylim_max + ylim_diff/50)
    
    if groupby is not None:
        xticks = positions.sum(axis=0)/2 if positions.ndim == 2 else positions
        plt.xticks(xticks, groups)
    else:
        plt.tick_params(
            axis="x",
            which="both",
            bottom=False,
            top=False,
            labelbottom=False
        )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax, bps
