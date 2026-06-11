#!/usr/bin/env python

from __future__ import annotations

import warnings
from collections.abc import Mapping as MappingABC
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure
from pandas import Series
from pandas.core.groupby.generic import SeriesGroupBy

from ..._compat import Literal
from .._typing import ScData, anndata_or_mudata_checker
from ._colors import (
    QUALITATIVE_COLORS,
    black,
    generate_colormap,
    gray,
    white,
)
from ._utils import colormap_colors, qualitative_color_values, set_window_title

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]
BoxItem = Literal["whiskers", "caps", "boxes", "medians", "fliers", "means"]
BoxPlots = Mapping[str, Any]
BoxplotReturn = Union[BoxPlots, Dict[object, BoxPlots]]


def __get_box_positions(
    widths: float,
    groups: Optional[Tuple[int, float]] = None,
    hues: Optional[Tuple[int, float]] = None,
) -> np.ndarray:

    if groups is None and hues is not None:
        raise ValueError(
            "invalid argument values for 'groups' and 'hues': "
            "expected either both specified or only 'groups' specified"
        )

    if groups is not None:
        if isinstance(groups, tuple):
            if not len(groups) == 2:
                raise ValueError(
                    f"invalid argument value for 'groups': "
                    f"expected 2-length tuple but received {len(groups)}-length tuple"
                )
        else:
            raise ValueError(
                f"invalid argument value for 'groups': "
                f"expected None or 2-length tuple but received {groups!r}"
            )

    if hues is not None:
        if isinstance(hues, tuple):
            if not len(hues) == 2:
                raise ValueError(
                    f"invalid argument value for 'hues': "
                    f"expected 2-length tuple but received {len(hues)}-length tuple"
                )
        else:
            raise ValueError(
                f"invalid argument value for 'hues': "
                f"expected None or 2-length tuple but received {hues!r}"
            )

    if groups is None and hues is None:
        return np.zeros(shape=(1,))
    elif groups is not None and hues is None:
        return np.array(np.arange(0, groups[0]) * (widths + groups[1]))
    else:
        assert groups is not None
        assert hues is not None
        within_widths = widths + hues[1]
        within_positions = np.array(np.arange(0, hues[0]) * within_widths)
        between_positions = within_widths + widths + groups[1]
        positions = cast(np.ndarray, np.tile(within_positions, (groups[0], 1)))
        for i, row in enumerate(positions):
            row[...] = row + between_positions * i
        return positions.transpose()


def __apply_box_colors(
    bp: BoxPlots,
    color: object,
    items: Sequence[BoxItem],
) -> None:
    for k in bp.keys():
        if k in items:
            plt.setp(bp.get(k), color=color)

    return None


def __add_points(
    series: Union[Series, SeriesGroupBy, Mapping[str, SeriesGroupBy]],
    positions: np.ndarray,
    scale: float,
    groups: Optional[Sequence[object]] = None,
    hues: Optional[Sequence[object]] = None,
    **kwargs: Any,
) -> None:
    ax = plt.gca()
    if groups is None:
        pos = positions[0]
        series_values = cast(Series, series)
        y = series_values.dropna()
        x = np.random.normal(pos, scale=scale, size=len(y))
        ax.scatter(x, y, **kwargs)
    elif groups is not None and hues is None:
        grouped_series = cast(SeriesGroupBy, series)

        for i, group in enumerate(groups):
            pos = positions[i]
            y = grouped_series.get_group(group).dropna()
            x = np.random.normal(pos, scale=scale, size=len(y))
            ax.scatter(x, y, **kwargs)
    elif groups is None and hues is not None:
        raise ValueError(
            "invalid argument values for 'groups' and 'hues': "
            "expected either both specified or only 'groups' specified"
        )
    else:
        grouped_by_hue = cast(Mapping[object, SeriesGroupBy], series)
        group_values = cast(Sequence[object], groups)
        hue_values = cast(Sequence[object], hues)

        for i, hue in enumerate(hue_values):
            for j, group in enumerate(group_values):
                pos = positions[i, j]
                y = grouped_by_hue[hue].get_group(group).dropna()
                x = np.random.normal(pos, scale=scale, size=len(y))
                ax.scatter(x, y, **cast(Any, kwargs)[hue])
    return None


@anndata_or_mudata_checker
def distribution(
    scdata: ScData,  # type: ignore
    obs: str,
    groupby: Optional[str] = None,
    hue: Optional[str] = None,
    notch: Optional[bool] = None,
    sym: Optional[str] = None,
    patch_artist: Optional[bool] = None,
    vert: Optional[bool] = None,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    sort: Optional[Literal["ascending", "descending"]] = None,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    box_colors: Optional[Colors] = None,
    point_colors: Optional[Colors] = None,
    boxitems_to_color: Tuple[BoxItem, ...] = ("whiskers", "caps", "boxes"),
    showmedians: bool = True,
    showmeans: bool = False,
    showcaps: bool = True,
    showbox: bool = True,
    showfliers: Optional[bool] = None,
    showpoints: Optional[bool] = None,
    showlegend: bool = True,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes, BoxplotReturn]]:
    """
    Draw a distribution plot for an observation-level variable.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Observation key to plot.
    groupby: str, optional
        Observation key used to group values.
    hue: str, optional
        Observation key used to split each group by hue.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    sort: {"ascending", "descending"}, optional
        Sort groups by median value.
    outfile: Path, optional
        If specified, save the figure instead of returning it.

    Returns
    -------
    tuple
        Matplotlib figure, axes and boxplot artists.

    Raises
    ------
    TypeError
        If `title` is neither a string nor a dictionary.
    ValueError
        If `groupby` and `hue` are inconsistently specified, or if `sort` is
        not `"ascending"` or `"descending"`.
    """

    fig = plt.figure()
    ax = fig.subplots()
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 6)
    fig.set_figwidth(
        kwargs["figwidth"] if "figwidth" in kwargs else 5 if groupby is None else 8
    )

    if title:
        if isinstance(title, str):
            set_window_title(fig, title)
            ax.set_title(title)
        elif isinstance(title, dict):
            set_window_title(fig, title["label"])
            ax.set_title(**title)
        else:
            raise TypeError(
                f"unsupported argument type for 'title': "
                f"expected {str} or {dict} but received {type(title)}"
            )

    if "boxplot" not in kwargs:
        kwargs["boxplot"] = {}

    if showfliers is None:
        showfliers = True if showpoints is not True else False
    if showpoints is None:
        showpoints = True if showfliers is not True else False

    if groupby is None:
        if hue is None:
            groups = None
            hues = None
            series = cast(Series, scdata.obs[obs])
        else:
            raise ValueError(
                "invalid argument values for 'groupby' and 'hue': "
                "expected either both specified or only 'groupby' specified"
            )
    else:
        if sort is None:
            groups = list(scdata.obs[groupby].cat.categories)
            series = cast(
                SeriesGroupBy, scdata.obs.groupby(by=groupby, observed=True)[obs]
            )
        elif sort in ["ascending", "descending"]:
            series = cast(
                SeriesGroupBy, scdata.obs.groupby(by=groupby, observed=True)[obs]
            )
            groups = list(
                series.median().sort_values(ascending=(sort == "ascending")).index
            )
        else:
            raise ValueError(
                f"invalid argument value for 'sort': "
                f"expected 'ascending' or 'descending' but received {sort!r}"
            )
        if hue is None:
            hues = None
        else:
            hues = list(scdata.obs[hue].cat.categories)
            series = {
                h: cast(
                    SeriesGroupBy,
                    scdata.obs[scdata.obs[hue] == h].groupby(
                        by=groupby,
                        observed=True,
                    )[obs],
                )
                for h in hues
            }

    positions = __get_box_positions(
        widths=widths,
        groups=None if groups is None else (len(groups), groupby_spacing),
        hues=None if hues is None else (len(hues), hue_spacing),
    )

    if hue is None:
        if point_colors is None:
            point_colors = gray
        bps = cast(
            BoxPlots,
            plt.boxplot(
                x=(
                    cast(Series, series).dropna()
                    if groupby is None
                    else [
                        cast(SeriesGroupBy, series).get_group(g).dropna()
                        for g in cast(Sequence[object], groups)
                    ]
                ),
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
                **kwargs["boxplot"],
            ),
        )
        if showmedians is False:
            for median in cast(Sequence[Any], bps["medians"]):
                median.set(linewidth=0)
    else:
        grouped_bps: Dict[object, BoxPlots] = {}
        hue_values = cast(Sequence[object], hues)

        if box_colors is None:
            if showpoints is True:
                box_colors = [black] * len(hue_values)
            else:
                box_colors = qualitative_color_values(
                    len(hue_values),
                    QUALITATIVE_COLORS,
                    generate_colormap,
                )
        if isinstance(box_colors, ListedColormap):
            box_colors = colormap_colors(box_colors)
        if not isinstance(box_colors, MappingABC):
            box_colors = {
                h: cast(Sequence[object], box_colors)[i]
                for i, h in enumerate(hue_values)
            }

        if point_colors is None:
            if showpoints is True:
                point_colors = qualitative_color_values(
                    len(hue_values),
                    QUALITATIVE_COLORS,
                    generate_colormap,
                )
            else:
                point_colors = [white] * len(hue_values)
        if isinstance(point_colors, ListedColormap):
            point_colors = colormap_colors(point_colors)
        if not isinstance(point_colors, MappingABC):
            point_colors = {
                h: cast(Sequence[object], point_colors)[i]
                for i, h in enumerate(hue_values)
            }

        if "medianprops" not in kwargs["boxplot"]:
            kwargs["boxplot"]["medianprops"] = {}
        if (
            "color" not in kwargs["boxplot"]["medianprops"]
            and "medians" not in boxitems_to_color
        ):
            kwargs["boxplot"]["medianprops"]["color"] = black

        positions_iterator = iter(positions)
        grouped_series_by_hue = cast(Mapping[object, SeriesGroupBy], series)

        for h, values in grouped_series_by_hue.items():
            grouped_bps[h] = cast(
                BoxPlots,
                plt.boxplot(
                    x=[
                        values.get_group(g).dropna()
                        for g in cast(Sequence[object], groups)
                    ],
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
                    **kwargs["boxplot"],
                ),
            )

        if showmedians is False:
            for bp in grouped_bps.values():
                for median in cast(Sequence[Any], bp["medians"]):
                    median.set(linewidth=0)

        if boxitems_to_color:
            for h in hue_values:
                __apply_box_colors(
                    bp=grouped_bps[h],
                    color=cast(Mapping[object, object], box_colors)[h],
                    items=boxitems_to_color,
                )

        if showlegend is True:
            handles = []
            for h in hue_values:
                handles.append(
                    mpatches.Patch(
                        edgecolor=cast(
                            Any,
                            cast(Mapping[object, object], box_colors)[h],
                        ),
                        facecolor=cast(
                            Any,
                            cast(Mapping[object, object], point_colors)[h],
                        ),
                        label=h,
                    )
                )
            plt.legend(handles=handles)

        bps = grouped_bps

    if showpoints:
        if "scatter" not in kwargs:
            kwargs["scatter"] = {}
        if hue is None:
            kwargs["scatter"].update(
                {
                    "s": 1 if "s" not in kwargs["scatter"] else kwargs["scatter"]["s"],
                    "alpha": (
                        0.7
                        if "alpha" not in kwargs["scatter"]
                        else kwargs["scatter"]["alpha"]
                    ),
                    "facecolors": (
                        point_colors
                        if "facecolors" not in kwargs["scatter"]
                        else kwargs["scatter"]["facecolors"]
                    ),
                    "edgecolors": (
                        "none"
                        if "edgecolors" not in kwargs["scatter"]
                        else kwargs["edgecolors"]["s"]
                    ),
                }
            )
        else:
            for h in cast(Sequence[object], hues):
                if h not in kwargs["scatter"]:
                    kwargs["scatter"][h] = {}
                kwargs["scatter"][h].update(
                    {
                        "s": (
                            1
                            if "s" not in kwargs["scatter"][h]
                            else kwargs["scatter"][h]["s"]
                        ),
                        "alpha": (
                            0.7
                            if "alpha" not in kwargs["scatter"][h]
                            else kwargs["scatter"][h]["alpha"]
                        ),
                        "facecolors": cast(Mapping[object, object], point_colors)[h],
                        "edgecolors": (
                            "none"
                            if "edgecolors" not in kwargs["scatter"][h]
                            else kwargs["edgecolors"][h]["s"]
                        ),
                    }
                )
        __add_points(
            series=series,
            positions=positions,
            scale=widths / 8,
            groups=groups,
            hues=hues,
            **kwargs["scatter"],
        )

    ylim_min, ylim_max = ax.get_ylim()
    ylim_diff = ylim_max - ylim_min
    plt.ylim(ylim_min - ylim_diff / 50, ylim_max + ylim_diff / 50)

    if groupby is not None:
        xticks = positions.sum(axis=0) / 2 if positions.ndim == 2 else positions
        plt.xticks(xticks, [str(group) for group in cast(Sequence[object], groups)])
    else:
        plt.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax, bps


def boxplot(*args: Any, **kwargs: Any) -> Optional[Tuple[Figure, Axes, BoxplotReturn]]:
    """
    Deprecated alias for `distribution`.
    """

    warnings.warn(
        "`bt.sct.pl.boxplot` is deprecated; use " "`bt.sct.pl.distribution` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return distribution(*args, **kwargs)
