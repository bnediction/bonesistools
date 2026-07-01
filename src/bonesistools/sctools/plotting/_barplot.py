#!/usr/bin/env python

from __future__ import annotations

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
    overload,
)

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure

from ..._compat import Literal
from .._typing import ScData, anndata_or_mudata_checker
from ._colors import QUALITATIVE_COLORS, generate_colormap
from ._utils import (
    _resolve_legend_argument,
    colormap_colors,
    colors_from_uns,
    figure_from_axes,
    qualitative_color_values,
    set_window_title,
)

Colors = Union[
    str,
    Sequence[object],
    Iterator[object],
    Colormap,
    Mapping[object, object],
]
Orientation = Literal["vertical", "horizontal"]


def __ordered_values(
    series: pd.Series,
    order: Optional[Sequence[object]],
) -> Sequence[object]:

    if order is not None:
        return list(order)

    if isinstance(series.dtype, pd.CategoricalDtype):
        return list(series.cat.categories)

    return list(pd.unique(series.dropna()))


def __index_from_values(values: Sequence[object], name: str) -> pd.Index:

    return pd.Index(np.asarray(list(values), dtype=object), name=name)


def __normalize_colors(
    colors: Optional[Colors],
    values: Sequence[object],
    scdata: ScData,
    uns_key: str,
) -> Sequence[object]:

    if colors is None:
        colors = colors_from_uns(scdata, uns_key, values)

    if colors is None:
        colors = qualitative_color_values(
            len(values),
            QUALITATIVE_COLORS,
            generate_colormap,
        )

    if isinstance(colors, ListedColormap):
        colors = colormap_colors(colors)

    if isinstance(colors, str):
        colors = plt.get_cmap(colors)

    if isinstance(colors, Colormap):
        positions = np.linspace(0, 1, len(values))
        return cast(Sequence[object], colors(positions))

    if isinstance(colors, MappingABC):
        missing = [value for value in values if value not in colors]

        if missing:
            raise KeyError(f"missing colors for values: {missing!r}")

        return [colors[value] for value in values]

    color_values = list(cast(Sequence[object], colors))

    if len(color_values) < len(values):
        raise ValueError(
            f"invalid argument value for 'colors': "
            f"expected at least {len(values)} colors but received {len(color_values)}"
        )

    return color_values[0 : len(values)]


def __apply_ticklabel_colors(
    ax: Axes,
    colors: Optional[Mapping[object, object]],
    orientation: Orientation,
) -> None:

    if colors is None:
        return None

    ticklabels = (
        ax.get_xticklabels() if orientation == "vertical" else ax.get_yticklabels()
    )

    for ticklabel in ticklabels:
        value = ticklabel.get_text()

        if value in colors:
            ticklabel.set_color(cast(Any, colors[value]))
            ticklabel.set_fontweight("bold")

    return None


@overload
def composition(
    scdata: ScData,
    obs: str,
    groupby: str,
    *,
    obs_order: Optional[Sequence[object]] = None,
    group_order: Optional[Sequence[object]] = None,
    normalize: bool = True,
    percent: bool = True,
    dropna: bool = True,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    orientation: Orientation = "vertical",
    width: float = 0.8,
    colors: Optional[Colors] = None,
    group_colors: Optional[Mapping[object, object]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]: ...


@overload
def composition(
    scdata: ScData,
    obs: str,
    groupby: str,
    *,
    obs_order: Optional[Sequence[object]] = None,
    group_order: Optional[Sequence[object]] = None,
    normalize: bool = True,
    percent: bool = True,
    dropna: bool = True,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    orientation: Orientation = "vertical",
    width: float = 0.8,
    colors: Optional[Colors] = None,
    group_colors: Optional[Mapping[object, object]] = None,
    ax: Optional[Axes] = None,
    outfile: Path,
    **kwargs: Any,
) -> None: ...


@overload
def composition(
    scdata: ScData,
    obs: str,
    groupby: str,
    *,
    obs_order: Optional[Sequence[object]] = None,
    group_order: Optional[Sequence[object]] = None,
    normalize: bool = True,
    percent: bool = True,
    dropna: bool = True,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    orientation: Orientation = "vertical",
    width: float = 0.8,
    colors: Optional[Colors] = None,
    group_colors: Optional[Mapping[object, object]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]: ...


@anndata_or_mudata_checker
def composition(
    scdata: ScData,  # type: ignore
    obs: str,
    groupby: str,
    *,
    obs_order: Optional[Sequence[object]] = None,
    group_order: Optional[Sequence[object]] = None,
    normalize: bool = True,
    percent: bool = True,
    dropna: bool = True,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    orientation: Orientation = "vertical",
    width: float = 0.8,
    colors: Optional[Colors] = None,
    group_colors: Optional[Mapping[object, object]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw a stacked barplot of observation composition across groups.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Observation column defining stacked segments.
    groupby: str
        Observation column defining bars.
    obs_order: sequence, optional
        Segment order. If None, categorical order or observed order is used.
    group_order: sequence, optional
        Bar order. If None, categorical order or observed order is used.
    normalize: bool (default: True)
        Plot proportions within each group instead of raw counts.
    percent: bool (default: True)
        Format the y-axis as percentages when `normalize=True`.
    dropna: bool (default: True)
        Drop observations with missing `obs` or `groupby` values.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    legend: bool or mapping (default: True)
        Legend configuration. False disables the legend. True draws the legend
        using default Matplotlib parameters. If a mapping is provided, it is
        forwarded as keyword arguments to `Axes.legend`.
    show_legend: bool, optional
        Deprecated. This parameter has no effect and will be removed in
        bonesistools 2.0.0. Use `legend` instead.
    orientation: {"vertical", "horizontal"} (default: "vertical")
        Draw vertical stacked bars or horizontal stacked bars.
    width: float (default: 0.8)
        Bar width.
    colors: colormap name, sequence, colormap or mapping, optional
        Segment colors. If None, `scdata.uns["<obs>_color"]`,
        `scdata.uns["<obs>_colors"]` or a qualitative palette is used.
    group_colors: mapping, optional
        Colors applied to group tick labels. If None,
        `scdata.uns["<groupby>_color"]` or `scdata.uns["<groupby>_colors"]` is
        used when available.
    ax: Axes, optional
        Existing axes used for drawing.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - xlim[tuple]: set x-axis limits with `Axes.set_xlim`
        - ylim[tuple]: set y-axis limits with `Axes.set_ylim`
        - labelsize[float]: set the tick-label size when `tick_params` is not
          specified
        - rotation[float]: rotate group tick labels
        - tick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines following the syntax of matplotlib.axes.Axes.tick_params
        - bar[dict]: keyword arguments passed to `DataFrame.plot`

    Returns
    -------
    tuple or None
        Figure and axes if `outfile` is None; otherwise None.

    Raises
    ------
    KeyError
        If explicit color mappings miss plotted segment values.
    ValueError
        If `orientation` is invalid, if too few colors are provided, or if one
        group has no observations when normalized proportions are requested.
    """

    if orientation not in ("vertical", "horizontal"):
        raise ValueError(
            f"invalid argument value for 'orientation': "
            f"expected 'vertical' or 'horizontal' but received {orientation!r}"
        )

    draw_legend, legend_kwargs = _resolve_legend_argument(legend, kwargs)

    data = scdata.obs[[groupby, obs]]
    data = data.dropna() if dropna else data

    groups = __ordered_values(scdata.obs[groupby], group_order)
    segments = __ordered_values(scdata.obs[obs], obs_order)

    counts = pd.crosstab(data[groupby], data[obs])
    counts = counts.reindex(
        index=__index_from_values(groups, groupby),
        columns=__index_from_values(segments, obs),
        fill_value=0,
    )

    if normalize:
        totals = counts.sum(axis=1)
        empty_groups = list(totals[totals == 0].index)

        if empty_groups:
            raise ValueError(
                f"cannot normalize composition for empty groups: {empty_groups!r}"
            )

        table = counts.div(totals, axis=0)
    else:
        table = counts

    color_values = __normalize_colors(colors, segments, scdata, obs)

    if ax is None:
        fig = plt.figure()
        ax = fig.subplots()
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 3)
        fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else 6)
    else:
        fig = figure_from_axes(ax)

    bar_kwargs = kwargs["bar"] if "bar" in kwargs else {}
    table.plot(
        kind="bar" if orientation == "vertical" else "barh",
        stacked=True,
        width=width,
        color=color_values,
        ax=ax,
        **bar_kwargs,
    )

    ax.set_xlabel(kwargs["xlabel"] if "xlabel" in kwargs else "")
    ax.set_ylabel(kwargs["ylabel"] if "ylabel" in kwargs else "")

    if normalize and percent:
        if orientation == "vertical":
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        else:
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    if "ylim" in kwargs:
        ax.set_ylim(*kwargs["ylim"])
    elif normalize and orientation == "vertical":
        ax.set_ylim(0, 1)

    if "xlim" in kwargs:
        ax.set_xlim(*kwargs["xlim"])
    elif normalize and orientation == "horizontal":
        ax.set_xlim(0, 1)

    if "tick_params" in kwargs:
        ax.tick_params(**kwargs["tick_params"])
    else:
        ax.tick_params(
            axis="both",
            labelsize=kwargs["labelsize"] if "labelsize" in kwargs else 12,
        )

    rotation = kwargs["rotation"] if "rotation" in kwargs else 0
    ticklabels = (
        ax.get_xticklabels() if orientation == "vertical" else ax.get_yticklabels()
    )
    for ticklabel in ticklabels:
        ticklabel.set_rotation(rotation)

    if group_colors is None:
        group_colors = colors_from_uns(scdata, groupby, groups)
    __apply_ticklabel_colors(ax, group_colors, orientation)

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

    existing_legend = ax.get_legend()

    if draw_legend:
        ax.legend(**legend_kwargs)
    elif existing_legend is not None:
        existing_legend.remove()

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None

    return fig, ax
