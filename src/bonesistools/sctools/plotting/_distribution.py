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

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure
from pandas import Series
from pandas.core.groupby.generic import SeriesGroupBy

from ..._compat import Literal
from ..._validation import _as_literal
from ..._warnings import _warn_deprecated, _warn_deprecated_argument
from .._typing import ScData, anndata_or_mudata_checker
from ._colors import (
    QUALITATIVE_COLORS,
    black,
    generate_colormap,
    gray,
    white,
)
from ._utils import (
    LegendArgument,
    _resolve_bool_kwarg,
    _resolve_legend_argument,
    colormap_colors,
    figure_from_axes,
    qualitative_color_values,
    set_window_title,
)

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]
BoxItem = Literal["whiskers", "caps", "boxes", "medians", "fliers", "means"]
BoxPlots = Mapping[str, Any]
BoxplotReturn = Union[BoxPlots, Dict[object, BoxPlots]]
Legend = LegendArgument
Points = Optional[Union[bool, Dict[str, Any]]]


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
    ax: Axes,
    groups: Optional[Sequence[object]] = None,
    hues: Optional[Sequence[object]] = None,
    **kwargs: Any,
) -> None:

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


def __resolve_points_argument(
    points: Points,
    kwargs: Dict[str, Any],
    show_fliers: Optional[bool],
) -> Tuple[bool, Dict[str, Any], bool]:

    for deprecated_name in ("showpoints", "show_points"):
        if deprecated_name not in kwargs:
            continue

        _warn_deprecated_argument(deprecated_name, "points", stacklevel=3)
        if points is not None:
            raise TypeError(
                f"invalid argument combination: use either '{deprecated_name}' "
                "or 'points', not both"
            )
        points = cast(bool, kwargs.pop(deprecated_name))

    if "scatter" in kwargs:
        _warn_deprecated_argument("scatter", "points", stacklevel=3)
        if points is not None:
            raise TypeError(
                "invalid argument combination: use either deprecated "
                "'scatter' or 'points', not both"
            )
        points = cast(Dict[str, Any], kwargs.pop("scatter"))

    if points is None:
        show_fliers = True if show_fliers is None else show_fliers
        return not show_fliers, {}, show_fliers

    if isinstance(points, bool):
        show_fliers = not points if show_fliers is None else show_fliers
        return points, {}, show_fliers

    if isinstance(points, dict):
        show_fliers = False if show_fliers is None else show_fliers
        return True, dict(points), show_fliers

    raise TypeError(
        f"unsupported argument type for 'points': "
        f"expected {bool}, {dict} or None but received {type(points)}"
    )


def __resolve_boxplot_argument(
    boxplot: Optional[Dict[str, Any]],
    kwargs: Dict[str, Any],
) -> Dict[str, Any]:

    if boxplot is None:
        boxplot_kwargs: Dict[str, Any] = {}
    elif isinstance(boxplot, dict):
        boxplot_kwargs = dict(boxplot)
    else:
        raise TypeError(
            f"unsupported argument type for 'boxplot': "
            f"expected {dict} or None but received {type(boxplot)}"
        )

    for name in ("notch", "sym", "patch_artist", "vert"):
        if name in kwargs:
            boxplot_kwargs[name] = kwargs.pop(name)

    return boxplot_kwargs


def __reject_removed_color_arguments(kwargs: Dict[str, Any]) -> None:

    replacements = {
        "box_colors": "boxplot={'colors': ...}",
        "point_colors": "points={'colors': ...}",
    }

    for name, replacement in replacements.items():
        if name not in kwargs:
            continue

        raise TypeError(
            f"unsupported argument {name!r}: use {replacement} instead"
        )


def __resolve_hue_colors(
    colors: Optional[Colors],
    values: Sequence[object],
    default: Colors,
) -> Mapping[object, object]:

    if colors is None:
        colors = default

    if isinstance(colors, ListedColormap):
        colors = colormap_colors(colors)
    elif isinstance(colors, Colormap):
        colors = cast(Sequence[object], colors(np.linspace(0, 1, len(values))))

    if isinstance(colors, MappingABC):
        return cast(Mapping[object, object], colors)

    color_values = cast(Sequence[object], colors)
    return {value: color_values[index] for index, value in enumerate(values)}


def __resolve_point_kwargs(
    point_kwargs: Dict[str, Any],
    point_colors: Colors,
    hues: Optional[Sequence[object]],
) -> Dict[Any, Any]:

    if hues is None:
        resolved = dict(point_kwargs)
        resolved.setdefault("s", 1)
        resolved.setdefault("alpha", 0.7)
        resolved.setdefault("facecolors", point_colors)
        resolved.setdefault("edgecolors", "none")
        return resolved

    hue_values = cast(Sequence[object], hues)
    shared_kwargs = {
        key: value for key, value in point_kwargs.items() if key not in hue_values
    }
    resolved_by_hue: Dict[Any, Dict[str, Any]] = {}

    for hue in hue_values:
        resolved = dict(shared_kwargs)
        hue_kwargs = cast(Mapping[object, Any], point_kwargs).get(hue, {})
        if isinstance(hue_kwargs, MappingABC):
            resolved.update(cast(Mapping[str, Any], hue_kwargs))
        resolved.setdefault("s", 1)
        resolved.setdefault("alpha", 0.7)
        resolved.setdefault(
            "facecolors",
            cast(Mapping[object, object], point_colors)[hue],
        )
        resolved.setdefault("edgecolors", "none")
        resolved_by_hue[hue] = resolved

    return resolved_by_hue


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: str,
    hue: str,
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Legend = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    points: Points = None,
    boxplot: Optional[Dict[str, Any]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes, Dict[object, BoxPlots]]: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: None = None,
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Legend = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    points: Points = None,
    boxplot: Optional[Dict[str, Any]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes, BoxPlots]: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: None = None,
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Legend = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    points: Points = None,
    boxplot: Optional[Dict[str, Any]] = None,
    ax: Optional[Axes] = None,
    outfile: Path,
    **kwargs: Any,
) -> None: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: Optional[str] = None,
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Legend = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    points: Points = None,
    boxplot: Optional[Dict[str, Any]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes, BoxplotReturn]]: ...


@anndata_or_mudata_checker
def distribution(
    scdata: ScData,  # type: ignore
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: Optional[str] = None,
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Legend = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    points: Points = None,
    boxplot: Optional[Dict[str, Any]] = None,
    ax: Optional[Axes] = None,
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
    sort: {"ascending", "descending", "preserve"} (default: "preserve")
        Sort groups by median value or preserve category order.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    legend: bool or dict (default: True)
        If True, draw the default hue legend. If False, do not draw a legend.
        If a dictionary, draw the legend and pass these keyword arguments to
        `Axes.legend`.
    widths: float (default: 0.5)
        Box width.
    groupby_spacing: float (default: 0.3)
        Spacing between groups.
    hue_spacing: float (default: 0.1)
        Spacing between hues within each group.
    points: bool or dict, optional
        If None, draw points only when boxplot fliers are hidden. If True, draw
        points with default scatter parameters. If False, do not draw points. If
        a dictionary, draw points and pass these keyword arguments to
        `Axes.scatter`. Use `points={"colors": ...}` to specify point colors.
    boxplot: dict, optional
        Keyword arguments passed to `Axes.boxplot`. Use
        `boxplot={"colors": ...}` to specify boxplot colors.
    ax: Axes, optional
        Existing axes used for drawing.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - show_medians[bool]: draw median lines in boxplots
        - show_means[bool]: draw mean markers in boxplots
        - show_caps[bool]: draw caps at the end of whiskers
        - show_box[bool]: draw box bodies
        - show_fliers[bool]: draw flier points
        - notch[bool]: forwarded to `Axes.boxplot`
        - sym[str]: forwarded to `Axes.boxplot`
        - patch_artist[bool]: forwarded to `Axes.boxplot`
        - vert[bool]: forwarded to `Axes.boxplot`
        - boxitems_to_color[sequence]: boxplot artist groups colored by
          `boxplot["colors"]`

    Returns
    -------
    tuple
        Matplotlib figure, axes and boxplot artists.

    Raises
    ------
    ValueError
        If `groupby` and `hue` are inconsistently specified, or if `sort` is
        not `"ascending"`, `"descending"` or `"preserve"`.
    """

    __reject_removed_color_arguments(kwargs)
    draw_legend, legend_kwargs = _resolve_legend_argument(legend, kwargs)
    boxplot_kwargs = __resolve_boxplot_argument(boxplot, kwargs)
    box_colors = cast(Optional[Colors], boxplot_kwargs.pop("colors", None))
    show_medians = bool(
        _resolve_bool_kwarg(kwargs, "show_medians", "showmedians", True)
    )
    show_means_specified = "show_means" in kwargs or "showmeans" in kwargs
    show_means = bool(_resolve_bool_kwarg(kwargs, "show_means", "showmeans", False))
    show_caps_specified = "show_caps" in kwargs or "showcaps" in kwargs
    show_caps = bool(_resolve_bool_kwarg(kwargs, "show_caps", "showcaps", True))
    show_box_specified = "show_box" in kwargs or "showbox" in kwargs
    show_box = bool(_resolve_bool_kwarg(kwargs, "show_box", "showbox", True))
    show_fliers_specified = "show_fliers" in kwargs or "showfliers" in kwargs
    show_fliers = _resolve_bool_kwarg(
        kwargs,
        "show_fliers",
        "showfliers",
        None,
    )
    if show_fliers is None and "showfliers" in boxplot_kwargs:
        show_fliers = cast(bool, boxplot_kwargs["showfliers"])

    show_points, point_kwargs, show_fliers = __resolve_points_argument(
        points,
        kwargs,
        show_fliers,
    )
    point_colors = cast(Optional[Colors], point_kwargs.pop("colors", None))

    if show_means_specified:
        boxplot_kwargs["showmeans"] = show_means
    else:
        boxplot_kwargs.setdefault("showmeans", show_means)
    if show_caps_specified:
        boxplot_kwargs["showcaps"] = show_caps
    else:
        boxplot_kwargs.setdefault("showcaps", show_caps)
    if show_box_specified:
        boxplot_kwargs["showbox"] = show_box
    else:
        boxplot_kwargs.setdefault("showbox", show_box)
    if show_fliers_specified:
        boxplot_kwargs["showfliers"] = show_fliers
    else:
        boxplot_kwargs.setdefault("showfliers", show_fliers)

    boxitems_to_color = cast(
        Optional[Sequence[BoxItem]],
        kwargs.pop("boxitems_to_color", ("whiskers", "caps", "boxes")),
    )

    if ax is None:
        fig, ax = plt.subplots()
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 6)
        fig.set_figwidth(
            kwargs["figwidth"]
            if "figwidth" in kwargs
            else 5 if groupby is None else 8
        )
    else:
        fig = figure_from_axes(ax)

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
        sort = _as_literal(
            sort,
            choices=("ascending", "descending", "preserve"),
            name="sort",
        )
        if sort == "preserve":
            groups = list(scdata.obs[groupby].astype("category").cat.categories)
            series = cast(
                SeriesGroupBy, scdata.obs.groupby(by=groupby, observed=True)[obs]
            )
        else:
            series = cast(
                SeriesGroupBy, scdata.obs.groupby(by=groupby, observed=True)[obs]
            )
            groups = list(
                series.median().sort_values(ascending=(sort == "ascending")).index
            )
        if hue is None:
            hues = None
        else:
            hues = list(scdata.obs[hue].astype("category").cat.categories)
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
            ax.boxplot(
                x=(
                    cast(Series, series).dropna()
                    if groupby is None
                    else [
                        cast(SeriesGroupBy, series).get_group(g).dropna()
                        for g in cast(Sequence[object], groups)
                    ]
                ),
                positions=positions,
                widths=widths,
                **boxplot_kwargs,
            ),
        )
        if show_medians is False:
            for median in cast(Sequence[Any], bps["medians"]):
                median.set(linewidth=0)
        if box_colors is not None and boxitems_to_color:
            __apply_box_colors(
                bp=bps,
                color=box_colors,
                items=boxitems_to_color,
            )
    else:
        grouped_bps: Dict[object, BoxPlots] = {}
        hue_values = cast(Sequence[object], hues)

        box_colors = __resolve_hue_colors(
            box_colors,
            hue_values,
            [black] * len(hue_values)
            if show_points is True
            else cast(
                Colors,
                qualitative_color_values(
                    len(hue_values),
                    QUALITATIVE_COLORS,
                    generate_colormap,
                ),
            ),
        )
        point_colors = __resolve_hue_colors(
            point_colors,
            hue_values,
            cast(
                Colors,
                qualitative_color_values(
                    len(hue_values),
                    QUALITATIVE_COLORS,
                    generate_colormap,
                ),
            )
            if show_points is True
            else [white] * len(hue_values),
        )

        if "medianprops" not in boxplot_kwargs:
            boxplot_kwargs["medianprops"] = {}
        if (
            "color" not in boxplot_kwargs["medianprops"]
            and boxitems_to_color is not None
            and "medians" not in boxitems_to_color
        ):
            boxplot_kwargs["medianprops"]["color"] = black

        positions_iterator = iter(positions)
        grouped_series_by_hue = cast(Mapping[object, SeriesGroupBy], series)

        for h, values in grouped_series_by_hue.items():
            grouped_bps[h] = cast(
                BoxPlots,
                ax.boxplot(
                    x=[
                        values.get_group(g).dropna()
                        for g in cast(Sequence[object], groups)
                    ],
                    positions=next(positions_iterator),
                    widths=widths,
                    **boxplot_kwargs,
                ),
            )

        if show_medians is False:
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

        if draw_legend is True:
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
            ax.legend(handles=handles, **legend_kwargs)

        bps = grouped_bps

    if show_points:
        resolved_point_kwargs = __resolve_point_kwargs(
            point_kwargs,
            cast(Colors, point_colors),
            hues,
        )
        __add_points(
            series=series,
            positions=positions,
            scale=widths / 8,
            ax=ax,
            groups=groups,
            hues=hues,
            **resolved_point_kwargs,
        )

    ylim_min, ylim_max = ax.get_ylim()
    ylim_diff = ylim_max - ylim_min
    ax.set_ylim(ylim_min - ylim_diff / 50, ylim_max + ylim_diff / 50)

    if groupby is not None:
        xticks = positions.mean(axis=0) if positions.ndim == 2 else positions
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(group) for group in cast(Sequence[object], groups)])
    else:
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax, bps


def boxplot(*args: Any, **kwargs: Any) -> Optional[Tuple[Figure, Axes, BoxplotReturn]]:
    """
    Deprecated alias for `distribution`.
    """

    _warn_deprecated(
        "`bt.sct.pl.boxplot`",
        replacement="`bt.sct.pl.distribution`",
        stacklevel=2,
    )
    return distribution(*args, **kwargs)
