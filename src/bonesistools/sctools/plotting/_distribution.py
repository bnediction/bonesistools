#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from collections.abc import Sequence as SequenceABC
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
from matplotlib import mlab
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap, is_color_like
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
    _resolve_bool_kwarg,
    _resolve_legend_argument,
    _resolve_toggle_mapping_argument,
    colormap_colors,
    figure_from_axes,
    qualitative_color_values,
    set_axis_label,
    set_window_title,
)

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]
BoxItem = Literal["whiskers", "caps", "boxes", "medians", "fliers", "means"]
BoxPlots = Dict[str, Any]
ViolinPlots = Dict[str, Any]
SingleDistributionReturn = Union[BoxPlots, ViolinPlots]
BoxplotReturn = Union[BoxPlots, Dict[object, BoxPlots]]
ViolinplotReturn = Union[ViolinPlots, Dict[object, ViolinPlots]]
DistributionReturn = Union[
    SingleDistributionReturn,
    Mapping[object, SingleDistributionReturn],
]
ViolinClip = Union[Literal["data"], Tuple[Optional[float], Optional[float]], None]


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


def __apply_violin_colors(
    vp: ViolinPlots,
    color: object,
    alpha: float,
) -> None:

    bodies = cast(Sequence[Any], vp.get("bodies", ()))
    if is_color_like(color):
        for body in bodies:
            body.set_alpha(None)
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha)
        for key in ("cbars", "cmins", "cmaxes"):
            if key in vp:
                vp[key].set_color(color)
        if "cmedians" in vp:
            vp["cmedians"].set_color(black)
            vp["cmedians"].set_linestyle("--")
            vp["cmedians"].set_linewidth(1.0)
        if "cmeans" in vp:
            vp["cmeans"].set_color("C2")
            vp["cmeans"].set_linestyle("--")
            vp["cmeans"].set_linewidth(1.0)
        return None

    colors = list(cast(Sequence[object], color))
    for body, body_color in zip(bodies, colors):
        body.set_alpha(None)
        body.set_facecolor(body_color)
        body.set_edgecolor(body_color)
        body.set_alpha(alpha)
    for key in ("cbars", "cmins", "cmaxes"):
        if key in vp:
            vp[key].set_color(colors)
    if "cmedians" in vp:
        vp["cmedians"].set_color(black)
        vp["cmedians"].set_linestyle("--")
        vp["cmedians"].set_linewidth(1.0)
    if "cmeans" in vp:
        vp["cmeans"].set_color("C2")
        vp["cmeans"].set_linestyle("--")
        vp["cmeans"].set_linewidth(1.0)

    return None


def __apply_summary_artist_kwargs(
    artists: Union[BoxPlots, ViolinPlots],
    *,
    median: Mapping[str, Any],
    mean: Mapping[str, Any],
) -> None:

    summary_kwargs = {
        "medians": median,
        "cmedians": median,
        "means": mean,
        "cmeans": mean,
    }

    for key, kwargs in summary_kwargs.items():
        if not kwargs or key not in artists:
            continue

        artist = artists[key]
        if isinstance(artist, SequenceABC) and not isinstance(artist, str):
            for item in cast(Sequence[Any], artist):
                item.set(**kwargs)
        else:
            cast(Any, artist).set(**kwargs)

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
    points: Optional[Union[bool, Mapping[str, Any]]],
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

    if isinstance(points, MappingABC):
        show_fliers = False if show_fliers is None else show_fliers
        return True, dict(points), show_fliers

    raise TypeError(
        f"unsupported argument type for 'points': "
        f"expected {bool}, {Mapping} or None but received {type(points)}"
    )


def __reject_removed_color_arguments(kwargs: Dict[str, Any]) -> None:

    replacements = {
        "box_colors": "colors=...",
        "point_colors": "points={'colors': ...}",
    }

    for name, replacement in replacements.items():
        if name not in kwargs:
            continue

        raise TypeError(f"unsupported argument {name!r}: use {replacement} instead")


def __reject_removed_distribution_arguments(kwargs: Dict[str, Any]) -> None:

    replacements = {
        "boxplot": "backend-specific keyword arguments",
        "violin": "backend-specific keyword arguments",
        "distribution": "backend-specific keyword arguments",
    }

    for name, replacement in replacements.items():
        if name not in kwargs:
            continue

        raise TypeError(f"unsupported argument {name!r}: use {replacement} instead")


def __reject_removed_legend_arguments(kwargs: Dict[str, Any]) -> None:

    replacements = {
        "legend_fontsize": "legend={'fontsize': ...}",
    }

    for name, replacement in replacements.items():
        if name not in kwargs:
            continue

        raise TypeError(f"unsupported argument {name!r}: use {replacement} instead")


def __reject_violin_boxplot_kwargs(kwargs: Dict[str, Any]) -> None:

    for name in (
        "notch",
        "sym",
        "patch_artist",
        "boxitems_to_color",
        "show_caps",
        "showcaps",
        "show_box",
        "showbox",
        "show_fliers",
        "showfliers",
    ):
        if name in kwargs:
            raise ValueError(
                f"unsupported keyword argument for kind='violin': {name!r}"
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


def __distribution_values(
    series: Union[Series, SeriesGroupBy],
    groups: Optional[Sequence[object]],
) -> Union[Series, Sequence[Series]]:

    if groups is None:
        return cast(Series, series).dropna()

    grouped_series = cast(SeriesGroupBy, series)
    empty = cast(Series, grouped_series.obj).iloc[0:0]
    values = []
    for group in groups:
        try:
            values.append(cast(Series, grouped_series.get_group(group)).dropna())
        except KeyError:
            values.append(empty)
    return values


def __draw_box_distribution(
    ax: Axes,
    values: Union[Series, Sequence[Series]],
    positions: np.ndarray,
    widths: float,
    kwargs: Dict[str, Any],
    labels: Optional[Sequence[object]] = None,
) -> BoxPlots:

    box_plots = cast(
        BoxPlots,
        ax.boxplot(
            x=values,
            positions=positions,
            widths=widths,
            **kwargs,
        ),
    )
    box_plots["positions"] = np.asarray(positions, dtype=float).ravel()
    if labels is not None:
        box_plots["groups"] = list(labels)
    return box_plots


def __draw_violin_distribution(
    ax: Axes,
    values: Union[Series, Sequence[Series]],
    positions: np.ndarray,
    widths: float,
    kwargs: Dict[str, Any],
    cut: float,
    clip: ViolinClip,
    labels: Optional[Sequence[object]] = None,
) -> ViolinPlots:

    violin_kwargs = dict(kwargs)
    violin_widths = violin_kwargs.pop("widths", widths)
    show_means = violin_kwargs.pop("showmeans", False)
    show_medians = violin_kwargs.pop("showmedians", False)
    show_extrema = violin_kwargs.pop("showextrema", True)
    points = int(violin_kwargs.pop("points", 100))
    bw_method = violin_kwargs.pop("bw_method", None)
    quantiles = violin_kwargs.pop("quantiles", None)
    (
        violin_values,
        violin_positions,
        violin_widths,
        violin_quantiles,
        violin_labels,
    ) = __nonempty_violin_inputs(
        values,
        positions=positions,
        widths=violin_widths,
        quantiles=quantiles,
        labels=labels,
    )
    violin_stats = __violin_stats(
        violin_values,
        points=points,
        bw_method=bw_method,
        quantiles=violin_quantiles,
        cut=cut,
        clip=clip,
    )
    violin_plots = cast(
        ViolinPlots,
        ax.violin(
            vpstats=violin_stats,
            positions=violin_positions,
            widths=violin_widths,
            showmeans=show_means,
            showmedians=show_medians,
            showextrema=show_extrema,
            **violin_kwargs,
        ),
    )
    __set_violin_summary_segments(
        violin_plots,
        violin_stats,
        positions=violin_positions,
        widths=violin_widths,
        orientation=cast(Optional[str], violin_kwargs.get("orientation")),
        vert=cast(Optional[bool], violin_kwargs.get("vert")),
    )
    violin_plots["positions"] = violin_positions
    if violin_labels is not None:
        violin_plots["groups"] = violin_labels
    return violin_plots


def __nonempty_violin_inputs(
    values: Union[Series, Sequence[Series]],
    *,
    positions: np.ndarray,
    widths: object,
    quantiles: Optional[Any],
    labels: Optional[Sequence[object]],
) -> Tuple[
    Sequence[Series],
    np.ndarray,
    np.ndarray,
    Sequence[np.ndarray],
    Optional[Sequence[object]],
]:

    value_sets = [values] if isinstance(values, Series) else list(values)
    position_values = np.asarray(positions, dtype=float).ravel()
    width_values = __violin_width_values(widths, n_values=len(position_values))
    quantile_sets = __violin_quantiles(quantiles, n_values=len(value_sets))
    label_values = None if labels is None else list(labels)

    if len(value_sets) != len(position_values):
        raise ValueError(
            "list of violin values and violin positions must have the same length"
        )
    if label_values is not None and len(label_values) != len(value_sets):
        raise ValueError(
            "list of violin labels and violin values must have the same length"
        )

    keep_indices = [
        index
        for index, value_set in enumerate(value_sets)
        if len(cast(Series, value_set).dropna()) > 0
    ]
    if not keep_indices:
        raise ValueError("cannot draw violin distribution from empty values")

    return (
        [cast(Series, value_sets[index]) for index in keep_indices],
        position_values[keep_indices],
        width_values[keep_indices],
        [quantile_sets[index] for index in keep_indices],
        None
        if label_values is None
        else [label_values[index] for index in keep_indices],
    )


def __set_violin_summary_segments(
    violin_plots: ViolinPlots,
    violin_stats: Sequence[Dict[str, Any]],
    *,
    positions: np.ndarray,
    widths: object,
    orientation: Optional[str],
    vert: Optional[bool],
) -> None:

    position_values = np.asarray(positions, dtype=float).ravel()
    width_values = __violin_width_values(widths, n_values=len(position_values))
    horizontal = orientation == "horizontal" or vert is False
    summary_keys = {
        "cmedians": "median",
        "cmeans": "mean",
    }

    for artist_key, stats_key in summary_keys.items():
        if artist_key not in violin_plots:
            continue

        segments = []
        for stats, position, width in zip(violin_stats, position_values, width_values):
            value = float(stats[stats_key])
            lower = position - width / 2
            upper = position + width / 2
            if horizontal:
                segments.append(np.array([[value, lower], [value, upper]]))
            else:
                segments.append(np.array([[lower, value], [upper, value]]))

        violin_plots[artist_key].set_segments(segments)
    return None


def __violin_width_values(
    widths: object,
    n_values: int,
) -> np.ndarray:

    width_values = np.asarray(widths, dtype=float)
    if width_values.ndim == 0:
        return np.full(n_values, float(width_values))
    if len(width_values) != n_values:
        raise ValueError(
            "list of violin widths and violin statistics must have the same length"
        )
    return width_values


def __violin_stats(
    values: Union[Series, Sequence[Series]],
    *,
    points: int,
    bw_method: Any,
    quantiles: Optional[Any],
    cut: float,
    clip: ViolinClip,
) -> Sequence[Dict[str, Any]]:

    value_sets = [values] if isinstance(values, Series) else list(values)
    quantile_sets = __violin_quantiles(quantiles, n_values=len(value_sets))
    stats = []

    for value_set, quantile_set in zip(value_sets, quantile_sets):
        observed = np.asarray(cast(Series, value_set).dropna(), dtype=float)
        if observed.size == 0:
            raise ValueError("cannot draw violin distribution from empty values")

        lower = float(np.nanmin(observed))
        upper = float(np.nanmax(observed))
        clip_lower, clip_upper = __violin_clip_bounds(
            clip,
            lower=lower,
            upper=upper,
        )
        if np.all(observed[0] == observed):
            coords = np.linspace(lower, upper, points)
            vals = (coords == observed[0]).astype(float)
        else:
            kde = mlab.GaussianKDE(observed, bw_method)
            bandwidth = float(np.sqrt(np.squeeze(kde.covariance)))
            domain_lower = lower - cut * bandwidth
            domain_upper = upper + cut * bandwidth
            if clip_lower is not None:
                domain_lower = max(domain_lower, clip_lower)
            if clip_upper is not None:
                domain_upper = min(domain_upper, clip_upper)
            if domain_lower > domain_upper:
                raise ValueError(
                    "invalid argument value for 'clip': "
                    "clipped violin range is empty"
                )
            coords = np.linspace(domain_lower, domain_upper, points)
            vals = kde.evaluate(coords)

        value_stats: Dict[str, Any] = {
            "coords": coords,
            "vals": vals,
            "mean": float(np.mean(observed)),
            "median": float(np.median(observed)),
            "min": lower,
            "max": upper,
            "quantiles": np.percentile(observed, 100 * quantile_set),
        }
        stats.append(value_stats)

    return stats


def __violin_clip_bounds(
    clip: ViolinClip,
    *,
    lower: float,
    upper: float,
) -> Tuple[Optional[float], Optional[float]]:

    if clip == "data":
        return lower, upper
    if clip is None:
        return None, None
    if isinstance(clip, str):
        raise ValueError(
            "invalid argument value for 'clip': "
            "expected 'data', None, or a 2-length tuple such as (0, None)"
        )
    if not isinstance(clip, tuple) or len(clip) != 2:
        raise TypeError(
            "unsupported argument type for 'clip': "
            "expected 'data', None, or a 2-length tuple such as (0, None)"
        )

    clip_lower = None if clip[0] is None else float(clip[0])
    clip_upper = None if clip[1] is None else float(clip[1])
    if clip_lower is not None and clip_upper is not None and clip_lower > clip_upper:
        raise ValueError(
            "invalid argument value for 'clip': "
            f"expected lower <= upper but received {clip}"
        )

    return clip_lower, clip_upper


def __violin_quantiles(
    quantiles: Optional[Any],
    n_values: int,
) -> Sequence[np.ndarray]:

    if quantiles is None:
        return [np.array([])] * n_values

    quantile_array = np.asarray(quantiles, dtype=float)
    if quantile_array.ndim == 0:
        quantile_array = quantile_array.reshape(1)
    if quantile_array.ndim == 1:
        return [quantile_array] * n_values
    if len(quantile_array) != n_values:
        raise ValueError(
            "list of violinplot statistics and quantile values must have "
            "the same length"
        )
    return [np.asarray(q, dtype=float) for q in quantile_array]


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: str,
    hue: str,
    kind: Literal["box"] = "box",
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    kind: Literal["box"] = "box",
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes, BoxPlots]: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: str,
    hue: str,
    kind: Literal["violin"],
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes, Dict[object, ViolinPlots]]: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: None = None,
    kind: Literal["violin"],
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes, ViolinPlots]: ...


@overload
def distribution(
    scdata: ScData,
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: None = None,
    kind: Literal["box", "violin"] = "box",
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    kind: Literal["box", "violin"] = "box",
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes, DistributionReturn]]: ...


@anndata_or_mudata_checker
def distribution(
    scdata: ScData,  # type: ignore
    obs: str,
    *,
    groupby: Optional[str] = None,
    hue: Optional[str] = None,
    kind: Literal["box", "violin"] = "box",
    sort: Literal["ascending", "descending", "preserve"] = "preserve",
    cut: float = 2,
    clip: ViolinClip = "data",
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    widths: float = 0.5,
    groupby_spacing: float = 0.3,
    hue_spacing: float = 0.1,
    colors: Optional[Colors] = None,
    alpha: float = 1.0,
    points: Optional[Union[bool, Mapping[str, Any]]] = None,
    median: Union[bool, Mapping[str, Any]] = True,
    mean: Union[bool, Mapping[str, Any]] = False,
    figwidth: Optional[float] = None,
    figheight: Optional[float] = None,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes, DistributionReturn]]:
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
    kind: {"box", "violin"} (default: "box")
        Distribution representation.
    sort: {"ascending", "descending", "preserve"} (default: "preserve")
        Sort groups by median value or preserve category order.
    cut: float (default: 2)
        Tail extension for `kind="violin"` only, expressed in KDE bandwidth
        units beyond the observed minimum and maximum values. Use `cut=0` to
        restrict violin bodies to the observed data range. Ignored for
        `kind="box"`.
    clip: {"data"}, tuple or None (default: "data")
        Bounds for the violin KDE evaluation range. With `"data"`, use the
        observed minimum and maximum values of each violin. Use a tuple such as
        `clip=(0, None)` for positive observation-level metrics. Use None to
        leave the KDE range unclipped. Ignored for `kind="box"`.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    legend: bool or mapping (default: True)
        Legend configuration. False disables the legend. True draws the legend
        using default Matplotlib parameters. If a mapping is provided, it is
        forwarded as keyword arguments to `Axes.legend`.
    widths: float (default: 0.5)
        Distribution width.
    groupby_spacing: float (default: 0.3)
        Spacing between groups.
    hue_spacing: float (default: 0.1)
        Spacing between hues within each group.
    colors: sequence, iterator, colormap or mapping, optional
        Distribution colors. When `hue` is specified, colors are mapped to hue
        values.
    alpha: float (default: 1.0)
        Violin body opacity. When points are shown, also used as the default
        point opacity. Explicit `points={"alpha": ...}` values take precedence.
    points: bool or mapping, optional
        If None, draw points only when boxplot fliers are hidden. If True, draw
        points with default scatter parameters. If False, do not draw points. If
        a mapping is provided, it is forwarded as keyword arguments to
        `Axes.scatter`. Use `points={"colors": ...}` to specify point colors.
    median: bool or mapping (default: True)
        Median line configuration. False disables the median line. True draws
        it with default parameters. For `kind="box"`, a mapping is forwarded to
        `medianprops`. For `kind="violin"`, it is applied to the `cmedians`
        artist.
    mean: bool or mapping (default: False)
        Mean line configuration. False disables the mean line. True draws it
        with default parameters. For `kind="box"`, a mapping is forwarded to
        `meanprops`. For `kind="violin"`, it is applied to the `cmeans` artist.
    figwidth: float, optional
        Figure width.
    figheight: float, optional
        Figure height.
    xlabel: str or mapping, optional
        X-axis label. If a mapping is provided, it is forwarded as keyword
        arguments to `Axes.set_xlabel` and must contain a `label` key.
    ylabel: str or mapping, optional
        Y-axis label. If a mapping is provided, it is forwarded as keyword
        arguments to `Axes.set_ylabel` and must contain a `label` key.
    ax: Axes, optional
        Existing axes used for drawing.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Additional keyword arguments forwarded to the underlying Matplotlib
        drawing function (`Axes.boxplot` or `Axes.violinplot`, depending on
        `kind`). For example, `flierprops={...}` is forwarded to `Axes.boxplot`
        to style boxplot fliers.

    Returns
    -------
    tuple
        Matplotlib figure, axes and distribution artists.

    Raises
    ------
    ValueError
        If `groupby` and `hue` are inconsistently specified, if `kind` is not
        `"box"` or `"violin"`, or if `sort` is not `"ascending"`,
        `"descending"` or `"preserve"`.
    """

    __reject_removed_color_arguments(kwargs)
    __reject_removed_distribution_arguments(kwargs)
    __reject_removed_legend_arguments(kwargs)
    if cut < 0:
        raise ValueError(
            f"invalid argument value for 'cut': "
            f"expected a non-negative number but received {cut}"
        )
    if not 0 <= alpha <= 1:
        raise ValueError(
            f"invalid argument value for 'alpha': "
            f"expected a value between 0 and 1 but received {alpha}"
        )
    kind = cast(
        Literal["box", "violin"],
        _as_literal(kind, choices=("box", "violin"), name="kind"),
    )
    if kind == "violin":
        __reject_violin_boxplot_kwargs(kwargs)
    legend = _resolve_legend_argument(legend, kwargs)

    draw_median, median_kwargs = _resolve_toggle_mapping_argument(
        median,
        kwargs,
        name="median",
        deprecated_names=("show_median", "show_medians", "showmedians"),
        default=True,
        stacklevel=3,
    )
    draw_mean, mean_kwargs = _resolve_toggle_mapping_argument(
        mean,
        kwargs,
        name="mean",
        deprecated_names=("show_mean", "show_means", "showmeans"),
        default=False,
        stacklevel=3,
    )
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
    boxitems_to_color = cast(
        Optional[Sequence[BoxItem]],
        kwargs.pop("boxitems_to_color", ("whiskers", "caps", "boxes")),
    )
    show_points, point_kwargs, show_fliers = __resolve_points_argument(
        points,
        kwargs,
        show_fliers,
    )
    point_colors = cast(Optional[Colors], point_kwargs.pop("colors", None))
    if show_points:
        point_kwargs.setdefault("alpha", alpha)
    distribution_kwargs = dict(kwargs)
    kwargs.clear()

    if kind == "box":
        distribution_kwargs.setdefault("showmeans", draw_mean)
        if median_kwargs:
            distribution_kwargs.setdefault("medianprops", {})
            distribution_kwargs["medianprops"].update(median_kwargs)
        if mean_kwargs:
            distribution_kwargs.setdefault("meanprops", {})
            distribution_kwargs["meanprops"].update(mean_kwargs)
        if show_caps_specified:
            distribution_kwargs["showcaps"] = show_caps
        else:
            distribution_kwargs.setdefault("showcaps", show_caps)
        if show_box_specified:
            distribution_kwargs["showbox"] = show_box
        else:
            distribution_kwargs.setdefault("showbox", show_box)
        if show_fliers_specified:
            distribution_kwargs["showfliers"] = show_fliers
        else:
            distribution_kwargs.setdefault("showfliers", show_fliers)
    else:
        distribution_kwargs.setdefault("widths", widths)
        distribution_kwargs["showmeans"] = draw_mean
        distribution_kwargs["showmedians"] = draw_median
        distribution_kwargs.setdefault("showextrema", False)

    if ax is None:
        fig, ax = plt.subplots()
        fig.set_figheight(6 if figheight is None else figheight)
        fig.set_figwidth(
            (5 if groupby is None else 8) if figwidth is None else figwidth
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
        distribution_values = __distribution_values(
            cast(Union[Series, SeriesGroupBy], series),
            None if groups is None else cast(Sequence[object], groups),
        )
        if kind == "box":
            distribution_artists = __draw_box_distribution(
                ax=ax,
                values=distribution_values,
                positions=positions,
                widths=widths,
                kwargs=distribution_kwargs,
                labels=None if groups is None else cast(Sequence[object], groups),
            )
            if draw_median is False:
                for median_artist in cast(
                    Sequence[Any],
                    distribution_artists["medians"],
                ):
                    median_artist.set(linewidth=0)
            __apply_summary_artist_kwargs(
                distribution_artists,
                median=median_kwargs,
                mean=mean_kwargs,
            )
            if colors is not None and boxitems_to_color:
                __apply_box_colors(
                    bp=distribution_artists,
                    color=colors,
                    items=boxitems_to_color,
                )
        else:
            distribution_artists = __draw_violin_distribution(
                ax=ax,
                values=distribution_values,
                positions=positions,
                widths=widths,
                kwargs=distribution_kwargs,
                cut=cut,
                clip=clip,
                labels=None if groups is None else cast(Sequence[object], groups),
            )
            violin_colors = colors
            if violin_colors is None:
                violin_colors = cast(
                    Colors,
                    qualitative_color_values(
                        1 if groups is None else len(cast(Sequence[object], groups)),
                        QUALITATIVE_COLORS,
                        generate_colormap,
                    ),
            )
            __apply_violin_colors(distribution_artists, violin_colors, alpha=alpha)
            __apply_summary_artist_kwargs(
                distribution_artists,
                median=median_kwargs,
                mean=mean_kwargs,
            )
            if groups is None:
                ax.set_xlim(-widths * 0.65, widths * 0.65)
    else:
        grouped_artists: Dict[object, SingleDistributionReturn] = {}
        hue_values = cast(Sequence[object], hues)

        artist_colors = __resolve_hue_colors(
            colors,
            hue_values,
            (
                [black] * len(hue_values)
                if kind == "box" and show_points is True
                else cast(
                    Colors,
                    qualitative_color_values(
                        len(hue_values),
                        QUALITATIVE_COLORS,
                        generate_colormap,
                    ),
                )
            ),
        )
        point_colors = __resolve_hue_colors(
            point_colors,
            hue_values,
            (
                (
                    [white] * len(hue_values)
                    if kind == "violin"
                    else cast(
                        Colors,
                        qualitative_color_values(
                            len(hue_values),
                            QUALITATIVE_COLORS,
                            generate_colormap,
                        ),
                    )
                )
                if show_points is True
                else [white] * len(hue_values)
            ),
        )

        if kind == "box":
            if "medianprops" not in distribution_kwargs:
                distribution_kwargs["medianprops"] = {}
            if (
                "color" not in distribution_kwargs["medianprops"]
                and boxitems_to_color is not None
                and "medians" not in boxitems_to_color
            ):
                distribution_kwargs["medianprops"]["color"] = black

        positions_iterator = iter(positions)
        grouped_series_by_hue = cast(Mapping[object, SeriesGroupBy], series)

        for h, values in grouped_series_by_hue.items():
            distribution_values = __distribution_values(
                values,
                cast(Sequence[object], groups),
            )
            if kind == "box":
                grouped_artists[h] = __draw_box_distribution(
                    ax=ax,
                    values=distribution_values,
                    positions=next(positions_iterator),
                    widths=widths,
                    kwargs=distribution_kwargs,
                    labels=cast(Sequence[object], groups),
                )
            else:
                grouped_artists[h] = __draw_violin_distribution(
                    ax=ax,
                    values=distribution_values,
                    positions=next(positions_iterator),
                    widths=widths,
                    kwargs=distribution_kwargs,
                    cut=cut,
                    clip=clip,
                    labels=cast(Sequence[object], groups),
                )

        if kind == "box":
            if draw_median is False:
                for bp in grouped_artists.values():
                    for median_artist in cast(
                        Sequence[Any],
                        cast(BoxPlots, bp)["medians"],
                    ):
                        median_artist.set(linewidth=0)
            for bp in grouped_artists.values():
                __apply_summary_artist_kwargs(
                    cast(BoxPlots, bp),
                    median=median_kwargs,
                    mean=mean_kwargs,
                )

            if boxitems_to_color:
                for h in hue_values:
                    __apply_box_colors(
                        bp=cast(BoxPlots, grouped_artists[h]),
                        color=cast(Mapping[object, object], artist_colors)[h],
                        items=boxitems_to_color,
                    )
        else:
            for h in hue_values:
                __apply_violin_colors(
                    cast(ViolinPlots, grouped_artists[h]),
                    cast(Mapping[object, object], artist_colors)[h],
                    alpha=alpha,
                )
                __apply_summary_artist_kwargs(
                    cast(ViolinPlots, grouped_artists[h]),
                    median=median_kwargs,
                    mean=mean_kwargs,
                )

        if legend is not False:
            handles = []
            for h in hue_values:
                if kind == "violin":
                    bodies = cast(
                        Sequence[Any],
                        cast(ViolinPlots, grouped_artists[h])["bodies"],
                    )
                    legend_edgecolor = bodies[0].get_edgecolor()[0]
                    legend_facecolor = bodies[0].get_facecolor()[0]
                else:
                    legend_edgecolor = cast(Mapping[object, object], artist_colors)[
                        h
                    ]
                    legend_facecolor = cast(Mapping[object, object], point_colors)[h]
                handles.append(
                    mpatches.Patch(
                        edgecolor=cast(Any, legend_edgecolor),
                        facecolor=cast(Any, legend_facecolor),
                        label=h,
                    )
                )
            if legend is True:
                ax.legend(handles=handles)
            else:
                ax.legend(handles=handles, **legend)

        distribution_artists = grouped_artists

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

    if xlabel is not None:
        set_axis_label(ax, "xlabel", xlabel)
    if ylabel is not None:
        set_axis_label(ax, "ylabel", ylabel)

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax, distribution_artists


def boxplot(
    *args: Any,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes, DistributionReturn]]:
    """
    Deprecated alias for `distribution`.
    """

    _warn_deprecated(
        "`bt.sct.pl.boxplot`",
        replacement="`bt.sct.pl.distribution`",
        stacklevel=2,
    )
    return distribution(*args, **kwargs)
