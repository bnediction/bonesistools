#!/usr/bin/env python

from __future__ import annotations

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
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from ..._warnings import _warn_deprecated, _warn_deprecated_argument
from .._typing import anndata_checker
from ..tools._utils import _UNSET
from ._colors import (
    CLASSIC_COLORS,
    QUALITATIVE_COLORS,
    blue,
    generate_colormap,
    gray,
)
from ._utils import (
    _resolve_legend_argument,
    colormap_colors,
    figure_from_axes,
    normalize_color,
    qualitative_color_values,
    set_axis_label,
    set_window_title,
)

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]


def _resolve_feature_argument(feature: Any, gene: Any, *, stacklevel: int) -> str:

    if gene is not _UNSET:
        _warn_deprecated_argument("gene", "feature", stacklevel=stacklevel)
        if feature is not _UNSET:
            raise TypeError(
                "received both 'feature' and deprecated 'gene'; "
                "please use only 'feature'"
            )
        feature = gene

    if feature is _UNSET:
        raise TypeError("missing required argument: 'feature'")

    return cast(str, feature)


def _resolve_expression_argument(
    expression: Optional[str],
    layer: Any,
    *,
    stacklevel: int,
) -> Optional[str]:

    if layer is not _UNSET:
        _warn_deprecated_argument("layer", "expression", stacklevel=stacklevel)
        if expression is not None:
            raise TypeError(
                "received both 'expression' and deprecated 'layer'; "
                "please use only 'expression'"
            )
        expression = cast(Optional[str], layer)

    return expression


def _resolve_clip_outliers_argument(clip_outliers: bool, clip: Any) -> bool:

    if clip is _UNSET:
        return clip_outliers

    _warn_deprecated_argument("clip", "clip_outliers", stacklevel=3)
    if clip_outliers is not False:
        raise TypeError(
            "received both 'clip_outliers' and deprecated 'clip'; "
            "please use only 'clip_outliers'"
        )

    return cast(bool, clip)


def _resolve_show_global_argument(show_global: bool, not_all: Any) -> bool:

    if not_all is _UNSET:
        return show_global

    _warn_deprecated_argument("not_all", "show_global", stacklevel=3)
    if show_global is not True:
        raise TypeError(
            "received both 'show_global' and deprecated 'not_all'; "
            "please use only 'show_global'"
        )

    return not cast(bool, not_all)


def _counts_vector(
    adata: AnnData,
    feature: str,
    expression: Optional[str],
) -> np.ndarray:

    counts = adata[:, feature].layers[expression] if expression else adata[:, feature].X

    if scipy.sparse.issparse(counts):
        counts = cast(Any, counts).toarray()

    return np.asarray(counts).squeeze()


@overload
def density(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    show_global: bool = True,
    clip_outliers: bool = False,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    clip: Any = _UNSET,
    not_all: Any = _UNSET,
    **kwargs: Any,
) -> Tuple[Figure, Axes]: ...


@overload
def density(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    show_global: bool = True,
    clip_outliers: bool = False,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Path,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    clip: Any = _UNSET,
    not_all: Any = _UNSET,
    **kwargs: Any,
) -> None: ...


@overload
def density(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    show_global: bool = True,
    clip_outliers: bool = False,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    clip: Any = _UNSET,
    not_all: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]: ...


@anndata_checker
def density(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    show_global: bool = True,
    clip_outliers: bool = False,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    clip: Any = _UNSET,
    not_all: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw feature-related density function using kernel density estimation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    feature: str
        Feature of interest to plot.
    expression: str, optional
        Layer to use instead of `adata.X`.
    obs: str, optional
        Observation column in `adata.obs` defining groups.
    colors: Colors (optional, default: None)
        Colors used for density curves.
    show_global: bool (default: True)
        If True, draw the density function using all observations.
    clip_outliers: bool (default: False)
        If True, clip density between the minimum value and the quantile at 99%.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    legend: bool or mapping (default: True)
        Legend configuration. False disables the legend. True draws the legend
        using default Matplotlib parameters. If a mapping is provided, it is
        forwarded as keyword arguments to `Axes.legend`.
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
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major
          formatter on x- and y-axis

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.

    Raises
    ------
    ValueError
        If `show_global=False` while `obs` is None.
    KeyError
        If `obs` is specified but not found in `adata.obs`.
    """

    feature = _resolve_feature_argument(feature, gene, stacklevel=2)
    expression = _resolve_expression_argument(expression, layer, stacklevel=2)
    clip_outliers = _resolve_clip_outliers_argument(clip_outliers, clip)
    show_global = _resolve_show_global_argument(show_global, not_all)

    if obs is None and show_global is False:
        raise ValueError(
            "invalid argument values for 'obs' and 'show_global': "
            "expected show_global=True when obs is None"
        )
    legend = _resolve_legend_argument(legend, kwargs)

    import seaborn as sns

    counts = pd.DataFrame(
        {"counting": _counts_vector(adata, feature, expression)},
        index=adata.obs.index,
    )

    if obs:
        counts[obs] = adata.obs[obs]
        if not colors:
            cluster_number = len(adata.obs[obs].astype("category").cat.categories)
            colors = qualitative_color_values(
                cluster_number,
                QUALITATIVE_COLORS,
                generate_colormap,
            )
        elif isinstance(colors, Mapping):
            colors = [
                colors[cluster]
                for cluster in adata.obs[obs].astype("category").cat.categories
            ]
        if isinstance(colors, Colormap):
            colors = colormap_colors(colors)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = figure_from_axes(ax)

    q = np.quantile(counts["counting"], 0.99)
    clip_range = cast(
        Optional[Tuple[float, float]],
        (float(min(counts["counting"])), float(q)) if clip_outliers is True else None,
    )

    if show_global is True:
        sns.kdeplot(
            data=cast(Any, counts["counting"]),
            ax=ax,
            color=cast(Any, gray),
            fill=True,
            clip=clip_range,
            label="all",
        )

    if obs is not None:
        color_values = cast(Sequence[object], colors)

        for _cluster, _color in zip(
            adata.obs[obs].astype("category").cat.categories, color_values
        ):
            _color = normalize_color(_color)
            sns.kdeplot(
                data=cast(Any, counts.loc[counts[obs] == _cluster]["counting"]),
                ax=ax,
                color=cast(Any, _color),
                fill=False,
                clip=clip_range,
                label=_cluster,
            )

    if xlabel is not None:
        set_axis_label(ax, "xlabel", xlabel)
    if ylabel is not None:
        set_axis_label(ax, "ylabel", ylabel)

    if title:
        if isinstance(title, str):
            set_window_title(fig, title)
            ax.set_title(title)
        elif isinstance(title, dict):
            set_window_title(fig, title["label"])
            ax.set_title(**title)

    if obs:
        if legend is True:
            ax.legend()
        elif legend is not False:
            ax.legend(**legend)

    (
        ax.xaxis.set_major_formatter(kwargs["formatter"])
        if "formatter" in kwargs
        else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    )
    (
        ax.yaxis.set_major_formatter(kwargs["formatter"])
        if "formatter" in kwargs
        else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    )

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax


@overload
def cdf(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: None = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    **kwargs: Any,
) -> Tuple[Figure, Axes]: ...


@overload
def cdf(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Path,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    **kwargs: Any,
) -> None: ...


@overload
def cdf(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]: ...


def cdf(
    adata: AnnData,
    feature: Any = _UNSET,
    *,
    expression: Optional[str] = None,
    obs: Optional[str] = None,
    colors: Optional[Colors] = None,
    legend: Union[bool, Mapping[str, Any]] = True,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    gene: Any = _UNSET,
    layer: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw feature-related cumulative density function.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    feature: str
        Feature of interest to plot.
    expression: str, optional
        Layer to use instead of `adata.X`.
    obs: str, optional
        Observation column in `adata.obs` defining groups.
    colors: Colors (optional, default: None)
        Colors used for cumulative density curves.
    legend: bool or mapping (default: True)
        Legend configuration. False disables the legend. True draws the legend
        using default Matplotlib parameters. If a mapping is provided, it is
        forwarded as keyword arguments to `Axes.legend`.
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
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major
          formatter on x- and y-axis

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.
    """
    feature = _resolve_feature_argument(feature, gene, stacklevel=2)
    expression = _resolve_expression_argument(expression, layer, stacklevel=2)

    legend = _resolve_legend_argument(legend, kwargs)

    def _ecdf(values):

        values = np.sort(np.asarray(values))
        y = np.arange(1, len(values) + 1) / len(values)
        return values, y

    counts = pd.DataFrame(
        {"counting": _counts_vector(adata, feature, expression)},
        index=adata.obs.index,
    )

    if obs:
        counts = pd.concat([counts, adata.obs[obs].astype("category")], axis=1)
        if colors is None:
            color_values = [
                gray,
                *CLASSIC_COLORS[
                    1 : len(adata.obs[obs].astype("category").cat.categories) + 1
                ],
            ]
        elif isinstance(colors, Mapping):
            color_values = [
                gray,
                *[
                    colors[cluster]
                    for cluster in adata.obs[obs].astype("category").cat.categories
                ],
            ]
        elif isinstance(colors, Colormap):
            color_values = [gray, *colormap_colors(colors)]
        else:
            color_values = list(colors)
    elif colors is None:
        color_values = [blue]
    elif isinstance(colors, Colormap):
        color_values = colormap_colors(colors)
    else:
        color_values = list(colors)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = figure_from_axes(ax)

    x, y = _ecdf(counts["counting"])
    ax.step(x, y, where="post", color=color_values[0], label="all")
    if obs is not None:
        for _cluster, _color in zip(counts[obs].cat.categories, color_values[1:]):
            _counts = counts.loc[counts[obs] == _cluster]["counting"]
            x, y = _ecdf(_counts)
            ax.step(x, y, where="post", color=_color, label=_cluster)

    if xlabel is not None:
        set_axis_label(ax, "xlabel", xlabel)
    if ylabel is not None:
        set_axis_label(ax, "ylabel", ylabel)

    if min(counts["counting"]) == 0:
        ax.set_xlim(min(counts["counting"]), max(counts["counting"]) * 1.1)
    if obs:
        if legend is True:
            ax.legend()
        elif legend is not False:
            ax.legend(**legend)

    (
        ax.xaxis.set_major_formatter(kwargs["formatter"])
        if "formatter" in kwargs
        else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    )
    (
        ax.yaxis.set_major_formatter(kwargs["formatter"])
        if "formatter" in kwargs
        else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    )

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax


def kde_plot(*args: Any, **kwargs: Any) -> Optional[Tuple[Figure, Axes]]:
    """
    Deprecated alias for `density()`.

    Use `density()` instead. This alias is kept temporarily for backward
    compatibility.
    """

    _warn_deprecated("`kde_plot()`", replacement="`density()`", stacklevel=2)

    return density(*args, **kwargs)


def ecdf_plot(*args: Any, **kwargs: Any) -> Optional[Tuple[Figure, Axes]]:
    """
    Deprecated alias for `cdf()`.

    Use `cdf()` instead. This alias is kept temporarily for backward
    compatibility.
    """

    _warn_deprecated("`ecdf_plot()`", replacement="`cdf()`", stacklevel=2)

    return cdf(*args, **kwargs)
