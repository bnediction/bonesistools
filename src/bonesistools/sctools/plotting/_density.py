#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from .._typing import anndata_checker
from ._colors import (
    COLORS,
    QUALITATIVE_COLORS,
    blue,
    generate_colormap,
    gray,
)

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]


def _counts_vector(adata: AnnData, gene: str, layer: Optional[str]) -> np.ndarray:
    counts = adata[:, gene].layers[layer] if layer else adata[:, gene].X

    if scipy.sparse.issparse(counts):
        counts = cast(Any, counts).toarray()

    return np.asarray(counts).squeeze()


def _figure_from_axes(ax: Axes) -> Figure:
    return cast(Figure, ax.figure)


def _set_window_title(fig: Figure, title: str) -> None:
    manager = fig.canvas.manager

    if manager is not None:
        manager.set_window_title(title)


def _colormap_colors(colors: Colors) -> Sequence[object]:
    if isinstance(colors, ListedColormap):
        return cast(Sequence[object], colors.colors)

    return cast(Sequence[object], colors)


def _normalize_color(color: object) -> object:
    if isinstance(color, str):
        return color

    normalized_color = (
        color.tolist() if isinstance(color, np.ndarray) else list(cast(Any, color))
    )
    if len(normalized_color) == 3:
        normalized_color.append(1)

    return normalized_color


@anndata_checker
def kde_plot(
    adata: AnnData,
    gene: str,
    layer: Optional[str] = None,
    obs: Optional[str] = None,
    not_all: bool = False,
    clip: bool = False,
    colors: Optional[Colors] = None,
    show_legend: bool = True,
    title: Optional[Union[str, dict[str, Any]]] = None,
    default_parameters: Optional[Callable[[], None]] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw gene-related density function using kernel density estimation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    gene: str
        Gene of interest to plot.
    layer: str, optional
        Layer to use instead of `adata.X`.
    obs: str, optional
        Observation column in `adata.obs` defining groups.
    not_all: bool (default: False)
        If True, do not draw density function using all barcodes. Raises an
        error if True and `obs` is not specified.
    clip: bool (default: False)
        If True, clip density between the minimum value and the quantile at 99%.
    colors: Colors (optional, default: None)
        Colors used for density curves.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    default_parameters: Callable (optional, default: None)
        Function specifying default figure parameters.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
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
        If `not_all=True` while `obs` is None.
    KeyError
        If `obs` is specified but not found in `adata.obs`.
    """

    if obs is None and not_all is True:
        raise ValueError(
            "invalid argument values for 'obs' and 'not_all': "
            "expected not_all=False when obs is None"
        )

    import seaborn as sns

    counts = pd.DataFrame(
        {"counting": _counts_vector(adata, gene, layer)},
        index=adata.obs.index,
    )

    if obs:
        if obs not in adata.obs:
            raise KeyError(f"key {obs!r} not found in adata.obs")
        counts[obs] = adata.obs[obs]
        if not colors:
            cluster_number = len(adata.obs[obs].astype("category").cat.categories)
            if len(QUALITATIVE_COLORS) >= cluster_number:
                colors = QUALITATIVE_COLORS[0:cluster_number]
            else:
                colors = generate_colormap(color_number=cluster_number)
        elif isinstance(colors, Mapping):
            colors = [
                colors[cluster]
                for cluster in adata.obs[obs].astype("category").cat.categories
            ]
        if isinstance(colors, Colormap):
            colors = _colormap_colors(colors)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = _figure_from_axes(ax)

    q = np.quantile(counts["counting"], 0.99)
    clip_range = cast(
        Optional[Tuple[float, float]],
        (float(min(counts["counting"])), float(q)) if clip is True else None,
    )

    if not_all is False:
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
            _color = _normalize_color(_color)
            sns.kdeplot(
                data=cast(Any, counts.loc[counts[obs] == _cluster]["counting"]),
                ax=ax,
                color=cast(Any, _color),
                fill=False,
                clip=clip_range,
                label=_cluster,
            )

    if "xlabel" in kwargs:
        ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
    if "ylabel" in kwargs:
        ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])

    if title:
        if isinstance(title, str):
            _set_window_title(fig, title)
            ax.set_title(title)
        elif isinstance(title, dict):
            _set_window_title(fig, title["label"])
            ax.set_title(**title)

    if obs and show_legend:
        ax.legend(loc="upper right")

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

    if default_parameters:
        default_parameters()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax


def ecdf_plot(
    adata: AnnData,
    gene: str,
    layer: Optional[str] = None,
    obs: Optional[str] = None,
    colors=None,
    show_legend: bool = True,
    default_parameters: Optional[Callable[[], None]] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw gene-related cumulative density function.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    gene: str
        Gene of interest to plot.
    layer: str, optional
        Layer to use instead of `adata.X`.
    obs: str, optional
        Observation column in `adata.obs` defining groups.
    colors: Colors (optional, default: None)
        Colors used for cumulative density curves.
    default_parameters: Callable (optional, default: None)
        Function specifying default figure parameters.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major
          formatter on x- and y-axis

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.
    """

    def _ecdf(values):
        values = np.sort(np.asarray(values))
        y = np.arange(1, len(values) + 1) / len(values)
        return values, y

    counts = pd.DataFrame(
        {"counting": _counts_vector(adata, gene, layer)},
        index=adata.obs.index,
    )

    if obs:
        counts = pd.concat([counts, adata.obs[obs].astype("category")], axis=1)
        if not colors:
            colors = [
                gray,
                *COLORS[1 : len(adata.obs[obs].astype("category").cat.categories) + 1],
            ]
        elif isinstance(colors, Mapping):
            colors = [
                gray,
                *[
                    colors[cluster]
                    for cluster in adata.obs[obs].astype("category").cat.categories
                ],
            ]
    elif not colors:
        colors = [blue]

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = _figure_from_axes(ax)

    x, y = _ecdf(counts["counting"])
    ax.step(x, y, where="post", color=colors[0], label="all")
    if obs is not None:
        for _cluster, _color in zip(counts[obs].cat.categories, colors[1:]):
            _counts = counts.loc[counts[obs] == _cluster]["counting"]
            x, y = _ecdf(_counts)
            ax.step(x, y, where="post", color=_color, label=_cluster)

    if "xlabel" in kwargs:
        ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
    if "ylabel" in kwargs:
        ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])

    if min(counts["counting"]) == 0:
        ax.set_xlim(min(counts["counting"]), max(counts["counting"]) * 1.1)
    if obs and show_legend:
        ax.legend(loc="upper right")

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

    if default_parameters:
        default_parameters()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax
