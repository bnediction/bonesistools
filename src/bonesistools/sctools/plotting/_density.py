#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Mapping,
    Sequence,
    Tuple,
    Callable,
    Any
)
from pathlib import Path
from anndata import AnnData
from ._typing import RGB
from .._typing import anndata_checker

import pandas as pd
import numpy as np

import scipy

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.ticker import FormatStrFormatter
from itertools import cycle
from . import _colors

Colors = Union[Sequence[RGB], cycle]

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
    title: Optional[Union[str, dict]] = None,
    default_parameters: Optional[Callable] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    **kwargs: Mapping[str, Any]
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw gene-related density function using kernel density estimation.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    gene: str
        Gene of interest for which display the density.
    layer: str (optional, default: None)
        If specified, the counting are retrieved from adata.layers['layer'], otherwise from adata.X.
    obs: str (optional, default: None)
        If specified, draw gene-related density w.r.t. categories.
        The classification is retrieved by adata.obs['obs'], which must be categorical/qualitative values.
    not_all: bool (default: False)
        If True, do not draw density function using all barcodes (raise Error if True and 'obs' is not specified).
    clip: bool (default: False)
        If True, clip density between the minimum value and the quantile at 99%.
    colors: Colors (optional, default: None)
        Density function are colored with respect to a list of color values.
    title: Union[str, dict] (optional, default: None)
        Add title to current axe.
    default_parameters: Callable (optional, default: None)
        Function specifying default figure parameters.
    outfile: Path (optional, default: None)
        If specified, save the figure.
    **kwargs
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major formatter on x- and y-axis

    Returns
    -------
    Depending on 'outfile', save figure or create matplotlib Figure and Axes object.
    """

    if obs is None and not_all is True:
        raise ValueError(f"invalid argument value for 'obs' and 'not_all': expected 'not_all' assigned to False when 'obs' is not specified")

    import seaborn as sns

    counts = adata[:,gene].layers[layer] if layer else adata[:,gene].X
    if scipy.sparse.issparse(counts):
        counts = pd.Series(counts.toarray().squeeze(), index=adata.obs.index, name="counting")
    if obs:
        counts = pd.concat([counts, adata.obs[obs].astype("category")], axis=1)

    if obs:
        if not colors:
            cluster_number = len(adata.obs[obs].astype("category").cat.categories)
            if len(_colors.QUALITATIVE_COLORS) >= cluster_number:
                colors = _colors.QUALITATIVE_COLORS[0:cluster_number]
            else:
                colors = _colors.generate_colormap(color_number=cluster_number)
        elif isinstance(colors, Mapping):
            colors = [colors[cluster] for cluster in adata.obs[obs].astype("category").cat.categories]
        if hasattr(colors, "colors"):
            colors = colors.colors
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    q = np.quantile(counts["counting"], 0.99)

    if not_all is False:
        sns.kdeplot(
            data=counts["counting"],
            ax=ax,
            color=_colors.gray,
            fill=True,
            clip=(min(counts["counting"]), q) if clip is True else None,
            label="all"
        )

    if obs:
        for _cluster, _color in zip(adata.obs[obs].astype("category").cat.categories, colors):
            if len(_color) == 4:
                if isinstance(_color, np.ndarray):
                    _color = _color.tolist()
                del _color[-1]
            sns.kdeplot(
                data=counts.loc[counts[obs] == _cluster]["counting"],
                ax=ax,
                color=_color,
                fill=False,
                clip=(min(counts["counting"]), q) if clip is True else None,
                label=_cluster
            )

    if "xlabel" in kwargs:
        ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
    if "ylabel" in kwargs:
        ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])
    
    if title:
        if isinstance(title, str):
            fig.canvas.manager.set_window_title(title)
            ax.set_title(title)
        elif isinstance(title, dict):
            fig.canvas.manager.set_window_title(title["label"])
            ax.set_title(**title)

    if obs and show_legend:
        ax.legend(loc="upper right")

    ax.xaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    
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
    colors = None,
    show_legend: bool = True,
    default_parameters: Optional[Callable] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    **kwargs: Mapping[str, Any]
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw gene-related cumulative density function.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    gene: str
        Gene of interest for which display the density.
    layer: str (optional, default: None)
        If specified, the counting are retrieved from adata.layers['layer'], otherwise from adata.X.
    obs: str (optional, default: None)
        If specified, draw gene-related density w.r.t. categories.
        The classification is retrieved by adata.obs['obs'], which must be categorical/qualitative values.
    colors: Colors (optional, default: None)
        Density function are colored with respect to a list of color values.
    default_parameters: Callable (optional, default: None)
        Function specifying default figure parameters.
    outfile: Path (optional, default: None)
        If specified, save the figure.
    **kwargs
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major formatter on x- and y-axis

    Returns
    -------
    Depending on 'outfile', save figure or create matplotlib Figure and Axes object.
    """

    counts = adata[:,gene].layers[layer] if layer else adata[:,gene].X
    if scipy.sparse.issparse(counts):
        counts = pd.Series(counts.toarray().squeeze(), index=adata.obs.index, name="counting")
    if obs:
        counts = pd.concat([counts, adata.obs[obs].astype("category")], axis=1)
        if not colors:
            colors = [_colors.gray, *_colors.COLORS[1:len(adata.obs[obs].astype("category").cat.categories)+1]]
        elif isinstance(colors, Mapping):
            colors = [_colors.gray, *[colors[cluster] for cluster in adata.obs[obs].astype("category").cat.categories]]
    elif not colors:
        colors = [_colors.blue]
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    ax.ecdf(
        counts["counting"],
        color=colors[0],
        label="all"
    )
    if obs:
        for _cluster, _color in zip(counts[obs].cat.categories, colors[1:]):
            _counts = counts.loc[counts[obs] == _cluster]["counting"]
            ax.ecdf(
                _counts,
                color=_color,
                label=_cluster
            )

    if "xlabel" in kwargs:
        ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
    if "ylabel" in kwargs:
        ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])

    if min(counts["counting"]) == 0:
        ax.set_xlim(min(counts["counting"]), max(counts["counting"])*1.1)
    if obs and show_legend:
        ax.legend(loc="upper right")

    ax.xaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    
    if default_parameters:
        default_parameters()
    
    if outfile:
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        return fig, ax
