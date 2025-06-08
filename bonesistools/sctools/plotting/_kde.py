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
from anndata import AnnData
from ._typing import RGB
from .._typing import anndata_checker

import pandas as pd

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
    colors: Optional[Colors] = None,
    default_parameters: Optional[Callable] = None,
    **kwargs: Mapping[str, Any]
) -> Tuple[Figure, Axes]:
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
    colors: Colors (optional, default: None)
        Density function are colored with respect to a list of color values.
    default_parameters: Callable (optional, default: None)
        Function specifying default figure parameters.
    **kwargs
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major formatter on x- and y-axis

    Returns
    -------
    Create matplotlib figure and axe.
    """

    import seaborn as sns

    counts = adata[:,gene].layers[layer] if layer else adata[:,gene].X
    if scipy.sparse.issparse(counts):
        counts = pd.Series(counts.toarray().squeeze(), index=adata.obs.index, name="counting")
    if obs:
        counts = pd.concat([counts, adata.obs[obs].astype("category")], axis=1)
        if not colors:
            colors = [_colors.gray, *_colors.COLORS[1:len(adata.obs[obs].astype("category").cat.categories)+1]]
    elif not colors:
        colors = [_colors.blue]
    
    fig, ax = plt.subplots()
    sns.kdeplot(
        data=counts["counting"],
        ax=ax,
        color=colors[0],
        fill=True,
        label="all"
    )
    if obs:
        for _cluster, _color in zip(sorted(counts[obs].cat.categories), colors[1:]):
            _counts = counts.loc[counts[obs] == _cluster]["counting"]
            sns.kdeplot(
                data=_counts,
                ax=ax,
                color=_color,
                fill=False,
                label=_cluster
            )

    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else fig.get_figheight())
    fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else fig.get_figwidth())

    if "xlabel" in kwargs:
        ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
    if "ylabel" in kwargs:
        ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])

    if min(counts["counting"]) == 0:
        plt.xlim(min(counts["counting"]),max(counts["counting"])*1.1)
    if obs:
        ax.legend(loc='upper right')

    ax.xaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    
    if default_parameters:
        default_parameters()
    
    return fig, ax
