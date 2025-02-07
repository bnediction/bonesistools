#!/usr/bin/env python

from collections.abc import Mapping
import types
from typing import (
    Optional,
    Union,
    Sequence,
    Tuple
)

from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Colormap
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from itertools import cycle
Colors = Union[Sequence[Tuple[str, str, str]], cycle, Colormap]

from . import _figure
from . import _colors

from .._check_anndata import _adata_arg_checking

def __default_plot(
    plot: types.FunctionType
):

    def wrapper(
        adata: ad.AnnData,
        obs: str,
        obsm: str,
        colors: Optional[Colors] = None,
        n_components: Optional[int] = 2,
        **kwargs
    ):

        if obs not in adata.obs:
            raise KeyError(f"adata.obs[`{obs}`] does not exist.")
        if obsm not in adata.obsm:
            raise KeyError(f"adata.obsm[`{obsm}`] does not exist.")

        fig, ax = plot(
            adata,
            obs,
            obsm,
            colors,
            n_components,
            **kwargs
        )

        if "xlabel" in kwargs:
            ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
        if "ylabel" in kwargs:
            ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])
        if "zlabel" in kwargs and n_components > 2:
            ax.set_zlabel("" if kwargs["zlabel"] is None else kwargs["zlabel"])
        
        if "tick_params" in kwargs:
            ax.tick_params(**kwargs["tick_params"])
        else:
            if "xtick_params" in kwargs:
                ax.tick_params(axis="x", **kwargs["xtick_params"])
            if "ytick_params" in kwargs:
                ax.tick_params(axis="y", **kwargs["ytick_params"])
            if n_components ==3 and "ztick_params" in kwargs:
                ax.tick_params(axis="z", **kwargs["ztick_params"])

        plt.sca(ax)
        ax.xaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.yaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        if n_components == 3:
            ax.zaxis.set_major_formatter(kwargs["formatter"]) if "formatter" in kwargs else ax.zaxis.set_major_formatter(FormatStrFormatter("%g"))

        if n_components == 3 and "background_visible" in kwargs:
            if kwargs["background_visible"] is False:
                ax.xaxis.pane.fill = False
                ax.yaxis.pane.fill = False
                ax.zaxis.pane.fill = False
                ax.xaxis.pane.set_edgecolor("w")
                ax.yaxis.pane.set_edgecolor("w")
                ax.zaxis.pane.set_edgecolor("w")

        return fig, ax
    
    return wrapper

@__default_plot
def __scatterplot_discrete(
    adata: ad.AnnData,
    obs: str,
    obsm: str,
    colors: Optional[Colors] = None,
    n_components: Optional[int] = 2,
    **kwargs
):

    if len(adata.obs[obs].unique()) < 2:
        raise ValueError(f"adata.obs[{obs}] specifies only one category, aborting")
    elif "add_legend" in kwargs:
        add_legend = kwargs["add_legend"]
    else:
        add_legend = False
    
    if not colors:
        colors = _colors.COLORS[0:len(adata.obs[obs].astype("category").cat.categories)]
    elif isinstance(colors, Mapping):
        colors = [colors[cluster] for cluster in adata.obs[obs].astype("category").cat.categories]

    fig = plt.figure()
    ax = plt.axes(projection = "rectilinear" if n_components == 2 else "3d")
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
    fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else 5)

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}

    if adata.obs[obs].isna().any():
        idx = adata.obs[obs].isna()
        if n_components==2:
            ax.scatter(
                adata.obsm[obsm][idx,0],
                adata.obsm[obsm][idx,1],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 3,
                facecolors=kwargs["nan"]["facecolors"] if "facecolors" in kwargs["nan"] else _colors.gray,
                edgecolors=kwargs["nan"]["facecolors"] if "edgecolors" in kwargs["nan"] else "none",
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3
            )
        elif n_components==3:
            ax.scatter3D(
                adata.obsm[obsm][idx,0],
                adata.obsm[obsm][idx,1],
                adata.obsm[obsm][idx,2],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=kwargs["nan"]["facecolors"] if "facecolors" in kwargs["nan"] else _colors.gray,
                edgecolors=kwargs["nan"]["facecolors"] if "edgecolors" in kwargs["nan"] else "none",
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3
            )
    
    for _cluster, _color in zip(adata.obs[obs].astype("category").cat.categories, colors):
        idx = np.where(adata.obs[obs] == _cluster)[0]
        if n_components==2:
            ax.scatter(
                adata.obsm[obsm][idx,0],
                adata.obsm[obsm][idx,1],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=_color,
                edgecolors="none",
                alpha=kwargs["alpha"] if "alpha" in kwargs and _color != _colors.lightgray else 0.3,
                label=_cluster
            )
        elif n_components==3:
            ax.scatter3D(
                adata.obsm[obsm][idx,0],
                adata.obsm[obsm][idx,1],
                adata.obsm[obsm][idx,2],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=_color,
                edgecolors="none",
                alpha=kwargs["alpha"] if "alpha" in kwargs and _color != _colors.lightgray else 0.3,
                label=_cluster,
            )

    if add_legend:
        fig.set_figwidth(kwargs["figwidth"]*1.25 if "figwidth" in kwargs else 6.25)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        if "legend_params" in kwargs:
            kwargs["lgd_params"] = kwargs["legend_params"]
        if "lgd_params" in kwargs:
            if "loc" not in kwargs["lgd_params"] and "bbox_to_anchor" not in kwargs["lgd_params"]:
                if n_components == 3:
                    fig.tight_layout()
                    fig.subplots_adjust(right=0.8)
                kwargs["lgd_params"]["loc"] = "center left"
                kwargs["lgd_params"]["bbox_to_anchor"] = (1.04, 0.5) if n_components == 2 else (1.09, 0.5)
            else:
                pass
        else:
            kwargs["lgd_params"] = {
                "loc":"center left",
                "bbox_to_anchor":(1.04, 0.5)
            }
    
        handles, labels = ax.get_legend_handles_labels()
        index = sorted(range(len(labels)), key=lambda idx: labels[idx])
        handles = [handles[i] for i in index]
        labels = [labels[i] for i in index]
        ax.legend(
            handles,
            labels,
            **kwargs["lgd_params"]
        )

    return fig, ax

@__default_plot
def __scatterplot_continuous(
    adata: ad.AnnData,
    obs: str,
    obsm: str,
    colors: Optional[Colormap] = None,
    n_components: Optional[int] = 2,
    **kwargs
):

    if colors:
        if hasattr(colors, "name"):
            _cmap = colors.name
        else:
            _cmap = colors
    else:
        _cmap = "autumn"

    fig = plt.figure()
    ax = plt.axes(projection = "rectilinear" if n_components == 2 else "3d")
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
    fig.set_figwidth(kwargs["figwidth"] if "figwidth" in kwargs else 5 if n_components == 2 else 6)
    
    if n_components == 2:
        sc = ax.scatter(
            adata.obsm[obsm][:,0],
            adata.obsm[obsm][:,1],
            s=kwargs["s"] if "s" in kwargs else 3,
            c=adata.obs[obs],
            cmap=_cmap,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1
        )
        cb = fig.colorbar(sc, shrink=kwargs["colorbar_scale"] if "colorbar_scale" in kwargs else 1)
        cb.update_ticks()
    elif n_components==3:
        ax.scatter3D(
            adata.obsm[obsm][:,0],
            adata.obsm[obsm][:,1],
            adata.obsm[obsm][:,2],
            s=kwargs["s"] if "s" in kwargs else 3,
            c=adata.obs[obs],
            cmap=_cmap,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1
        )

    return fig, ax

@_adata_arg_checking
def __graph_to_plot(
    adata: ad.AnnData,
    ax: Optional[Axes] = None,
    dim: Optional[int] = 2,
    **kwargs
    ):

    if "nx" not in locals():
        import networkx as nx

    if ax is None:
        ax = plt.gca()

    epg = adata.uns["epg"]
    flat_tree = adata.uns["flat_tree"]
    epg_node_pos = nx.get_node_attributes(epg,"pos")

    traces = set()
    for node in flat_tree:
        for trace in flat_tree.adj[node].values():
            traces.add(tuple(trace["nodes"]))

    edge_curves = list()
    for trace in traces:
        _edge_curve = list()
        for node in trace:
            _edge_curve.append(np.array([epg_node_pos[node]]))
        edge_curves.append(np.concatenate(_edge_curve))

    for edge_curve in edge_curves:
        if dim == 2:
            x, y = edge_curve[:,0], edge_curve[:, 1]
            line = plt.Line2D(xdata=x, ydata=y, color=_colors.black, **kwargs)
        elif dim == 3:
            x, y, z = edge_curve[:,0], edge_curve[:, 1], edge_curve[:, 2]
            line = mplot3d.art3d.Line3D(xs=x, ys=y, zs=z, color=_colors.black, **kwargs)
        ax.add_line(line)

@_adata_arg_checking
def __text_to_plot(
    adata: ad.AnnData,
    ax: Optional[Axes] = None,
    dim: Optional[int] = 2,
    **kwargs
    ):

    if "nx" not in locals():
        import networkx as nx

    if ax is None:
        ax = plt.gca()

    flat_tree = adata.uns["flat_tree"]
    flat_tree_node_label = nx.get_node_attributes(flat_tree, "label")
    flat_tree_node_pos = nx.get_node_attributes(flat_tree, "pos")

    if dim == 2:
        for node in flat_tree.nodes:
            plt.text(
                x=flat_tree_node_pos[node][0],
                y=flat_tree_node_pos[node][1],
                s=flat_tree_node_label[node],
                **kwargs
            )
    elif dim == 3:
        for node in flat_tree.nodes:
            ax.text(
                x=flat_tree_node_pos[node][0],
                y=flat_tree_node_pos[node][1],
                z=flat_tree_node_pos[node][2],
                s=flat_tree_node_label[node],
                **kwargs
            )

@_adata_arg_checking
def embedding_plot(
    adata: ad.AnnData,
    obs: str,
    obsm: str,
    colors: Optional[Colormap] = None,
    n_components: Optional[int] = 2,
    outfile: Optional[Path] = None,
    add_graph: Optional[bool] = None,
    add_text: Optional[bool] = None,
    default_parameters: Optional[types.FunctionType] = None,
    **kwargs
) -> Tuple[Figure, Axes]:
    """
    Draw a scatterplot between the `n_components` first columns of .obsm[`obsm`]
    by using a classification/clusterization with respect to .obs[`obs`].

    Parameters
    ----------
    adata
        Annotated data matrix
    obs
        The classification is retrieved by .obs[`obs`], which must be categorical/qualitative values
    obsm
        The data points are retrieved by the first columns in .obsm[`obsm`]
    colors
        Visualization of the mapping from a list of color values
    n_components
        Number of plotted dimensions (default: 2)
    outfile
        If specified, save the figure
    add_graph
        Draw elastic principal graph
    add_text
        Add node labels of elastic principal graph
    default_parameters
        Function specifying default figure parameters
    **kwargs
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - zlabel[str]: set the label for the z-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major formatter on x-, y- and z-axis
        - add_legend[bool]: when .obs[`obs`] are discrete values, specify whether to draw legend
        - lgd_params[dict]: when add_legend is True, modify legend following the syntax of matplotlib.pyplot.legend
        - tick_params[dict]: change the appearance of ticks, tick labels, and gridlines following the syntax of matplotlib.axes.Axes.tick_params
        - xtick_params[dict]: change the appearance of ticks, tick labels, and gridlines on x-axis following the syntax of matplotlib.axes.Axes.tick_params
        - ytick_params[dict]: change the appearance of ticks, tick labels, and gridlines on y-axis following the syntax of matplotlib.axes.Axes.tick_params
        - ztick_params[dict]: change the appearance of ticks, tick labels, and gridlines on z-axis following the syntax of matplotlib.axes.Axes.tick_params
        - text[dict]: change the appearance of text in figure following the syntax of matplotlib.text
        - background_visible[bool]: specify if background color is visible or not in case of 3D plotting

    Returns
    -------
    Depending on `outfile`, save figure or create figure and axe.
    """

    if n_components not in [2,3]:
        raise ValueError(f"`n_components` parameter value is {n_components}, please set it to 2 or 3.")
    elif n_components == 3 and adata.obsm[obsm].shape[1] < 3:
        raise ValueError(f"incoherent value: `n_components` parameter value is {n_components} while number of dimension in .obsm[{obsm}] is {adata.obsm[obsm].shape[1]}.")

    if pd.api.types.is_float_dtype(adata.obs[obs]):
        fig, ax = __scatterplot_continuous(
            adata,
            obs,
            obsm,
            colors,
            n_components,
            **kwargs
        )
    elif pd.api.types.is_integer_dtype(adata.obs[obs]) or \
         pd.api.types.is_bool_dtype(adata.obs[obs]) or \
         pd.api.types.is_string_dtype(adata.obs[obs]) or \
         pd.api.types.is_categorical_dtype(adata.obs[obs]):
        fig, ax = __scatterplot_discrete(
            adata,
            obs,
            obsm,
            colors,
            n_components,
            **kwargs
        )
    
    if add_graph:
        ax = plt.gca()
        __graph_to_plot(adata, ax=ax, dim=n_components)
    
    if add_text:
        ax = plt.gca()
        if "text" not in kwargs:
            kwargs["text"] = dict()
            kwargs["text"]["verticalalignment"] = "bottom" if add_graph else "center"
        elif "verticalalignment" not in kwargs["text"]:
            kwargs["text"]["verticalalignment"] = "bottom" if add_graph else "center"
        else:
            pass
        __text_to_plot(
            adata,
            ax,
            dim=n_components,
            **kwargs["text"]
        )
    
    if default_parameters:
        default_parameters()
    
    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax