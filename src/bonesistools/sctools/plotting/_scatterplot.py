#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

from ._typing import RGB
from .._typing import ScData, anndata_or_mudata_checker

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Colormap, ListedColormap
from mpl_toolkits.mplot3d import Axes3D, art3d

Colors = Union[Sequence[object], Iterator[object], Colormap, Mapping[object, object]]

from ..tools._utils import choose_representation

from ._colors import (
    QUALITATIVE_COLORS,
    black,
    generate_colormap,
    gray,
    lightgray,
)
from ..tools import barycenters

import networkx as nx


def __figure_from_axes(ax: Axes) -> Figure:
    return cast(Figure, ax.figure)


def __set_window_title(fig: Figure, title: str) -> None:
    manager = fig.canvas.manager

    if manager is not None:
        manager.set_window_title(title)


def __colormap_colors(colors: Colors) -> Sequence[object]:
    if isinstance(colors, ListedColormap):
        return cast(Sequence[object], colors.colors)

    return cast(Sequence[object], colors)


def __default_plot(plot: Callable[..., Tuple[Figure, Axes]]):

    def wrapper(
        scdata: ScData,  # type: ignore
        obs: str,
        use_rep: str,
        colors: Optional[Colors] = None,
        n_components: Optional[int] = 2,
        ax: Optional[Axes] = None,
        **kwargs: Any,
    ) -> Tuple[Figure, Axes]:

        if obs not in scdata.obs:
            raise KeyError(f"key '{obs}' not found in scdata.obs")
        if use_rep not in scdata.obsm:
            raise KeyError(f"key '{use_rep}' not found in scdata.obsm")

        n_components = 2 if n_components is None else n_components
        fig, ax = cast(
            Tuple[Figure, Axes],
            plot(scdata, obs, use_rep, colors, n_components, ax=ax, **kwargs),
        )

        if "xlabel" in kwargs:
            ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
        if "ylabel" in kwargs:
            ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])
        if "zlabel" in kwargs and n_components > 2:
            cast(Any, ax).set_zlabel(
                "" if kwargs["zlabel"] is None else kwargs["zlabel"]
            )

        if "tick_params" in kwargs:
            ax.tick_params(**kwargs["tick_params"])
        else:
            if "xtick_params" in kwargs:
                ax.tick_params(axis="x", **kwargs["xtick_params"])
            if "ytick_params" in kwargs:
                ax.tick_params(axis="y", **kwargs["ytick_params"])
            if n_components == 3 and "ztick_params" in kwargs:
                cast(Any, ax).tick_params(axis="z", **kwargs["ztick_params"])

        plt.sca(ax)
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
        if n_components == 3:
            ax3d = cast(Any, ax)
            (
                ax3d.zaxis.set_major_formatter(kwargs["formatter"])
                if "formatter" in kwargs
                else ax3d.zaxis.set_major_formatter(FormatStrFormatter("%g"))
            )

        if n_components == 3 and "background_visible" in kwargs:
            if kwargs["background_visible"] is False:
                ax3d = cast(Any, ax)
                ax3d.xaxis.pane.fill = False
                ax3d.yaxis.pane.fill = False
                ax3d.zaxis.pane.fill = False
                ax3d.xaxis.pane.set_edgecolor("w")
                ax3d.yaxis.pane.set_edgecolor("w")
                ax3d.zaxis.pane.set_edgecolor("w")

        return fig, ax

    return wrapper


@__default_plot
def __scatterplot_discrete(
    scdata: ScData,  # type: ignore
    obs: str,
    use_rep: str,
    colors: Optional[Colors] = None,
    n_components: int = 2,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]:

    if "add_legend" in kwargs:
        add_legend = kwargs["add_legend"]
    else:
        add_legend = False

    categories = set(scdata.obs[obs].unique()) & set(
        scdata.obs[obs].astype("category").cat.categories
    )

    if not colors:
        cluster_number = len(scdata.obs[obs].astype("category").cat.categories)
        if len(QUALITATIVE_COLORS) >= cluster_number:
            colors = QUALITATIVE_COLORS[0:cluster_number]
        else:
            colors = generate_colormap(color_number=cluster_number)
    elif isinstance(colors, Mapping):
        colors = [
            colors[cluster]
            for cluster in scdata.obs[obs].astype("category").cat.categories
        ]
    if isinstance(colors, ListedColormap):
        colors = __colormap_colors(colors)

    X = cast(
        np.ndarray,
        choose_representation(scdata, use_rep=use_rep, n_components=n_components),
    )

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(
            111, projection="rectilinear" if n_components == 2 else "3d"
        )
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
        fig.set_figwidth(
            kwargs["figwidth"]
            if "figwidth" in kwargs
            else 6 if n_components == 2 else 8
        )
    else:
        fig = __figure_from_axes(ax)

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}

    if scdata.obs[obs].isna().any():
        idx = scdata.obs[obs].isna()
        if n_components == 2:
            ax.scatter(
                X[idx, 0],
                X[idx, 1],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 3,
                facecolors=(
                    kwargs["nan"]["facecolor"]
                    if "facecolor" in kwargs["nan"]
                    else gray
                ),
                edgecolors=(
                    kwargs["nan"]["edgecolor"]
                    if "edgecolor" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                X[idx, 0],
                X[idx, 1],
                X[idx, 2],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=(
                    kwargs["nan"]["facecolors"]
                    if "facecolors" in kwargs["nan"]
                    else gray
                ),
                edgecolors=(
                    kwargs["nan"]["facecolors"]
                    if "edgecolors" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3,
            )

    color_values = cast(Sequence[object], colors)

    for _cluster, _color in zip(
        scdata.obs[obs].astype("category").cat.categories, color_values
    ):
        _color = cast(Any, _color)

        if _cluster not in categories:
            continue
        if len(_color) == 4:
            if isinstance(_color, np.ndarray):
                _color = _color.tolist()
            del _color[-1]
        idx = np.where(scdata.obs[obs] == _cluster)[0]
        if n_components == 2:
            ax.scatter(
                X[idx, 0],
                X[idx, 1],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=_color,
                edgecolors="none",
                alpha=(
                    kwargs["alpha"]
                    if "alpha" in kwargs and _color != lightgray
                    else 0.3
                ),
                label=_cluster,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                X[idx, 0],
                X[idx, 1],
                X[idx, 2],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=_color,
                edgecolors="none",
                alpha=(
                    kwargs["alpha"]
                    if "alpha" in kwargs and _color != lightgray
                    else 0.3
                ),
                label=_cluster,
            )

    if add_legend:
        box = ax.get_position()
        ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))
        if "legend_params" in kwargs:
            kwargs["lgd_params"] = kwargs["legend_params"]
        if "lgd_params" in kwargs:
            if (
                "loc" not in kwargs["lgd_params"]
                and "bbox_to_anchor" not in kwargs["lgd_params"]
            ):
                if n_components == 3:
                    fig.tight_layout()
                    fig.subplots_adjust(right=0.8)
                kwargs["lgd_params"]["loc"] = "center left"
                kwargs["lgd_params"]["bbox_to_anchor"] = (
                    (1.04, 0.5) if n_components == 2 else (1.09, 0.5)
                )
            else:
                pass
        else:
            kwargs["lgd_params"] = {"loc": "center left", "bbox_to_anchor": (1.04, 0.5)}

        handles, labels = ax.get_legend_handles_labels()
        index = sorted(range(len(labels)), key=lambda idx: labels[idx])
        handles = [handles[i] for i in index]
        labels = [labels[i] for i in index]
        ax.legend(handles, labels, **kwargs["lgd_params"])

    return fig, ax


@__default_plot
def __scatterplot_continuous(
    scdata: ScData,  # type: ignore
    obs: str,
    use_rep: str,
    colors: Optional[Colormap] = None,
    n_components: int = 2,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]:

    cmap = colors if colors is not None else plt.get_cmap("autumn_r")

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}
    if "facecolor" in kwargs["nan"] and not "color" in kwargs["nan"]:
        kwargs["nan"]["color"] = kwargs["nan"]["facecolor"]
    if not np.all(cmap.get_bad() != 0) and not "facecolor" in kwargs["nan"]:
        kwargs["nan"]["color"] = cmap.get_bad()

    X = cast(
        np.ndarray,
        choose_representation(scdata, use_rep=use_rep, n_components=n_components),
    )

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(
            111, projection="rectilinear" if n_components == 2 else "3d"
        )
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
        fig.set_figwidth(
            kwargs["figwidth"]
            if "figwidth" in kwargs
            else 6 if n_components == 2 else 8
        )
    else:
        fig = __figure_from_axes(ax)

    if scdata.obs[obs].isna().any():
        idx = scdata.obs[obs].isna()
        if n_components == 2:
            ax.scatter(
                X[idx, 0],
                X[idx, 1],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 3,
                facecolors=(
                    kwargs["nan"]["facecolor"]
                    if "facecolor" in kwargs["nan"]
                    else lightgray
                ),
                edgecolors=(
                    kwargs["nan"]["edgecolor"]
                    if "edgecolor" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                X[idx, 0],
                X[idx, 1],
                X[idx, 2],
                s=kwargs["s"] if "s" in kwargs else 3,
                facecolors=(
                    kwargs["nan"]["facecolors"]
                    if "facecolors" in kwargs["nan"]
                    else gray
                ),
                edgecolors=(
                    kwargs["nan"]["facecolors"]
                    if "edgecolors" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 0.3,
            )

    if n_components == 2:
        sc = ax.scatter(
            X[:, 0],
            X[:, 1],
            s=kwargs["s"] if "s" in kwargs else 3,
            c=scdata.obs[obs],
            cmap=cmap,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1,
        )
        cb = fig.colorbar(
            sc, shrink=kwargs["colorbar_scale"] if "colorbar_scale" in kwargs else 1
        )
        cb.update_ticks()
    elif n_components == 3:
        ax3d = cast(Any, ax)
        ax3d.scatter3D(
            X[:, 0],
            X[:, 1],
            X[:, 2],
            s=kwargs["s"] if "s" in kwargs else 3,
            c=scdata.obs[obs],
            cmap=cmap,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1,
        )

    return fig, ax


def __add_labels(
    scdata: ScData,  # type: ignore
    obs: str,
    use_rep: str,
    ax: Optional[Axes] = None,
    dim: int = 2,
    **kwargs: Any,
) -> None:

    barycenter_values = cast(
        Dict[object, np.ndarray],
        barycenters(scdata=scdata, obs=obs, use_rep=use_rep),
    )
    keys_to_remove = []
    for k, v in barycenter_values.items():
        if np.isnan(v).any():
            keys_to_remove.append(k)
    for k in keys_to_remove:
        barycenter_values.pop(k, None)

    if ax is None:
        ax = plt.gca()

    if "verticalalignment" not in kwargs:
        kwargs["verticalalignment"] = "center"
    else:
        pass

    if dim == 2:
        for label, value in barycenter_values.items():
            ax.text(x=value[0], y=value[1], s=str(label), **kwargs)
    elif dim == 3:
        ax3d = cast(Any, ax)
        for label, value in barycenter_values.items():
            ax3d.text(x=value[0], y=value[1], z=value[2], s=label, **kwargs)


def __graph_to_plot(
    scdata: ScData,  # type: ignore
    ax: Optional[Axes] = None,
    dim: int = 2,
    z_offset: float = 0.0,
    **kwargs: Any,
) -> None:

    if ax is None:
        ax = plt.gca()

    epg = scdata.uns["epg"]
    flat_tree = scdata.uns["flat_tree"]
    epg_node_pos = nx.get_node_attributes(epg, "pos")

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
            x, y = edge_curve[:, 0], edge_curve[:, 1]
            line = Line2D(
                xdata=x,
                ydata=y,
                color=cast(Any, black),
                **kwargs,
            )
        elif dim == 3:
            x, y, z = edge_curve[:, 0], edge_curve[:, 1], edge_curve[:, 2] + z_offset
            line = art3d.Line3D(xs=x, ys=y, zs=z, color=black, **kwargs)
        else:
            raise ValueError(
                f"invalid argument value for 'dim': "
                f"expected 2 or 3 but received {dim!r}"
            )

        ax.add_line(line)


def __add_labels_to_graph(
    scdata: ScData,  # type: ignore
    ax: Optional[Axes] = None,
    dim: int = 2,
    **kwargs: Any,
) -> None:

    if ax is None:
        ax = plt.gca()

    flat_tree = scdata.uns["flat_tree"]
    flat_tree_node_label = nx.get_node_attributes(flat_tree, "label")
    flat_tree_node_pos = nx.get_node_attributes(flat_tree, "pos")

    if dim == 2:
        for node in flat_tree.nodes:
            ax.text(
                x=flat_tree_node_pos[node][0],
                y=flat_tree_node_pos[node][1],
                s=flat_tree_node_label[node],
                **kwargs,
            )
    elif dim == 3:
        ax3d = cast(Any, ax)
        for node in flat_tree.nodes:
            ax3d.text(
                x=flat_tree_node_pos[node][0],
                y=flat_tree_node_pos[node][1],
                z=flat_tree_node_pos[node][2],
                s=flat_tree_node_label[node],
                **kwargs,
            )


@anndata_or_mudata_checker
def embedding_plot(
    scdata: ScData,  # type: ignore
    obs: str,
    use_rep: str,
    colors: Optional[Colormap] = None,
    n_components: int = 2,
    title: Optional[Union[str, dict[str, Any]]] = None,
    add_labels: bool = False,
    add_graph: bool = False,
    add_labels_to_graph: bool = False,
    automatic_resize: bool = False,
    default_parameters: Optional[Callable[[], None]] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw a scatterplot from an embedding stored in `scdata.obsm`.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Observation column in `scdata.obs` defining groups or values.
    use_rep: str
        Representation key in `scdata.obsm`.
    colors: matplotlib.colors.Colormap, optional
        Colormap or colors used to draw observations.
    n_components: 2 or 3 (default: 2)
        Number of plotted dimensions.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    add_labels: bool (default: False)
        Add labels retrieved from `scdata.obs[obs]` to the embedding.
    add_graph: bool (default: False)
        Draw the elastic principal graph stored in `scdata.uns["epg"]`.
        This is useful for STREAM outputs [1].
    add_labels_to_graph: bool (default: False)
        Add node labels from the elastic principal graph.
    automatic_resize: bool (default: False)
        Resize figure to accommodate large legends.
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
        - zlabel[str]: set the label for the z-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major
          formatter on x-, y- and z-axis
        - add_legend[bool]: when .obs['obs'] are discrete values, specify
          whether to draw legend
        - lgd_params[dict]: when add_legend is True, modify legend following
          the syntax of matplotlib.pyplot.legend
        - tick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines following the syntax of matplotlib.axes.Axes.tick_params
        - xtick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on x-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - ytick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on y-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - ztick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on z-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - text[dict]: change the appearance of text in figure following the
          syntax of matplotlib.text
        - graph[dict]: change the appearance of elastic principal graph lines
        - graph_z_offset[float]: vertical offset applied to graph lines in 3D
        - background_visible[bool]: specify if background color is visible or
          not in case of 3D plotting

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.

    Raises
    ------
    ValueError
        If `n_components` is not 2 or 3.
    TypeError
        If `title` is neither a string nor a dictionary.

    References
    ----------
    [1] Chen et al. (2019). Single-cell trajectories reconstruction,
    exploration and mapping of omics data with STREAM. Nature communications,
    10(1), 1903 (https://www.nature.com/articles/s41467-019-09670-4)
    """

    if n_components not in [2, 3]:
        raise ValueError(
            f"invalid argument value for 'n_components': "
            f"expected 2 or 3 but received {n_components!r}"
        )

    component_number = n_components

    if pd.api.types.is_float_dtype(scdata.obs[obs]):
        fig, ax = __scatterplot_continuous(
            scdata, obs, use_rep, colors, component_number, ax=ax, **kwargs
        )
    elif (
        pd.api.types.is_integer_dtype(scdata.obs[obs])
        or pd.api.types.is_bool_dtype(scdata.obs[obs])
        or pd.api.types.is_string_dtype(scdata.obs[obs])
        or isinstance(scdata.obs[obs].dtype, pd.CategoricalDtype)
    ):
        fig, ax = __scatterplot_discrete(
            scdata, obs, use_rep, colors, component_number, ax=ax, **kwargs
        )
    else:
        raise TypeError(
            f"unsupported dtype for observation {obs!r}: "
            f"{scdata.obs[obs].dtype!r}"
        )

    if add_labels:
        _kwargs = {} if "text" not in kwargs else kwargs["text"]
        __add_labels(
            scdata, obs=obs, use_rep=use_rep, ax=ax, dim=component_number, **_kwargs
        )

    if add_graph:
        _kwargs = {"linewidth": 2.5, "zorder": 10}
        _kwargs.update({} if "graph" not in kwargs else kwargs["graph"])
        __graph_to_plot(
            scdata,
            ax=ax,
            dim=component_number,
            z_offset=kwargs["graph_z_offset"] if "graph_z_offset" in kwargs else 0.0,
            **_kwargs,
        )

    if add_labels_to_graph:
        _kwargs = (
            {"verticalalignment": "bottom" if add_graph else "center"}
            if "text" not in kwargs
            else kwargs["text"]
        )
        __add_labels_to_graph(scdata, ax=ax, dim=component_number, **_kwargs)

    if automatic_resize:
        fig.tight_layout(pad=1.2, rect=(0, 0, 0.84, 1))

    if default_parameters:
        default_parameters()

    if title:
        if isinstance(title, str):
            __set_window_title(fig, title)
            ax.set_title(title)
        elif isinstance(title, dict):
            __set_window_title(fig, title["label"])
            ax.set_title(**title)
        else:
            raise TypeError(
                f"unsupported argument type for 'title': "
                f"expected {str} or {dict} but received {type(title)}"
            )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax
