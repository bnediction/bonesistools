#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from collections.abc import Sequence as SequenceABC
from numbers import Number
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    Mapping,
    Optional,
    Tuple,
    Union,
    cast,
    overload,
)

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from anndata import AnnData
from matplotlib.axes._axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import art3d

from ..._warnings import _warn_deprecated
from .._typing import ScData, anndata_checker, anndata_or_mudata_checker
from ..tools import extract_paga_graph
from ._colors import black
from ._scatterplot import Colors, embedding


@anndata_or_mudata_checker
def graph_overlay(
    scdata: ScData,  # type: ignore
    graph_key: str = "epg",
    ax: Optional[Axes] = None,
    n_components: int = 2,
    show_labels: bool = False,
    z_offset: float = 0.0,
    label_key: str = "label",
    label_kwargs: Optional[Dict[str, Any]] = None,
    **kwargs: Any,
) -> Axes:
    """
    Overlay a principal graph on an existing embedding axes.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    graph_key: str (default: "epg")
        Key pointing to the principal graph stored in `scdata.uns[graph_key]`.
    ax: Axes, optional
        Matplotlib axes containing the embedding. If None, use current axes.
    n_components: 2 or 3 (default: 2)
        Number of plotted dimensions.
    show_labels: bool (default: False)
        Draw graph node labels.
    z_offset: float (default: 0.0)
        Vertical offset applied to graph lines in 3D.
    label_key: str (default: "label")
        Node attribute used as label when `show_labels=True`.
    label_kwargs: dict, optional
        Keyword arguments passed to Matplotlib text drawing.
    **kwargs: Any
        Keyword arguments passed to graph line drawing.

    Returns
    -------
    Axes
        Matplotlib axes with the graph overlay.

    Raises
    ------
    KeyError
        If `graph_key` or graph node positions are missing.
    ValueError
        If `n_components` is not 2 or 3.

    Chen et al. (2019). Single-cell trajectories reconstruction, exploration
    and mapping of omics data with STREAM. Nature Communications, 10(1), 1903.
    """

    if n_components not in (2, 3):
        raise ValueError(
            f"invalid argument value for 'n_components': "
            f"expected 2 or 3 but received {n_components!r}"
        )

    if ax is None:
        ax = plt.gca()

    if graph_key not in scdata.uns:
        raise KeyError(f"key {graph_key!r} not found in scdata.uns")

    graph = cast(nx.Graph, scdata.uns[graph_key])
    positions = cast(Mapping[object, np.ndarray], nx.get_node_attributes(graph, "pos"))

    if not positions:
        raise KeyError(f"node positions not found in scdata.uns[{graph_key!r}]")

    if (
        graph_key == "epg"
        and graph.number_of_edges() == 0
        and "flat_tree" in scdata.uns
    ):
        edge_graph = cast(nx.Graph, scdata.uns["flat_tree"])
    else:
        edge_graph = graph

    line_kwargs = {"linewidth": 2.5, "zorder": 10}
    line_kwargs.update(kwargs)
    edge_curves = []

    for source, target, edge_data in edge_graph.edges(data=True):
        edge_nodes = cast(Mapping[str, Any], edge_data).get("nodes", (source, target))
        curve_points = [
            np.asarray(positions[node])
            for node in cast(Iterable[object], edge_nodes)
            if node in positions
        ]

        if len(curve_points) >= 2:
            edge_curves.append(np.asarray(curve_points))

    for edge_curve in edge_curves:
        if n_components == 2:
            x, y = edge_curve[:, 0], edge_curve[:, 1]
            line = Line2D(xdata=x, ydata=y, color=cast(Any, black), **line_kwargs)
        else:
            x, y = edge_curve[:, 0], edge_curve[:, 1]
            z = edge_curve[:, 2] + z_offset
            line = art3d.Line3D(xs=x, ys=y, zs=z, color=black, **line_kwargs)

        ax.add_line(line)

    if show_labels:
        text_kwargs = {"verticalalignment": "bottom"}
        text_kwargs.update({} if label_kwargs is None else label_kwargs)
        label_graph = (
            cast(nx.Graph, scdata.uns["flat_tree"])
            if graph_key == "epg" and "flat_tree" in scdata.uns
            else graph
        )
        label_positions = nx.get_node_attributes(label_graph, "pos")
        labels = nx.get_node_attributes(label_graph, label_key)

        for node in label_graph.nodes:
            if node not in label_positions:
                continue

            label = labels[node] if node in labels else node
            position = label_positions[node]

            if n_components == 2:
                cast(Any, ax).text(position[0], position[1], str(label), **text_kwargs)
            else:
                ax3d = cast(Any, ax)
                ax3d.text(
                    position[0],
                    position[1],
                    position[2],
                    str(label),
                    **text_kwargs,
                )

    return ax


@overload
def trajectory(
    scdata: ScData,
    obs: str,
    use_rep: str,
    graph_key: str = "epg",
    colors: Optional[Colors] = None,
    n_components: int = 2,
    title: Optional[Any] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    automatic_resize: bool = False,
    default_parameters: Optional[Any] = None,
    outfile: None = None,
    ax: Optional[Axes] = None,
    graph: Optional[Dict[str, Any]] = None,
    graph_z_offset: float = 0.0,
    label_key: str = "label",
    **kwargs: Any,
) -> Tuple[Figure, Axes]: ...


@overload
def trajectory(
    scdata: ScData,
    obs: str,
    use_rep: str,
    graph_key: str = "epg",
    colors: Optional[Colors] = None,
    n_components: int = 2,
    title: Optional[Any] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    automatic_resize: bool = False,
    default_parameters: Optional[Any] = None,
    *,
    outfile: Path,
    ax: Optional[Axes] = None,
    graph: Optional[Dict[str, Any]] = None,
    graph_z_offset: float = 0.0,
    label_key: str = "label",
    **kwargs: Any,
) -> None: ...


@overload
def trajectory(
    scdata: ScData,
    obs: str,
    use_rep: str,
    graph_key: str = "epg",
    colors: Optional[Colors] = None,
    n_components: int = 2,
    title: Optional[Any] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    automatic_resize: bool = False,
    default_parameters: Optional[Any] = None,
    *,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    graph: Optional[Dict[str, Any]] = None,
    graph_z_offset: float = 0.0,
    label_key: str = "label",
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]: ...


@anndata_or_mudata_checker
def trajectory(
    scdata: ScData,  # type: ignore
    obs: str,
    use_rep: str,
    graph_key: str = "epg",
    colors: Optional[Colors] = None,
    n_components: int = 2,
    title: Optional[Any] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    automatic_resize: bool = False,
    default_parameters: Optional[Any] = None,
    outfile: Optional[Path] = None,
    ax: Optional[Axes] = None,
    graph: Optional[Dict[str, Any]] = None,
    graph_z_offset: float = 0.0,
    label_key: str = "label",
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw an embedding and overlay a principal trajectory graph.

    The graph is read from `scdata.uns[graph_key]`, then drawn on the same axes
    as the cell embedding. Use `show_labels=True` to draw graph node labels.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Observation column in `scdata.obs` defining groups or values.
    use_rep: str
        Representation key in `scdata.obsm`.
    graph_key: str (default: "epg")
        Key pointing to the principal graph stored in `scdata.uns[graph_key]`.
    colors: matplotlib.colors.Colormap, optional
        Colormap or colors used to draw observations.
    n_components: 2 or 3 (default: 2)
        Number of plotted dimensions.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    show_legend: bool (default: True)
        Draw the legend when `scdata.obs[obs]` contains discrete values.
    show_labels: bool (default: False)
        Draw graph node labels.
    automatic_resize: bool (default: False)
        Resize figure to accommodate large legends.
    default_parameters: Callable, optional
        Function specifying default figure parameters.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    ax: Axes, optional
        Existing axes used for drawing.
    graph: dict, optional
        Keyword arguments passed to graph line drawing.
    graph_z_offset: float (default: 0.0)
        Vertical offset applied to graph lines in 3D.
    label_key: str (default: "label")
        Node attribute used as label when `show_labels=True`.
    **kwargs: Any
        Keyword arguments passed to `embedding`.

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.

    Raises
    ------
    KeyError
        If `graph_key` is missing from `scdata.uns`.

    References
    ----------
    Chen et al. (2019). Single-cell trajectories reconstruction, exploration
    and mapping of omics data with STREAM. Nature Communications, 10(1), 1903.
    """

    result = embedding(
        scdata,
        obs=obs,
        use_rep=use_rep,
        colors=colors,
        n_components=n_components,
        title=title,
        show_legend=show_legend,
        automatic_resize=automatic_resize,
        default_parameters=default_parameters,
        outfile=None,
        ax=ax,
        **kwargs,
    )
    fig, drawn_ax = cast(Tuple[Figure, Axes], result)

    graph_overlay(
        scdata,
        graph_key=graph_key,
        ax=drawn_ax,
        n_components=n_components,
        show_labels=show_labels,
        z_offset=graph_z_offset,
        label_key=label_key,
        label_kwargs=cast(Optional[Dict[str, Any]], kwargs.get("text")),
        **({} if graph is None else graph),
    )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        return None

    return fig, drawn_ax


@anndata_checker
def paga(
    adata: AnnData,
    obs: str,
    use_rep: str,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
    ax: Optional[Axes] = None,
    with_labels: bool = False,
    node_color: Any = black,
    outfile: Optional[Path] = None,
    **kwargs: Any,
) -> Union[Axes, None]:
    """
    Draw the partition-based graph abstraction (PAGA) graph with Matplotlib.

    To compute the PAGA matrix, run `scanpy.tl.paga`. PAGA can also be
    computed using scVelo [1]. In this case, use
    `adata.uns[edges] = adata.uns["paga"][edges]` before calling this function.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    obs: str
        Observation column in `adata.obs` defining groups.
    use_rep: str
        Representation key in `adata.obsm`.
    edges: str (default: "transitions_confidence")
        Key in `adata.uns` containing the adjacency matrix.
    threshold: float (default: 0.01)
        Confidence threshold used to prune adjacency values.
    ax: Axes (optional, default: None)
        Draw the paga graph in the specified Matplotlib axes.
    with_labels: bool (default: False)
        Add labels on the nodes.
    node_color: Union[Sequence[RGB], cycle, Mapping] (optional, default: None)
        Color specification for graph nodes.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Keyword arguments passed to `networkx.draw_networkx`.

    Returns
    -------
    Axes or None
        Matplotlib axes if `outfile` is None; otherwise None after saving the
        figure.

    References
    ----------
    Bergen et al. (2020). Generalizing RNA velocity to transient cell states
    through dynamical modeling. Nature Biotechnology, 38(12), 1408-1414.
    """

    if ax is None:
        ax = plt.gca()

    paga_graph = extract_paga_graph(
        adata=adata, obs=obs, use_rep=use_rep, edges=edges, threshold=threshold
    )
    barycenters = {
        node: position[:2]
        for node, position in nx.get_node_attributes(paga_graph, "pos").items()
    }

    if isinstance(node_color, MappingABC):
        node_color = [node_color[node] for node in paga_graph]
    elif (
        isinstance(node_color, SequenceABC)
        and not isinstance(node_color, str)
        and len(node_color) in [3, 4]
        and all(isinstance(channel, Number) for channel in node_color)
    ):
        node_color = [node_color for _ in paga_graph]

    nx.draw_networkx(
        paga_graph,
        pos=barycenters,
        ax=ax,
        with_labels=with_labels,
        node_color=cast(Any, node_color),
        **kwargs,
    )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None

    return ax


def add_graph(*args: Any, **kwargs: Any) -> Axes:
    """
    Deprecated alias for `graph_overlay`.
    """

    _warn_deprecated(
        "`bt.sct.pl.add_graph`",
        replacement="`bt.sct.pl.graph_overlay`",
        stacklevel=2,
    )
    return graph_overlay(*args, **kwargs)


def draw_paga(*args: Any, **kwargs: Any) -> Union[Axes, None]:
    """
    Deprecated alias for `paga()`.

    Use `paga()` instead. This alias is kept temporarily for backward
    compatibility.
    """

    _warn_deprecated("`draw_paga()`", replacement="`paga()`", stacklevel=2)

    return paga(*args, **kwargs)
