#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping, Sequence
from numbers import Number
from pathlib import Path
from typing import (
    Any,
    Optional,
    Union,
    cast,
)

import matplotlib.pyplot as plt
import networkx as nx
from anndata import AnnData
from matplotlib.axes._axes import Axes

from .._typing import anndata_checker
from ..tools import get_paga_graph
from ._colors import black


@anndata_checker
def draw_paga(
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
    [1] Bergen et al. (2020). Generalizing RNA velocity to transient cell
    states through dynamical modeling. Nature biotechnology, 38(12), 1408-1414
    (https://doi.org/10.1038/s41587-020-0591-3)
    """

    if ax is None:
        ax = plt.gca()

    paga = get_paga_graph(
        adata=adata, obs=obs, use_rep=use_rep, edges=edges, threshold=threshold
    )
    barycenters = {
        node: position[:2]
        for node, position in nx.get_node_attributes(paga, "pos").items()
    }

    if isinstance(node_color, Mapping):
        node_color = [node_color[node] for node in paga]
    elif (
        isinstance(node_color, Sequence)
        and not isinstance(node_color, str)
        and len(node_color) in [3, 4]
        and all(isinstance(channel, Number) for channel in node_color)
    ):
        node_color = [node_color for _ in paga]

    nx.draw_networkx(
        paga,
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
    else:
        return ax
