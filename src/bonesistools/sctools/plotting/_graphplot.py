#!/usr/bin/env python

from typing import Optional, Union, Sequence, Mapping
from pathlib import Path
from anndata import AnnData
from ._typing import RGB
from .._typing import anndata_checker

import matplotlib.pyplot as plt
from matplotlib.axes._axes import Axes
from itertools import cycle

from . import _colors

import networkx as nx

from ..tools import get_paga_graph


@anndata_checker
def draw_paga(
    adata: AnnData,
    obs: str,
    use_rep: str,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
    ax: Optional[Axes] = None,
    with_labels: bool = False,
    node_color: Optional[Union[Sequence[RGB], cycle, Mapping]] = _colors.black,
    outfile: Optional[Path] = None,
    **kwargs,
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
    **kwargs: Mapping[str, Any]
        Keyword arguments passed to `networkx.draw_networkx`.

    Returns
    -------
    Axes or None
        Matplotlib axes if `outfile` is None; otherwise None after saving the
        figure.

    References
    ----------
    [1] Bergen et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling.
    Nature biotechnology, 38(12), 1408-1414 (https://doi.org/10.1038/s41587-020-0591-3)
    """

    if ax is None:
        ax = plt.gca()

    paga = get_paga_graph(
        adata=adata, obs=obs, use_rep=use_rep, edges=edges, threshold=threshold
    )
    barycenters = nx.get_node_attributes(paga, "pos")

    if isinstance(node_color, Mapping):
        node_color = [node_color[node] for node in paga]

    nx.draw_networkx(
        paga,
        pos=barycenters,
        ax=ax,
        with_labels=with_labels,
        node_color=node_color,
        **kwargs,
    )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return ax
