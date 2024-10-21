#!/usr/bin/env python

from collections.abc import Mapping
import typing

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes._axes import Axes
from itertools import cycle

from . import _colors

import networkx as nx

from ..tools import get_paga_graph
from .._check_anndata import _adata_arg_checking

@_adata_arg_checking
def draw_paga(
    adata,
    obs,
    obsm,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
    ax: Axes = None,
    with_labels: bool = False,
    node_color: typing.Optional[typing.Union[typing.Sequence[typing.Sequence[str]], cycle, Mapping]] = _colors.black,
    outfile: typing.Optional[Path] = None,
    **kwargs
):
    """
    Draw the paga graph with Matplotlib.
    To compute the PAGA matrix, please use `scanpy.tl.paga`.
    PAGA can also be computed using `scvelo`. In this case, please run previously
    adata.uns[`edges`] = adata.uns["paga"][`edges`]

    Parameters
    ----------
    adata
        Annotated data matrix
    obs
        The classification is retrieved by .obs[`obs`], which must be categorical/qualitative values
    obsm
        The data points are retrieved by the first columns in .obsm[`obsm`]
    edges
        The adjacency matrix-based data structure is retrieved by .uns[`edges`] (default: transitions_confidence)
        Please refer to 
    threshold
        Confidence threshold (default: 0.01)
    ax
        Draw the paga graph in the specified Matplotlib axes
    with_labels
        Set to True to draw labels on the nodes
    node_color
        Node color, array of node colors or mapping of node colors
    outfile
        If specified, save the figure
    **kwargs
        keyword arguments passed to networkx.draw_networkx() function

    Returns
    -------
    Depending on `outfile`, save figure or return axe.
    """
    
    if ax is None:
        ax = plt.gca()
    
    paga = get_paga_graph(
        adata=adata,
        obs=obs,
        obsm=obsm,
        edges=edges,
        threshold=threshold
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
        **kwargs
    )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return ax
