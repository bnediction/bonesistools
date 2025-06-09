#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Sequence,
)
from collections.abc import Mapping
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
    **kwargs
) -> Union[Axes, None]:
    """
    Draw the paga graph with Matplotlib. To compute the paga matrix, please run 'scanpy.tl.paga'.
    Paga graph can also be computed using scvelo [1]. In this case, please use
    adata.uns['edges'] = adata.uns["paga"]['edges'] before.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    obs: str
        The classification is retrieved by adata.obs['obs'], which must be categorical/qualitative values.
    use_rep: str
        The data points are retrieved by the first columns in adata.obsm['use_rep'].
    edges: str (default: 'transitions_confidence')
        The adjacency matrix-based data structure is retrieved by adata.uns['edges'].
    threshold: float (default: 0.01)
        Confidence threshold used.
    ax: Axes (optional, default: None)
        Draw the paga graph in the specified Matplotlib axes.
    with_labels: bool (default: False)
        Add labels on the nodes.
    node_color: Union[Sequence[RGB], cycle, Mapping] (optional, default: None)
        Specify color for each nodes.
    outfile: Path (optional, default: None)
        If specified, save the figure.
    **kwargs
        keyword arguments passed to function networkx.draw_networkx.

    Returns
    -------
    Depending on 'outfile', save figure or return axe.

    References
    ----------
    [1] Bergen et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling.
    Nature biotechnology, 38(12), 1408-1414 (https://doi.org/10.1038/s41587-020-0591-3)
    """
    
    if ax is None:
        ax = plt.gca()
    
    paga = get_paga_graph(
        adata=adata,
        obs=obs,
        use_rep=use_rep,
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
