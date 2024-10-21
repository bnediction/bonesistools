#!/usr/bin/env python

from anndata import AnnData

import numpy as np

import networkx as nx

from .. import _adata_arg_checking

@_adata_arg_checking
def get_paga_graph(
    adata: AnnData,
    obs: str,
    obsm: str,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
):
    """
    Get the partition-based graph abstraction (PAGA) graph stored in `adata`
    with a graph-based data structure.
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
    threshold
        confidence threshold (default: 0.01)

    Returns
    -------
    nx.DiGraph (directed graph-based data structure)
    """

    clusters = adata.obs[obs].cat.categories
    barycenters = {cluster: np.nanmean(adata.obsm[obsm][adata.obs[obs] == cluster], axis=0) for cluster in clusters}

    if edges not in adata.uns:
        if "connectivities" in adata.uns:
            edges = "connectivities"
        else:
            raise ValueError(f"edges `{edges}` not in adata.uns.")
    adjacency = adata.uns[edges].copy()

    if not isinstance(threshold, float):
        raise TypeError(f"`threshold` is not a float.")
    elif threshold < 0:
        raise ValueError(f"`threshold` is not positive.")
    elif threshold > 0:
        adjacency.data[adjacency.data < threshold] = 0
        adjacency.eliminate_zeros()
        adjacency = adjacency.todense()

    paga = nx.DiGraph()
    for node in clusters:
        paga.add_node(node, pos=barycenters[node])
    for i in range(len(adjacency)):
        for j in range(i+1, len(adjacency)):
            if adjacency[i,j] > 0:
                paga.add_edge(clusters[j], clusters[i])
            if adjacency[j, i] > 0:
                paga.add_edge(clusters[i], clusters[j])

    return paga
