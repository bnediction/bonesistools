#!/usr/bin/env python

from anndata import AnnData
from .._typing import anndata_checker

import numpy as np

import networkx as nx

@anndata_checker
def get_paga_graph(
    adata: AnnData,
    obs: str,
    use_rep: str,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
) -> nx.DiGraph:
    """
    Get the partition-based graph abstraction (paga) graph stored in 'adata' with a graph-based data structure.
    To compute the paga matrix, please run 'scanpy.tl.paga'. Paga can also be computed using scvelo [1].
    In this case, please use adata.uns['edges'] = adata.uns["paga"]['edges'] before.
    
    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    obs: str
        The classification is retrieved by adata.obs['obs'], which must be categorical/qualitative values.
    use_rep: str
        The data points are retrieved by the first columns in .obsm['use_rep'].
    edges: str (default: transitions_confidence)
        The adjacency matrix-based data structure is retrieved by adata.uns['edges'].
    threshold: float (default: 0.01)
        Confidence threshold.

    Returns
    -------
    Return DiGraph corresponding to the paga.

    References
    ----------
    [1] Bergen et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling.
    Nature biotechnology, 38(12), 1408-1414 (https://doi.org/10.1038/s41587-020-0591-3)
    """

    clusters = adata.obs[obs].cat.categories
    barycenters = {cluster: np.nanmean(adata.obsm[use_rep][adata.obs[obs] == cluster], axis=0) for cluster in clusters}

    if edges not in adata.uns:
        if "connectivities" in adata.uns:
            edges = "connectivities"
        else:
            raise KeyError(f"key '{edges}' not found in adata.uns")
    adjacency = adata.uns[edges].copy()

    if not isinstance(threshold, float):
        raise TypeError(f"unsupported argument type for 'threshold': expected {float} but received {type(threshold)}")
    elif threshold < 0:
        raise ValueError(f"invalid argument value for 'threshold': expected positive value but received '{threshold}'")
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
