#!/usr/bin/env python

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
from anndata import AnnData

from ..._validation import _as_non_negative_number
from ..._warnings import _warn_deprecated
from .._typing import anndata_checker


@anndata_checker
def extract_paga_graph(
    adata: AnnData,
    obs: str,
    use_rep: str,
    edges: str = "transitions_confidence",
    threshold: float = 0.01,
) -> nx.DiGraph[Any]:
    """
    Extract the partition-based graph abstraction (PAGA) graph stored in `adata`.

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

    Returns
    -------
    nx.DiGraph
        PAGA graph with cluster barycenters stored as node positions.

    Raises
    ------
    KeyError
        If `edges` is not found in `adata.uns` and no fallback
        connectivities matrix is available.
    ValueError
        If `threshold` is negative.

    References
    ----------
    Bergen et al. (2020). Generalizing RNA velocity to transient cell states
    through dynamical modeling. Nature Biotechnology, 38(12), 1408-1414.
    """

    clusters = adata.obs[obs].cat.categories
    barycenters = {
        cluster: np.nanmean(adata.obsm[use_rep][adata.obs[obs] == cluster], axis=0)
        for cluster in clusters
    }

    if edges not in adata.uns:
        if "connectivities" in adata.uns:
            edges = "connectivities"
        else:
            raise KeyError(f"key '{edges}' not found in adata.uns")
    adjacency = adata.uns[edges].copy()

    threshold = _as_non_negative_number(threshold, "threshold")
    if threshold > 0:
        adjacency.data[adjacency.data < threshold] = 0
        adjacency.eliminate_zeros()
        adjacency = adjacency.todense()

    paga = nx.DiGraph()
    for node in clusters:
        paga.add_node(node, pos=barycenters[node])
    for i in range(len(adjacency)):
        for j in range(i + 1, len(adjacency)):
            if adjacency[i, j] > 0:
                paga.add_edge(clusters[j], clusters[i])
            if adjacency[j, i] > 0:
                paga.add_edge(clusters[i], clusters[j])

    return paga


def get_paga_graph(*args: Any, **kwargs: Any) -> nx.DiGraph[Any]:
    """
    Deprecated alias for `extract_paga_graph`.
    """

    _warn_deprecated(
        "`bt.sct.tl.get_paga_graph`",
        replacement="`bt.sct.tl.extract_paga_graph`",
        stacklevel=2,
    )
    return extract_paga_graph(*args, **kwargs)
