#!/usr/bin/env python

from typing import Optional, Union

import numpy as np
import math
from scipy.sparse import csr_matrix, issparse, diags
from sklearn.metrics import pairwise_distances

import anndata as ad

from .._check_anndata import _adata_arg_checking

@_adata_arg_checking
def _shared_nearest_neighbors_graph(
    adata: ad.AnnData,
    cluster_key: str,
    prune_snn: float
) -> csr_matrix:

    k_neighbors = adata.uns[cluster_key]["params"]["n_neighbors"] - 1
    if prune_snn < 0:
        raise ValueError("`prune_snn` parameter must be positive, aborting")
    elif prune_snn < 1:
        prune_snn = math.ceil(k_neighbors*prune_snn)
    elif prune_snn >= k_neighbors:
        raise ValueError("`prune_snn` parameter must be smaller than `n_neighbors` used for KNN computation, aborting")

    n_cells = adata.n_obs
    distances_key = adata.uns[cluster_key]["distances_key"]

    neighborhood_graph = adata.obsp[distances_key].copy()
    if not issparse(neighborhood_graph):
        neighborhood_graph = csr_matrix(neighborhood_graph)
    neighborhood_graph.data[neighborhood_graph.data > 0] = 1

    neighborhood_graph = neighborhood_graph * neighborhood_graph.transpose()
    neighborhood_graph -= (k_neighbors * diags(np.ones(n_cells), offsets=0, shape=(n_cells, n_cells)))
    neighborhood_graph.sort_indices()
    neighborhood_graph = neighborhood_graph.astype(dtype=np.int8)

    if prune_snn:
        mask = (neighborhood_graph.data <= prune_snn)
        neighborhood_graph.data[mask] = 0
        neighborhood_graph.eliminate_zeros()

    return neighborhood_graph

@_adata_arg_checking
def shared_neighbors(
    adata: ad.AnnData,
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[float] = 1/15,
    metric: Optional[str] = "euclidean",
    normalize_similarities: bool = True,
    distances_key: Optional[str] = None,
    similarities_key: Optional[str] = None,
    copy: bool = False
) -> Union[ad.AnnData, None]:
    """Compute a shared neighborhood (SNN) graph of observations.

    The neighbor search relies on a previously computed neighborhood graph
    (such as kNN algorithm).

    Parameters
    ----------
    adata
        Annotated data matrix.
    knn_key
        If not specified, the used neighbors data are retrieved from .uns['neighbors'].
        If specified, the used neighbors data are retrieved from .uns[key_added].
    snn_key
        If not specified, the shared neighbors data are stored in .uns['shared_neighbors'].
        If specified, the shared neighbors data are added to .uns[key_added].
    prune_snn
        If zero value, no prunning is performed. If strictly positive, removes edge between two neighbors
        in the shared neighborhood graph who have a number of neighbors less than the specified value.
        Value can be relative (float between 0 and 1) or absolute (integer between 1 and k).
    metric
        Metric used for computing distances between two neighbors by using .obsm['X_pca'].
    normalize_similarities
        If false, similarities provide the absolute number of shared neighbors (integer between 0 and k),
        otherwise provide the relative number of shared neighbors (float between 0 and k).
    distances_key
        If specified, distances are stored in .obsp[distances_key],
        otherwise in .obsp[snn_key+'_distances'].
    similarities_key
        If specified, distances are stored in .obsp[similarities_key],
        otherwise in .obsp[snn_key+'_similarities'].
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:

    See `snn_key` parameter description for the storage path of
    similarities and distances.

    **similarities** : sparse matrix.
        Weighted adjacency matrix of the shared neighborhood graph.
        Weights should be interpreted as number of shared neighbors.
    **distances** : sparse matrix of dtype `float64`.
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """

    if knn_key not in adata.uns:
        raise KeyError((
            "Neighborhood graph not already computed or not finding. "
            "Please use `scanpy.pp.neighbors` function before or "
            "specify `key_added` parameter when scanpy.pp.neighbors has been called, aborting"
    ))
    if prune_snn is None:
        prune_snn = 0
    if metric is None:
        metric = "euclidean"
    if distances_key is None:
        distances_key = f"{snn_key}_distances"
    if similarities_key is None:
        similarities_key = f"{snn_key}_similarities"
    n_neighbors = adata.uns[knn_key]["params"]["n_neighbors"]

    adata = adata.copy() if copy else adata
    
    snn_graph = _shared_nearest_neighbors_graph(adata, cluster_key=knn_key, prune_snn = prune_snn)

    n_pcs = adata.uns[knn_key]["params"]["n_pcs"]

    X = adata.obsm["X_pca"][:,0:n_pcs]
    zeros_ones = snn_graph.toarray()
    zeros_ones[zeros_ones > 0] = 1
    
    distances_matrix = pairwise_distances(X, metric=metric)
    distances_matrix = np.multiply(zeros_ones, distances_matrix)
    distances_matrix = csr_matrix(distances_matrix)
    similarities_matrix = snn_graph.copy()
    if normalize_similarities:
        similarities_matrix = similarities_matrix.astype(float)
        similarities_matrix.data /= n_neighbors

    adata.obsp[distances_key] = distances_matrix
    adata.obsp[similarities_key] = similarities_matrix

    adata.uns[snn_key] = dict()
    adata.uns[snn_key]["distances_key"] = distances_key
    adata.uns[snn_key]["similarities_key"] = similarities_key
    adata.uns[snn_key]["params"] = {
        "knn_base": f"adata.uns['{knn_key}']",
        "prune_snn": prune_snn if prune_snn >= 1 else math.ceil(n_neighbors*prune_snn),
        "metric": metric
    }

    return adata if copy else None
