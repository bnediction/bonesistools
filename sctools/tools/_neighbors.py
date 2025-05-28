#!/usr/bin/env python

from typing import Optional, Union
from .._typing import anndata_or_mudata_checker, ScData

import numpy as np
import math
from scipy.sparse import csr_matrix, issparse, diags
from sklearn.metrics import pairwise_distances

import anndata as ad

@anndata_or_mudata_checker
def _shared_nearest_neighbors_graph(
    scdata: ScData, # type: ignore
    cluster_key: str,
    prune_snn: float
) -> csr_matrix:

    k_neighbors = scdata.uns[cluster_key]["params"]["n_neighbors"] - 1
    if prune_snn < 0:
        raise ValueError("'prune_snn' parameter must be positive, aborting")
    elif prune_snn < 1:
        prune_snn = math.ceil(k_neighbors*prune_snn)
    elif prune_snn >= k_neighbors:
        raise ValueError("'prune_snn' parameter must be smaller than 'n_neighbors' used for KNN computation, aborting")

    n_cells = scdata.n_obs
    distances_key = scdata.uns[cluster_key]["distances_key"]

    neighborhood_graph = scdata.obsp[distances_key].copy()
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

@anndata_or_mudata_checker
def shared_neighbors(
    scdata: ScData, # type: ignore
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[float] = 1/15,
    metric: Optional[str] = "euclidean",
    normalize_connectivities: bool = True,
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False
) -> Union[ad.AnnData, None]:
    """Compute a shared neighborhood (SNN) graph of observations.

    The neighbor search relies on a previously computed neighborhood graph
    (such as kNN algorithm).

    Parameters
    ----------
    scdata
        Annotated data matrix
    knn_key
        If not specified, the used neighbors data are retrieved from .uns['neighbors'],
        otherwise the used neighbors data are retrieved from .uns['key_added']
    snn_key
        If not specified, the shared neighbors data are stored in .uns['shared_neighbors']
        If specified, the shared neighbors data are added to .uns['key_added']
    prune_snn
        If zero value, no prunning is performed. If strictly positive, removes edge between two neighbors
        in the shared neighborhood graph who have a number of neighbors less than the specified value
        Value can be relative (float between 0 and 1) or absolute (integer between 1 and k)
    metric
        Metric used for computing distances between two neighbors by using scdata.obsm
    normalize_connectivities
        If false, connectivities provide the absolute number of shared neighbors (integer between 0 and k),
        otherwise provide the relative number of shared neighbors (float between 0 and k)
    distances_key
        If specified, distances are stored in .obsp[distances_key],
        otherwise in .obsp[snn_key+'_distances']
    connectivities_key
        If specified, distances are stored in .obsp[connectivities_key],
        otherwise in .obsp[snn_key+'_connectivities']
    copy
        Return a copy instead of writing to scdata.

    Returns
    -------
    Depending on 'copy', updates or returns 'scdata' with the following:

    See 'snn_key' parameter description for the storage path of
    connectivities and distances.

    **connectivities** : sparse matrix.
        Weighted adjacency matrix of the shared neighborhood graph.
        Weights should be interpreted as number of shared neighbors.
    **distances** : sparse matrix of dtype 'float64'.
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """

    if knn_key not in scdata.uns:
        raise KeyError((
            "Neighborhood graph not already computed or not finding. "
            "Please use 'scanpy.pp.neighbors' function before or "
            "specify 'key_added' parameter when scanpy.pp.neighbors has been called, aborting"
    ))
    if prune_snn is None:
        prune_snn = 0
    if metric is None:
        metric = "euclidean"
    if distances_key is None:
        distances_key = f"{snn_key}_distances"
    if connectivities_key is None:
        connectivities_key = f"{snn_key}_connectivities"
    n_neighbors = scdata.uns[knn_key]["params"]["n_neighbors"]

    scdata = scdata.copy() if copy else scdata
    
    snn_graph = _shared_nearest_neighbors_graph(scdata, cluster_key=knn_key, prune_snn = prune_snn)

    n_pcs = scdata.uns[knn_key]["params"]["n_pcs"]
    obsm = scdata.uns[knn_key]["params"]["use_rep"]

    X = scdata.obsm[obsm][:,0:n_pcs]
    zeros_ones = snn_graph.toarray()
    zeros_ones[zeros_ones > 0] = 1
    
    distances_matrix = pairwise_distances(X, metric=metric)
    distances_matrix = np.multiply(zeros_ones, distances_matrix)
    distances_matrix = csr_matrix(distances_matrix)
    connectivities_matrix = snn_graph.copy()
    if normalize_connectivities:
        connectivities_matrix = connectivities_matrix.astype(float)
        connectivities_matrix.data /= n_neighbors

    scdata.obsp[distances_key] = distances_matrix
    scdata.obsp[connectivities_key] = connectivities_matrix

    scdata.uns[snn_key] = dict()
    scdata.uns[snn_key]["distances_key"] = distances_key
    scdata.uns[snn_key]["connectivities_key"] = connectivities_key
    scdata.uns[snn_key]["params"] = {
        "knn_base": f"scdata.uns['{knn_key}']",
        "prune_snn": prune_snn if prune_snn >= 1 else math.ceil(n_neighbors*prune_snn),
        "metric": metric
    }

    return scdata if copy else None
