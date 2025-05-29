#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Mapping,
    Literal,
    Any,
    get_args
)
from .._typing import (
    ScData,
    Metric,
    anndata_or_mudata_checker
)

import warnings
import importlib

import math
import numpy as np
import pandas as pd
import anndata as ad

import networkx as nx
from networkx import Graph

from scipy.sparse import csr_matrix, issparse, diags
from sklearn.metrics import pairwise_distances

from ._maths import barycenters
from ._utils import choose_representation

@anndata_or_mudata_checker
def kneighbors_graph(
    scdata: ScData, # type: ignore
    n_neighbors: int,
    n_components: Optional[int] = None,
    use_rep: Optional[str] = None,
    metric: Metric = "euclidean",
    create_using: Optional[Graph] = nx.DiGraph,
    edge_attr: Optional[str] = "distance",
    index_or_name: Literal["index","name"] = "index",
    n_jobs: int = 1,
    **metric_kwds: Mapping[str, Any]
) -> Graph:
    """
    Compute the k-nearest neighbors-based graph using an embedding space.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    n_neighbors
        number of closest neighbors
    n_components
        Number of principal components or dimensions in the embedding space
        taken into account for each observation
    use_rep
        Use the indicated representation in adata.obsm
    metric
        Metric used when calculating pairwise distances between observations
    create_using
        Graph type to return
        If graph instance, then cleared it before computing k-nearest neighbors
    edge_attr
        Attribute to which the distance values are assigned on each edge
    index_or_name
        node names are referring either to index number ('index') or index name ('name')
    n_jobs
        Number of allocated processors
    **metric_kwds
        Any further parameters passed to the distance function

    Returns
    -------
    return Graph object.
    """

    from sklearn import neighbors
    
    X = choose_representation(
        scdata,
        use_rep=use_rep,
        n_components=n_components
    )
    weighted_adjacency_matrix = neighbors.kneighbors_graph(
        X=X,
        n_neighbors=n_neighbors,
        mode="distance",
        metric=metric,
        n_jobs=n_jobs,
        **metric_kwds
    )
    kneighbors_graph = nx.from_numpy_array(
        weighted_adjacency_matrix,
        create_using=create_using,
        edge_attr=edge_attr
    )

    if index_or_name == "index":
        return kneighbors_graph
    elif index_or_name == "name":
        return nx.relabel_nodes(kneighbors_graph,dict(zip(list(range(len(scdata.obs.index))), scdata.obs.index)))
    else:
        raise ValueError(f"invalid parameter value for 'index_or_name': expected 'index' or 'name' but received '{index_or_name}'")

class Knnbs(object):

    def __init__(
        self,
        n_neighbors: int,
        use_rep: Optional[str]="X_pca",
        n_components: Optional[str]=None,
        metric: Metric="euclidean",
        **metric_kwds: Mapping[str, Any]
    ):
    
        if isinstance(n_neighbors, int):
            if n_neighbors > 0:
                self.n_neighbors = n_neighbors
            else:
                raise ValueError(f"invalid argument value for 'n_neighbors': expected non-null positive but received '{n_neighbors}'")
        elif isinstance(n_neighbors, float):
            if n_neighbors.is_integer():
                self.n_neighbors = n_neighbors
            else:
                raise ValueError(f"invalid parameter value for 'n_neighbors': expected integer but received '{n_neighbors}'")
        else:
            raise TypeError(f"unsupported parameter type for 'n_neighbors': expected '{int}' but received '{type(n_neighbors)}'")

        if n_components is None:
            self.n_components = None
        elif isinstance(n_components, int):
            if n_components > 0:
                self.n_components = n_components
            else:
                raise ValueError(f"invalid argument value for 'n_components': expected non-null positive but received '{n_components}'")
        elif isinstance(n_components, float):
            if n_components.is_integer():
                self.n_components = n_components
            else:
                raise ValueError(f"invalid parameter value for 'n_components': expected integer but received '{n_components}'")
        else:
            raise TypeError(f"unsupported parameter type for 'n_components': expected '{int}' but received '{type(n_components)}'")
        
        if isinstance(use_rep, str):
            self.use_rep = use_rep
        else:
            raise TypeError(f"unsupported parameter type for 'use_rep': expected '{str}' but received '{type(use_rep)}'")
        
        if metric in get_args(Metric):
            self.metric = metric
        else:
            raise ValueError(f"unsupported parameter value for 'metric': '{metric}'")
        
        self.metric_kwds = metric_kwds
    
    def __repr__(self) -> str:
    
        return f"{self.__class__.__name__}(n_neighbors={self.n_neighbors},use_rep={self.use_rep},n_components={self.n_components},metric={self.metric},metric_kwds={self.metric_kwds})"

    def get(self, attribute: str) -> Any:
    
        if hasattr(self, attribute):
            return getattr(self, attribute)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{attribute}'")

    def fit(
        self,
        adata: ad.AnnData,
        obs: str,
        n_jobs: int=1
    ) -> None:
        
        from sklearn.metrics import pairwise_distances

        X = choose_representation(adata, use_rep=self.use_rep, n_components=self.n_components)

        _kneighbors_graph = kneighbors_graph(
            adata,
            n_neighbors=self.n_neighbors,
            n_components=self.n_components,
            use_rep=self.use_rep,
            metric=self.metric,
            create_using=nx.Graph,
            index_or_name="name",
            n_jobs=n_jobs
        )

        _barycenters = barycenters(
            adata,
            obs=obs,
            use_rep=self.use_rep,
            n_components=self.n_components
        )
        for key, value in _barycenters.items():
            barycenter_coordinate = value.reshape(1,-1)
            distances = pairwise_distances(
                X,
                barycenter_coordinate,
                metric=self.metric,
                n_jobs=n_jobs
            ).reshape(1,-1)
            _kneighbors_graph.add_node(key)
            knn_indices = list(np.argpartition(distances, kth=self.n_neighbors, axis=1)[:, :self.n_neighbors].reshape(-1))
            distances = list(distances[0,knn_indices].reshape(-1))
            for obs_name, distance in zip(adata.obs.index[knn_indices], distances):
                _kneighbors_graph.add_edge(key, obs_name, distance=distance)
        
        if nx.number_connected_components(_kneighbors_graph) > 1:
            raise nx.NetworkXAlgorithmError(f"'kneighbors_graph' not weakly connected")
        else:
            self.kneighbors_graph = _kneighbors_graph
            self.obs = adata.obs[obs]
            return None

    def shortest_path_lengths(
        self,
        method: str="dijkstra",
        n_jobs: int=1
    ) -> None:
            
        def shortest_path_lengths_from(
            source: str,
        ):
            distances = pd.Series(data=self.obs.index,index=self.obs.index)
            distances = distances.apply(
                lambda x: nx.shortest_path_lengths(
                    self.kneighbors_graph,
                    source=source,
                    target=x,
                    weight="distance",
                    method=method
                )
            )
            return (source, distances)

        try:
            _multiprocess_is_available = importlib.util.find_spec("multiprocess") is not None
        except:
            _multiprocess_is_available = importlib.find_loader("multiprocess") is not None

        if _multiprocess_is_available:
            from multiprocess import Pool
            with Pool(processes=n_jobs) as pool:
                shortest_path_lengths_ls = pool.map(shortest_path_lengths_from, self.obs.cat.categories)
        else:
            warnings.warn("module multiprocess not available: not using CPU parallelization")
            shortest_path_lengths_ls = []
            for obs in self.obs.cat.categories:
                shortest_path_lengths_ls.append((obs,shortest_path_lengths_from(obs)))
        
        shortest_path_lengths_dict = {}
        for k, v in shortest_path_lengths_ls:
            shortest_path_lengths_dict[k] = v
        self.shortest_path_lengths_df = pd.DataFrame.from_dict(data=shortest_path_lengths_dict, orient="columns")
    
    def subclusters(
            self,
            key: str="knnbs",
            size: int=30
    ) -> pd.Series:
        # TODO
        return

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
