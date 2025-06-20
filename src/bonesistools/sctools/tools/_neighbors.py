#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Sequence,
    Any,
    Type
)
from collections.abc import Mapping
try:
    from typing import Literal, get_args
except ImportError:
    from typing_extensions import Literal, get_args
from .._typing import (
    AnnData,
    ScData,
    Metric,
    Shortest_Path_Method,
    anndata_or_mudata_checker
)

import warnings
import importlib

import math
from itertools import combinations

import numpy as np
import pandas as pd

import networkx as nx
from networkx import Graph, DiGraph

from scipy.sparse import csr_matrix, issparse, diags
from sklearn.metrics import pairwise_distances

from ._maths import barycenters
from ._utils import choose_representation

@anndata_or_mudata_checker
def kneighbors_graph(
    scdata: ScData, # type: ignore
    n_neighbors: int,
    use_rep: Optional[str] = None,
    n_components: Optional[int] = None,
    metric: Metric = "euclidean",
    create_using: Type = DiGraph,
    edge_attr: str = "distance",
    index_or_name: Literal["index","name"] = "index",
    n_jobs: int = 1,
    **metric_kwds: Mapping[str, Any]
) -> Graph:
    """
    Compute the k-nearest neighbors-based graph using an embedding space.

    Parameters
    ----------
    scdata: ad.AnnData | md.MuData
        Unimodal or multimodal annotated data matrix.
    n_neighbors: int
        number of closest neighbors.
    use_rep: str (optional, default: None)
        Use the indicated representation in scdata.obsm.
    n_components: int (optional, default: None)
        Number of principal components or dimensions in the embedding space
        taken into account for each observation.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    create_using: Type (default: nx.DiGraph)
        Graph type to return. If graph instance, then cleared it before computing k-nearest neighbors.
    edge_attr: str (default: 'distance')
        Attribute to which the distance values are assigned on each edge.
    index_or_name: 'index' | 'name' (default: 'index')
        node names are referring either to index number ('index') or index name ('name').
    n_jobs: int (default: 1)
        Number of allocated processors.
    **metric_kwds
        Any further parameters passed to the distance function.

    Returns
    -------
    Return Graph object.
    The graph stores weighted edges that connects two nodes.
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
        raise ValueError(f"invalid argument value for 'index_or_name': expected 'index' or 'name' but received '{index_or_name}'")


class Knnbs(object):
    """
    Class for using k-nearest neighbors-based subclusters (knnbs) algorithm.
    A k-nearest neighbors graph is constructed in order to compute distance
    between cells and clusters' barycenters. Two methods are used for finding subclusters:
    (1) searching for cell manifolds maximizing distances to other clusters' barycenters
    (2) searching for cell manifolds minimizing distances to self barycenter

    Parameters
    ----------
    n_neighbors: int
        Number of closest neighbors.
    use_rep: str (default: 'X_pca')
        Use the indicated representation in adata.obsm.
    n_components: int (optional, default: None)
        Number of principal components or dimensions in the embedding space
        taken into account for each observation.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    **metric_kwds
        Any further parameters passed to the distance function.
    """

    def __init__(
        self,
        n_neighbors: int,
        use_rep: str = "X_pca",
        n_components: Optional[int] = None,
        metric: Metric = "euclidean",
        **metric_kwds: Mapping[str, Any]
    ):
    
        if isinstance(n_neighbors, int):
            if n_neighbors > 0:
                self.n_neighbors = n_neighbors
            else:
                raise ValueError(f"invalid argument value for 'n_neighbors': expected non-null positive value but received '{n_neighbors}'")
        elif isinstance(n_neighbors, float):
            if n_neighbors.is_integer():
                self.n_neighbors = n_neighbors
            else:
                raise ValueError(f"invalid argument value for 'n_neighbors': expected integer but received '{n_neighbors}'")
        else:
            raise TypeError(f"unsupported argument type for 'n_neighbors': expected {int} but received {type(n_neighbors)}")

        if n_components is None:
            self.n_components = None
        elif isinstance(n_components, int):
            if n_components > 0:
                self.n_components = n_components
            else:
                raise ValueError(f"invalid argument value for 'n_components': expected non-null positive value but received '{n_components}'")
        elif isinstance(n_components, float):
            if n_components.is_integer():
                self.n_components = n_components
            else:
                raise ValueError(f"invalid argument value for 'n_components': expected integer but received '{n_components}'")
        else:
            raise TypeError(f"unsupported argument type for 'n_components': expected {int} but received {type(n_components)}")
        
        if isinstance(use_rep, str):
            self.use_rep = use_rep
        else:
            raise TypeError(f"unsupported argument type for 'use_rep': expected {str} but received {type(use_rep)}")
        
        if metric in get_args(Metric):
            self.metric = metric
        else:
            raise ValueError(f"invalid argument value for 'metric': {metric}")
        
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
        adata: AnnData,
        obs: str,
        n_jobs: int = 1
    ) -> None:
        """
        Compute the k-nearest neighbors-based graph using an embedding space.
        If parameters set at the knnbs instanciation are not desired,
        first update knnbs object-related attributes with desired options.

        Parameters
        ----------
        adata: ad.AnnData
            Unimodal annotated data matrix.
        obs: str
            Column name in 'adata.obs' used as reference for clusters.
        n_jobs: int (default: 1)
            Number of allocated processors.

        Returns
        -------
        Update Knnbs object.
        Add attribute 'shortest_path_lengths_df' to Knnbs instance.
        """

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
            warnings.warn(f"'kneighbors_graph' not weakly connected: add edges for joining connected components")
            scc = list(nx.connected_components(_kneighbors_graph))
            scc = [list(cc - set(_barycenters.keys())) for cc in scc]
            for paired_scc in combinations(scc, 2):
                X = choose_representation(adata[paired_scc[0],:], use_rep=self.use_rep)
                Y = choose_representation(adata[paired_scc[1],:], use_rep=self.use_rep)
                dists = pairwise_distances(
                    X,
                    Y,
                    metric=self.metric,
                    n_jobs=n_jobs
                )
                i, j  = np.unravel_index(np.argmin(dists), shape=dists.shape, order="C")
                _kneighbors_graph.add_edge(paired_scc[0][i], paired_scc[1][j], distance=dists[i,j])
        
        self.kneighbors_graph = _kneighbors_graph
        self.obs = adata.obs[obs]
        return None

    def shortest_path_lengths(
        self,
        method: Shortest_Path_Method = "dijkstra",
        n_jobs: int = 1
    ) -> None:
        """
        Compute pairwise shortest path lengths between cells and barycenters.

        Parameters
        ----------
        method: 'dijkstra' | 'bellman-ford' (default: 'dijkstra')
            Algorithm used to compute the shortest path lengths.
        n_jobs: int (default: 1)
            Number of allocated processors.
    
        Returns
        -------
        Update Knnbs object.
        Add attribute 'shortest_path_lengths_df' to Knnbs instance.
        """
        
        def shortest_path_lengths_from(
            source: str,
        ):
            distances = pd.Series(
                data=self.obs.index,
                index=self.obs.index
            )
            distances = distances.apply(
                lambda x: nx.shortest_path_length(
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
            for category in self.obs.cat.categories:
                shortest_path_lengths_ls.append(shortest_path_lengths_from(category))
        
        shortest_path_lengths_dict = {}
        for k, v in shortest_path_lengths_ls:
            shortest_path_lengths_dict[k] = v

        self.shortest_path_lengths_df = pd.DataFrame.from_dict(data=shortest_path_lengths_dict, orient="columns")
        return None
    
    def find_furthest_cells_to_other_barycenters(
            self,
            size: int = 30,
            key: str = "knnbs",
            clusters: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Find cluster related-cell manifolds maximizing distances to other clusters' barycenters.

        Parameters
        ----------
        size: int (default: 30)
            Number of cells in each macrostate.
            If 'size' is superior to cluster size, cluster related-cell manifolds are equal to its cluster.
        key: str (default: 'knnbs')
            Pandas serie name.
        clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed.
    
        Returns
        -------
        Return Series object.
        Series stores furthest cells to other barycenters.
        """
        
        if clusters is None:
            clusters = self.obs.cat.categories

        subclusters_series = pd.Series(
            data=np.nan,
            index=self.obs.index,
            dtype="category",
            copy=True
        ).cat.add_categories(self.obs.cat.categories)
        if key is not None:
            subclusters_series.name = key
        
        shortest_path_length_free_from_self_cluster_df = self.shortest_path_lengths_df.copy(deep=True)
        for idx, obs in self.obs.items():
            if not pd.isna(obs):
                shortest_path_length_free_from_self_cluster_df.at[idx,obs] = np.nan

        _min_dists = shortest_path_length_free_from_self_cluster_df.min(axis=1, skipna=True)
        _min_dists.name = "min_dist"
        min_dists = _min_dists.to_frame().merge(
            right=self.obs.to_frame(),
            left_index=True,
            right_index=True
        )
        for cluster in clusters:
            _min_dists_cluster = min_dists[min_dists[self.obs.name] == cluster]["min_dist"]
            if len(_min_dists_cluster) < size:
                _obs = _min_dists_cluster.index
            else:
                _idx = np.argpartition(_min_dists_cluster, kth=len(_min_dists_cluster) - size)[-size:]
                _obs = _min_dists_cluster.iloc[_idx].index
            subclusters_series.loc[_obs] = cluster
        
        return subclusters_series
    
    def find_closest_cells_to_self_barycenter(
        self,
        size: int = 30,
        key: str = "knnbs",
        clusters: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Find cluster related-cell manifolds minimizing distances to self barycenter.

        Parameters
        ----------
        size: int (default: 30)
            Number of cells in each macrostate.
            If 'size' is superior to cluster size, cluster related-cell manifolds are equal to its cluster.
        key: str (default: 'knnbs')
            Pandas serie name.
        clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed.
    
        Returns
        -------
        Return Series object.
        Series stores closest cells to self barycenter.
        """
        
        if clusters is None:
            clusters = self.obs.cat.categories

        subclusters_series = pd.Series(
            data=np.nan,
            index=self.obs.index,
            dtype="category",
            copy=True
        ).cat.add_categories(self.obs.cat.categories)
        if key is not None:
            subclusters_series.name = key

            _dists = pd.Series(
                data=np.nan,
                index=self.obs.index
            )
            _dists.name = "dist"
            for idx, obs in self.obs.items():
                if not pd.isna(obs):
                    _dists.at[idx] = self.shortest_path_lengths_df.loc[idx,obs]

        _dists = _dists.to_frame().merge(
            right=self.obs.to_frame(),
            left_index=True,
            right_index=True
        )
        for cluster in clusters:
            _dists_cluster = _dists[_dists[self.obs.name] == cluster]["dist"]
            if len(_dists_cluster) < size:
                _obs = _dists_cluster.index
            else:
                _idx = np.argpartition(_dists_cluster, kth=size)[:size]
                _obs = _dists_cluster.iloc[_idx].index
            subclusters_series.loc[_obs] = cluster
        
        return subclusters_series
    
    def knnbs(
        self,
        size: int = 30,
        key: str = "knnbs",
        subclusters_maximizing_distances: Optional[Sequence[str]] = None,
        subclusters_minimizing_distances: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Find cluster related-cell manifolds using k-nearest neighbors-based subclusters algorithm.

        Parameters
        ----------
        size
            Number of cells in each macrostate.
            If 'size' is superior to cluster size, cluster related-cell manifolds are equal to its cluster.
        key: str (default: 'knnbs')
            Pandas serie name.
        subclusters_maximizing_distances: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed
            by maximizing distances to other clusters' barycenters.
        subclusters_minimizing_distances: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed
            by minimizing distances to self barycenter .
    
        Returns
        -------
        Return Series object.
        Series stores subclusters derived from k-nearest neighbors-based subclusters algorithm.
        """
        
        if subclusters_maximizing_distances is None and subclusters_minimizing_distances is None:
            subclusters_maximizing_distances = set(self.obs.cat.categories)
            subclusters_minimizing_distances = set()
        else:
            subclusters_maximizing_distances = set(subclusters_maximizing_distances)
            subclusters_minimizing_distances = set(subclusters_minimizing_distances)
        
        if subclusters_maximizing_distances & subclusters_minimizing_distances:
            raise RuntimeError("'subclusters_maximizing_distances' and 'subclusters_minimizing_distances' are not disjoint")

        if subclusters_maximizing_distances:
            max_dists = self.find_furthest_cells_to_other_barycenters(
                size=size,
                key="max_dists",
                clusters=subclusters_maximizing_distances
            )
        else:
            max_dists = pd.Series(
                data=np.nan,
                index=self.obs.index,
                dtype="category",
                copy=True
            ).cat.add_categories(self.obs.cat.categories)
            max_dists.name="max_dists"

        if subclusters_minimizing_distances:
            min_dists = self.find_closest_cells_to_self_barycenter(
                size=size,
                key="min_dists",
                clusters=subclusters_minimizing_distances
            )
        else:
            min_dists = pd.Series(
                data=np.nan,
                index=self.obs.index,
                dtype="category",
                copy=True
            ).cat.add_categories(self.obs.cat.categories)
            min_dists.name="min_dists"
        
        dists = max_dists.to_frame().merge(
            right=min_dists.to_frame(),
            left_index=True,
            right_index=True
        )
        subclusters = dists.bfill(axis=1).iloc[:, 0]
        subclusters.name = key

        return subclusters

@anndata_or_mudata_checker
def _shared_nearest_neighbors_graph(
    scdata: ScData, # type: ignore
    cluster_key: str,
    prune_snn: float
) -> csr_matrix:

    n_neighbors = scdata.uns[cluster_key]["params"]["n_neighbors"] - 1
    if prune_snn < 0:
        raise ValueError(f"invalid argument value for 'prune_snn': expected positive value but received '{prune_snn}'")
    elif prune_snn < 1:
        prune_snn = math.ceil(n_neighbors*prune_snn)
    elif prune_snn >= n_neighbors:
        raise ValueError(f"invalid argument values for 'prune_snn' and 'n_neighbors: 'prune_snn' is higher than 'n_neighbors'")

    n_cells = scdata.n_obs
    distances_key = scdata.uns[cluster_key]["distances_key"]

    neighborhood_graph = scdata.obsp[distances_key].copy()
    if not issparse(neighborhood_graph):
        neighborhood_graph = csr_matrix(neighborhood_graph)
    neighborhood_graph.data[neighborhood_graph.data > 0] = 1

    neighborhood_graph = neighborhood_graph * neighborhood_graph.transpose()
    neighborhood_graph -= (n_neighbors * diags(np.ones(n_cells), offsets=0, shape=(n_cells, n_cells)))
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
    prune_snn: Union[float, int] = 1/15,
    metric: Metric = "euclidean",
    normalize_connectivities: bool = True,
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False
) -> Union[ScData, None]: # type: ignore
    """
    Compute a shared neighborhood (SNN) graph of observations.
    The neighbor search relies on a previously computed neighborhood graph (such as kNN algorithm).
    Be careful: connectivity are not well computed. This issue has to be fixed.

    Parameters
    ----------
    scdata: ad.AnnData | md.MuData
        Unimodal or multimodal annotated data matrix.
    knn_key: str (default: 'neighbors')
        If not specified, the used neighbors data are retrieved from scdata.uns['neighbors'],
        otherwise the used neighbors data are retrieved from scdata.uns['key_added'].
    snn_key: str (default: 'shared_neighbors')
        If not specified, the shared neighbors data are stored in scdata.uns['shared_neighbors']
        If specified, the shared neighbors data are added to scdata.uns['key_added'].
    prune_snn: float | int (default: 1/15)
        If zero value, no prunning is performed. If strictly positive, removes edge between two neighbors
        in the shared neighborhood graph who have a number of neighbors less than the specified value.
        Value can be relative (float between 0 and 1) or absolute (integer between 1 and k).
    metric: Metric (default: 'euclidean')
        Metric used for computing distances between two neighbors by using scdata.obsm.
    normalize_connectivities: bool (default: True)
        If false, connectivities provide the absolute number of shared neighbors (integer between 0 and k),
        otherwise provide the relative number of shared neighbors (float between 0 and k).
    distances_key: str (optional, default: None)
        If specified, distances are stored in scdata.obsp[distances_key],
        otherwise in scdata.obsp[snn_key+'_distances'].
    connectivities_key: str (optional, default: None)
        If specified, distances are stored in scdata.obsp[connectivities_key],
        otherwise in scdata.obsp[snn_key+'_connectivities'].
    copy: bool (default: False)
        Return a copy instead of writing to scdata.

    Returns
    -------
    Depending on 'copy', update 'scdata' or return ScData object.

    See 'snn_key' parameter description for the storage path of
    connectivities and distances.

    ScData contains two keys in scdata.uns['snn_key']:
        - connectivities: sparse matrix.
            Weighted adjacency matrix of the shared neighborhood graph.
            Weights should be interpreted as number of shared neighbors.
        - distances: sparse matrix of dtype 'float64'.
            Instead of decaying weights, this stores distances for each pair of neighbors.
    """

    if knn_key not in scdata.uns:
        raise KeyError("neighborhood graph not found in 'scdata': please run 'scanpy.pp.neighbors' or specify 'knn_key")
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
