#!/usr/bin/env python

from __future__ import annotations

import math
import warnings
from importlib import import_module
from importlib import util as importlib_util
from itertools import combinations
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Hashable,
    Iterable,
    Optional,
    Sequence,
    Set,
    Type,
    Union,
    cast,
)

import networkx as nx
import numpy as np
import pandas as pd
from networkx import DiGraph, Graph
from scipy.sparse import (
    csr_matrix,
    diags,
    issparse,
)

from ..._compat import Literal, get_args
from .._dependencies import require_sklearn
from .._typing import (
    AnnData,
    Metric,
    ScData,
    Shortest_Path_Method,
    anndata_or_mudata_checker,
)
from ._maths import barycenters
from ._utils import choose_representation


@require_sklearn
@anndata_or_mudata_checker
def kneighbors_graph(
    scdata: ScData,  # type: ignore
    n_neighbors: int,
    use_rep: Optional[str] = None,
    n_components: Optional[int] = None,
    metric: Metric = "euclidean",
    create_using: Type[Any] = DiGraph,
    edge_attr: str = "distance",
    index_or_name: Literal["index", "name"] = "index",
    n_jobs: int = 1,
    **metric_kwargs: Any,
) -> Graph[Any]:
    """
    Compute a k-nearest-neighbor graph from an embedding space.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    n_neighbors: int
        Number of nearest neighbors.
    use_rep: str, optional
        Representation key in `scdata.obsm`.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    create_using: Type (default: nx.DiGraph)
        Graph type to return.
    edge_attr: str (default: 'distance')
        Attribute to which the distance values are assigned on each edge.
    index_or_name: 'index' | 'name' (default: 'index')
        Whether graph nodes use integer positions or observation names.
    n_jobs: int (default: 1)
        Number of allocated processors.
    **metric_kwargs: Any
        Additional keyword arguments passed to the distance function.

    Returns
    -------
    Graph
        Weighted graph storing nearest-neighbor distances.

    Raises
    ------
    ValueError
        If `index_or_name` is not `"index"` or `"name"`.
    """

    from sklearn import neighbors

    X = choose_representation(scdata, use_rep=use_rep, n_components=n_components)
    weighted_adjacency_matrix = neighbors.kneighbors_graph(
        X=X,
        n_neighbors=n_neighbors,
        mode="distance",
        metric=metric,
        n_jobs=n_jobs,
        **metric_kwargs,
    )
    try:
        kneighbors_graph = nx.from_scipy_sparse_array(
            weighted_adjacency_matrix,
            create_using=create_using,
            edge_attribute=edge_attr,
        )
    except AttributeError:
        from_scipy_sparse_matrix = getattr(nx, "from_scipy_sparse_matrix")
        kneighbors_graph = from_scipy_sparse_matrix(
            weighted_adjacency_matrix,
            create_using=create_using,
            edge_attribute=edge_attr,
        )

    if index_or_name == "index":
        return kneighbors_graph
    elif index_or_name == "name":
        return nx.relabel_nodes(
            kneighbors_graph,
            dict(zip(list(range(len(scdata.obs.index))), scdata.obs.index)),
        )
    else:
        raise ValueError(
            f"invalid argument value for 'index_or_name': "
            f"expected 'index' or 'name' but received {index_or_name!r}"
        )


class Knnbs:
    """
    K-nearest-neighbor-based subcluster detection.

    A k-nearest-neighbor graph is constructed to compute distances between
    cells and cluster barycenters. Subclusters can be selected by maximizing
    distances to other clusters' barycenters or by minimizing distances to
    their own cluster barycenter.

    Parameters
    ----------
    n_neighbors: int or float
        Integer-valued number of nearest neighbors.
    use_rep: str (default: "X_pca")
        Representation key in `scdata.obsm`.
    n_components: int or float, optional
        Integer-valued number of dimensions to use. If None, use all dimensions.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    **metric_kwargs: Any
        Additional keyword arguments passed to the distance function.

    Raises
    ------
    TypeError
        If `n_neighbors`, `n_components` or `use_rep` has an unsupported type.
    ValueError
        If `n_neighbors`, `n_components` or `metric` has an unsupported value.
    """

    if TYPE_CHECKING:
        kneighbors_graph: Graph[Any]
        cluster_key: str
        obs: pd.Series
        shortest_path_lengths_df: pd.DataFrame

    def __init__(
        self,
        n_neighbors: Union[int, float],
        use_rep: str = "X_pca",
        n_components: Optional[Union[int, float]] = None,
        metric: Metric = "euclidean",
        **metric_kwargs: Any,
    ):

        if not isinstance(n_neighbors, (int, float)):
            raise TypeError(
                f"unsupported argument type for 'n_neighbors': "
                f"expected {int} but received {type(n_neighbors)}"
            )

        if n_neighbors <= 0:
            raise ValueError(
                f"invalid argument value for 'n_neighbors': "
                f"expected non-null positive value but received {n_neighbors!r}"
            )

        if isinstance(n_neighbors, float) and not n_neighbors.is_integer():
            raise ValueError(
                f"invalid argument value for 'n_neighbors': "
                f"expected integer but received {n_neighbors!r}"
            )

        if n_components is not None:
            if not isinstance(n_components, (int, float)):
                raise TypeError(
                    f"unsupported argument type for 'n_components': "
                    f"expected {int} but received {type(n_components)}"
                )

            if n_components <= 0:
                raise ValueError(
                    f"invalid argument value for 'n_components': "
                    f"expected non-null positive value but received {n_components!r}"
                )

            if isinstance(n_components, float) and not n_components.is_integer():
                raise ValueError(
                    f"invalid argument value for 'n_components': "
                    f"expected integer but received {n_components!r}"
                )

        if isinstance(use_rep, str):
            self.use_rep = use_rep
        else:
            raise TypeError(
                f"unsupported argument type for 'use_rep': "
                f"expected {str} but received {type(use_rep)}"
            )

        if metric in get_args(Metric):
            self.metric: Metric = metric
        else:
            raise ValueError(
                f"invalid argument value for 'metric': "
                f"expected one of {get_args(Metric)} but received {metric!r}"
            )

        self.n_neighbors: int = int(n_neighbors)
        self.n_components: Optional[int] = (
            None if n_components is None else int(n_components)
        )
        self.metric_kwargs: Dict[str, Any] = metric_kwargs

    def __repr__(self) -> str:
        """
        Return a compact representation of the estimator parameters.

        Returns
        -------
        str
            String representation of the Knnbs configuration.
        """

        return (
            f"{self.__class__.__name__}("
            f"n_neighbors={self.n_neighbors}, "
            f"use_rep={self.use_rep}, "
            f"n_components={self.n_components}, "
            f"metric={self.metric}, "
            f"metric_kwargs={self.metric_kwargs})"
        )

    @property
    def metric_kwds(self) -> Dict[str, Any]:
        """
        Deprecated alias for `metric_kwargs`.
        """

        warnings.warn(
            "`metric_kwds` is deprecated; use `metric_kwargs` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.metric_kwargs

    @metric_kwds.setter
    def metric_kwds(self, value: Dict[str, Any]) -> None:
        warnings.warn(
            "`metric_kwds` is deprecated; use `metric_kwargs` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.metric_kwargs = value

    @require_sklearn
    def fit(
        self,
        adata: AnnData,
        cluster_key: Optional[str] = None,
        n_jobs: int = 1,
        *,
        obs: Optional[str] = None,
    ) -> None:
        """
        Compute the k-nearest neighbors-based graph using an embedding space.

        The estimator parameters stored in the object are used for graph
        construction.

        Parameters
        ----------
        adata: AnnData
            Unimodal annotated data matrix.
        cluster_key: str
            Observation column in `adata.obs` defining clusters.
        n_jobs: int (default: 1)
            Number of allocated processors.
        obs: str, optional
            Deprecated alias for `cluster_key`.

        Notes
        -----
        The method updates the Knnbs object in place and adds the
        `shortest_path_lengths_df` attribute.
        """

        from sklearn.metrics import pairwise_distances

        if obs is not None:
            warnings.warn(
                "`obs` is deprecated; use `cluster_key` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if cluster_key is not None:
                raise TypeError(
                    "received both 'cluster_key' and deprecated 'obs'; "
                    "please use only 'cluster_key'"
                )
            cluster_key = obs

        if cluster_key is None:
            raise TypeError("missing required argument: 'cluster_key'")

        representation = choose_representation(
            adata, use_rep=self.use_rep, n_components=self.n_components
        )

        _kneighbors_graph = kneighbors_graph(
            adata,
            n_neighbors=self.n_neighbors,
            n_components=self.n_components,
            use_rep=self.use_rep,
            metric=self.metric,
            create_using=nx.Graph,
            index_or_name="name",
            n_jobs=n_jobs,
        )

        _barycenters = barycenters(
            adata, obs=cluster_key, use_rep=self.use_rep, n_components=self.n_components
        )
        for key, value in _barycenters.items():
            barycenter_coordinate = value.reshape(1, -1)
            distances = pairwise_distances(
                representation,
                barycenter_coordinate,
                metric=self.metric,
                n_jobs=n_jobs,
            ).reshape(1, -1)
            _kneighbors_graph.add_node(key)
            knn_indices = np.argpartition(distances, kth=self.n_neighbors, axis=1)[
                :, : self.n_neighbors
            ].reshape(-1)
            distances = list(distances[0, knn_indices].reshape(-1))
            for obs_name, distance in zip(adata.obs.index.take(knn_indices), distances):
                _kneighbors_graph.add_edge(key, obs_name, distance=distance)

        if nx.number_connected_components(_kneighbors_graph) > 1:
            warnings.warn(
                "'kneighbors_graph' not weakly connected: add edges for "
                "joining connected components"
            )
            scc = list(nx.connected_components(_kneighbors_graph))
            scc = [list(cc - set(_barycenters.keys())) for cc in scc]
            for paired_scc in combinations(scc, 2):
                adata_any = cast(Any, adata)
                x_matrix = choose_representation(
                    adata_any[paired_scc[0], :],
                    use_rep=self.use_rep,
                )
                y_matrix = choose_representation(
                    adata_any[paired_scc[1], :],
                    use_rep=self.use_rep,
                )
                dists = pairwise_distances(
                    x_matrix,
                    y_matrix,
                    metric=self.metric,
                    n_jobs=n_jobs,
                )
                i, j = np.unravel_index(np.argmin(dists), shape=dists.shape, order="C")
                _kneighbors_graph.add_edge(
                    paired_scc[0][i], paired_scc[1][j], distance=dists[i, j]
                )

        self.kneighbors_graph = _kneighbors_graph
        self.cluster_key = cluster_key
        self.obs = cast(pd.Series, adata.obs[cluster_key])
        return None

    def compute_shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:
        """
        Compute pairwise shortest path lengths between cells and barycenters.

        Parameters
        ----------
        method: 'dijkstra' | 'bellman-ford' (default: 'dijkstra')
            Algorithm used to compute the shortest path lengths.
        n_jobs: int (default: 1)
            Number of allocated processors.

        Notes
        -----
        The method updates the Knnbs object in place and adds the
        `shortest_path_lengths_df` attribute.
        """

        def shortest_path_lengths_from(
            source: str,
        ):
            distances = pd.Series(data=self.obs.index, index=self.obs.index)
            distances = distances.apply(
                lambda x: nx.shortest_path_length(
                    self.kneighbors_graph,
                    source=source,
                    target=x,
                    weight="distance",
                    method=method,
                )
            )
            return (source, distances)

        _multiprocess_is_available = (
            importlib_util.find_spec("multiprocess") is not None
        )

        if n_jobs == 1:
            shortest_path_lengths_ls = [
                shortest_path_lengths_from(category)
                for category in self.obs.cat.categories
            ]

        elif _multiprocess_is_available:
            Pool = import_module("multiprocess").Pool

            with Pool(processes=n_jobs) as pool:
                shortest_path_lengths_ls = pool.map(
                    shortest_path_lengths_from,
                    self.obs.cat.categories,
                )

        else:
            warnings.warn(
                "module multiprocess not available: not using CPU parallelization"
            )

            shortest_path_lengths_ls = [
                shortest_path_lengths_from(category)
                for category in self.obs.cat.categories
            ]

        shortest_path_lengths_dict = {}
        for k, v in shortest_path_lengths_ls:
            shortest_path_lengths_dict[k] = v

        self.shortest_path_lengths_df = pd.DataFrame.from_dict(
            data=shortest_path_lengths_dict, orient="columns"
        )
        return None

    def shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:
        """
        Deprecated alias for `compute_shortest_path_lengths`.
        """

        warnings.warn(
            "`shortest_path_lengths` is deprecated; use "
            "`compute_shortest_path_lengths` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.compute_shortest_path_lengths(method=method, n_jobs=n_jobs)

    def select_peripheral_cells(
        self,
        subcluster_size: int = 30,
        key: Optional[str] = "knnbs",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Find cells farthest from other clusters' barycenters.

        Parameters
        ----------
        subcluster_size: int (default: 30)
            Number of cells in each macrostate.
            If `subcluster_size` is greater than the cluster size, the corresponding
            subcluster contains the full cluster.
        key: str, optional (default: 'knnbs')
            Pandas Series name. If None, leave the returned Series unnamed.
        clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed.

        Returns
        -------
        pd.Series
            Subcluster labels for cells farthest from other barycenters.
        """

        clusters_iterable = cast(
            Iterable[str],
            self.obs.cat.categories if clusters is None else clusters,
        )

        subclusters_series = pd.Series(
            data=np.nan, index=self.obs.index, dtype="category", copy=True
        ).cat.add_categories(self.obs.cat.categories)
        if key is not None:
            subclusters_series.name = key

        shortest_path_length_free_from_self_cluster_df = (
            self.shortest_path_lengths_df.copy(deep=True)
        )
        for idx, obs in self.obs.items():
            if not bool(pd.isna(obs)):
                shortest_path_length_free_from_self_cluster_df.at[idx, obs] = np.nan

        _min_dists = cast(
            pd.Series,
            shortest_path_length_free_from_self_cluster_df.min(
                axis=1,
                skipna=True,
            ),
        )
        _min_dists.name = "min_dist"
        min_dists = _min_dists.to_frame().merge(
            right=self.obs.to_frame(), left_index=True, right_index=True
        )
        for cluster in clusters_iterable:
            _min_dists_cluster = cast(
                pd.Series,
                min_dists[min_dists[self.obs.name] == cluster]["min_dist"],
            )
            if len(_min_dists_cluster) < subcluster_size:
                _obs = _min_dists_cluster.index
            else:
                _idx = np.argpartition(
                    _min_dists_cluster,
                    kth=len(_min_dists_cluster) - subcluster_size,
                )[-subcluster_size:]
                _obs = _min_dists_cluster.iloc[_idx].index
            subclusters_series.loc[_obs] = cluster

        return subclusters_series

    def find_furthest_cells_to_other_barycenters(
        self,
        size: int = 30,
        key: Optional[str] = "knnbs",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Deprecated alias for `select_peripheral_cells`.
        """

        warnings.warn(
            "`find_furthest_cells_to_other_barycenters` is deprecated; use "
            "`select_peripheral_cells` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.select_peripheral_cells(
            subcluster_size=size,
            key=key,
            clusters=clusters,
        )

    def select_central_cells(
        self,
        subcluster_size: int = 30,
        key: Optional[str] = "knnbs",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Find cells closest to their own cluster barycenter.

        Parameters
        ----------
        subcluster_size: int (default: 30)
            Number of cells in each macrostate.
            If `subcluster_size` is greater than the cluster size, the corresponding
            subcluster contains the full cluster.
        key: str, optional (default: 'knnbs')
            Pandas Series name. If None, leave the returned Series unnamed.
        clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed.

        Returns
        -------
        pd.Series
            Subcluster labels for cells closest to their own barycenter.
        """

        clusters_iterable = cast(
            Iterable[str],
            self.obs.cat.categories if clusters is None else clusters,
        )

        subclusters_series = pd.Series(
            data=np.nan, index=self.obs.index, dtype="category", copy=True
        ).cat.add_categories(self.obs.cat.categories)
        if key is not None:
            subclusters_series.name = key

        _central_scores = pd.Series(data=np.nan, index=self.obs.index)
        _central_scores.name = "dist"
        for idx, obs in self.obs.items():
            if not bool(pd.isna(obs)):
                row = cast(Hashable, idx)
                column = cast(Hashable, obs)
                _central_scores.at[row] = self.shortest_path_lengths_df.at[row, column]

        _central_scores = _central_scores.to_frame().merge(
            right=self.obs.to_frame(), left_index=True, right_index=True
        )
        for cluster in clusters_iterable:
            _central_scores_cluster = cast(
                pd.Series,
                _central_scores[_central_scores[self.obs.name] == cluster]["dist"],
            )
            if len(_central_scores_cluster) < subcluster_size:
                _obs = _central_scores_cluster.index
            else:
                _idx = np.argpartition(_central_scores_cluster, kth=subcluster_size)[
                    :subcluster_size
                ]
                _obs = _central_scores_cluster.iloc[_idx].index
            subclusters_series.loc[_obs] = cluster

        return subclusters_series

    def find_closest_cells_to_self_barycenter(
        self,
        size: int = 30,
        key: Optional[str] = "knnbs",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Deprecated alias for `select_central_cells`.
        """

        warnings.warn(
            "`find_closest_cells_to_self_barycenter` is deprecated; use "
            "`select_central_cells` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.select_central_cells(
            subcluster_size=size,
            key=key,
            clusters=clusters,
        )

    def predict(
        self,
        subcluster_size: int = 30,
        key: str = "knnbs",
        peripheral_clusters: Optional[Sequence[str]] = None,
        central_clusters: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Find cluster-related cell manifolds with the KNNBS algorithm.

        Parameters
        ----------
        subcluster_size: int (default: 30)
            Number of cells in each macrostate.
            If `subcluster_size` is greater than the cluster size, the corresponding
            subcluster contains the full cluster.
        key: str (default: 'knnbs')
            Pandas Series name.
        peripheral_clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed
            by maximizing distances to other clusters' barycenters.
        central_clusters: Sequence[str] (optional, default: None)
            List of clusters for which cell subpopulations are computed
            by minimizing distances to their own barycenter.

        Returns
        -------
        pd.Series
            Subcluster labels derived from the k-nearest-neighbor-based
            subcluster algorithm.

        Raises
        ------
        RuntimeError
            If the peripheral and central cluster sets overlap.
        """

        if peripheral_clusters is None and central_clusters is None:
            peripheral_cluster_set: Set[str] = set(self.obs.cat.categories)
            central_cluster_set: Set[str] = set()
        else:
            peripheral_cluster_set = (
                set() if peripheral_clusters is None else set(peripheral_clusters)
            )
            central_cluster_set = (
                set() if central_clusters is None else set(central_clusters)
            )

        if peripheral_cluster_set & central_cluster_set:
            raise RuntimeError(
                "'peripheral_clusters' and 'central_clusters' are not disjoint"
            )

        if peripheral_cluster_set:
            peripheral_subclusters = self.select_peripheral_cells(
                subcluster_size=subcluster_size,
                key="peripheral_subclusters",
                clusters=peripheral_cluster_set,
            )
        else:
            peripheral_subclusters = pd.Series(
                data=np.nan, index=self.obs.index, dtype="category", copy=True
            ).cat.add_categories(self.obs.cat.categories)
            peripheral_subclusters.name = "peripheral_subclusters"

        if central_cluster_set:
            central_subclusters = self.select_central_cells(
                subcluster_size=subcluster_size,
                key="central_subclusters",
                clusters=central_cluster_set,
            )
        else:
            central_subclusters = pd.Series(
                data=np.nan, index=self.obs.index, dtype="category", copy=True
            ).cat.add_categories(self.obs.cat.categories)
            central_subclusters.name = "central_subclusters"

        subcluster_candidates = peripheral_subclusters.to_frame().merge(
            right=central_subclusters.to_frame(), left_index=True, right_index=True
        )
        subclusters = subcluster_candidates.bfill(axis=1).iloc[:, 0]
        subclusters.name = key

        return subclusters

    def knnbs(
        self,
        size: int = 30,
        key: str = "knnbs",
        subclusters_maximizing_distances: Optional[Sequence[str]] = None,
        subclusters_minimizing_distances: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Deprecated alias for `predict`.
        """

        warnings.warn(
            "`knnbs` is deprecated; use `predict` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.predict(
            subcluster_size=size,
            key=key,
            peripheral_clusters=subclusters_maximizing_distances,
            central_clusters=subclusters_minimizing_distances,
        )


@anndata_or_mudata_checker
def _shared_nearest_neighbors_graph(
    scdata: ScData, cluster_key: str, prune_snn: float  # type: ignore
) -> csr_matrix:

    n_neighbors = scdata.uns[cluster_key]["params"]["n_neighbors"] - 1
    if prune_snn < 0:
        raise ValueError(
            f"invalid argument value for 'prune_snn': "
            f"expected non-negative value but received {prune_snn!r}"
        )
    elif prune_snn < 1:
        prune_snn = math.ceil(n_neighbors * prune_snn)
    elif prune_snn >= n_neighbors:
        raise ValueError(
            f"invalid argument values for 'prune_snn' and 'n_neighbors': "
            f"expected prune_snn < n_neighbors but received "
            f"prune_snn={prune_snn!r} and n_neighbors={n_neighbors!r}"
        )

    n_cells = scdata.n_obs
    distances_key = scdata.uns[cluster_key]["distances_key"]

    neighborhood_graph = scdata.obsp[distances_key].copy()
    if not issparse(neighborhood_graph):
        neighborhood_graph = csr_matrix(neighborhood_graph)
    neighborhood_graph.data[neighborhood_graph.data > 0] = 1

    neighborhood_graph = neighborhood_graph * neighborhood_graph.transpose()
    neighborhood_graph -= n_neighbors * diags(
        np.ones(n_cells), offsets=0, shape=(n_cells, n_cells)
    )
    neighborhood_graph.sort_indices()
    neighborhood_graph = neighborhood_graph.astype(dtype=np.int8)

    if prune_snn:
        mask = neighborhood_graph.data <= prune_snn
        neighborhood_graph.data[mask] = 0
        neighborhood_graph.eliminate_zeros()

    return neighborhood_graph


@require_sklearn
@anndata_or_mudata_checker
def shared_neighbors(
    scdata: ScData,  # type: ignore
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[Union[float, int]] = 1 / 15,
    metric: Metric = "euclidean",
    normalize_connectivities: bool = True,
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False,
) -> Union[ScData, None]:  # type: ignore
    """
    Compute a shared-nearest-neighbor graph of observations.

    The neighbor search relies on a previously computed neighborhood graph,
    such as one produced by `scanpy.pp.neighbors`.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    knn_key: str (default: 'neighbors')
        Key in `scdata.uns` containing the source neighborhood graph metadata.
    snn_key: str (default: 'shared_neighbors')
        Key used to store shared-neighbor graph metadata in `scdata.uns`.
    prune_snn: float | int (default: 1/15)
        If zero, no pruning is performed. If strictly positive, remove edges
        whose number of shared neighbors is less than or equal to the threshold.
        Value can be relative (float between 0 and 1) or absolute
        (integer between 1 and k).
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    normalize_connectivities: bool (default: True)
        If False, connectivities store the absolute number of shared neighbors.
        Otherwise, connectivities are normalized.
    distances_key: str (optional, default: None)
        Key used to store distances in `scdata.obsp`.
    connectivities_key: str (optional, default: None)
        Key used to store connectivities in `scdata.obsp`.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Returns
    -------
    ScData or None
        Copy of `scdata` with shared-neighbor results if `copy=True`;
        otherwise None after updating `scdata` in place.

    The resulting object stores metadata in `scdata.uns[snn_key]` and matrices
    in `scdata.obsp[distances_key]` and `scdata.obsp[connectivities_key]`.

    Raises
    ------
    KeyError
        If `knn_key` is not found in `scdata.uns`.
    ValueError
        If `prune_snn` is negative or greater than or equal to the number of
        neighbors.
    """

    from sklearn.metrics import pairwise_distances

    if knn_key not in scdata.uns:
        raise KeyError(
            "neighborhood graph not found in 'scdata': "
            "please run 'scanpy.pp.neighbors' or specify 'knn_key'"
        )
    if prune_snn is None:
        prune_snn = 0
    if distances_key is None:
        distances_key = f"{snn_key}_distances"
    if connectivities_key is None:
        connectivities_key = f"{snn_key}_connectivities"
    n_neighbors = scdata.uns[knn_key]["params"]["n_neighbors"]

    scdata = scdata.copy() if copy else scdata

    snn_graph = _shared_nearest_neighbors_graph(
        scdata, cluster_key=knn_key, prune_snn=prune_snn
    )

    n_pcs = scdata.uns[knn_key]["params"]["n_pcs"]
    obsm = scdata.uns[knn_key]["params"]["use_rep"]

    X = scdata.obsm[obsm][:, 0:n_pcs]
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
        "prune_snn": (
            prune_snn if prune_snn >= 1 else math.ceil(n_neighbors * prune_snn)
        ),
        "metric": metric,
    }

    return scdata if copy else None
