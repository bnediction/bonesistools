#!/usr/bin/env python

from __future__ import annotations

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
    Tuple,
    cast,
)

import networkx as nx
import numpy as np
import pandas as pd
from networkx import Graph

from ..._compat import get_args
from ..._validation import (
    _as_literal,
    _as_non_negative_integer,
    _as_positive_integer,
    _as_string,
)
from ..._warnings import _warn_deprecated, _warn_deprecated_argument
from .._typing import AnnData, Metric, Shortest_Path_Method
from ._neighbors import knn_graph
from ._utils import (
    _UNSET,
    _resolve_representation_argument,
    get_representation,
)


class KNNSC:
    """
    K-nearest-neighbor-based subcluster detection.

    A k-nearest-neighbor graph is constructed to compute distances between
    cells and cluster barycenters. Subclusters can be selected by maximizing
    distances to other clusters' barycenters or by minimizing distances to
    their own cluster barycenter.

    Configure and run the algorithm with `fit()`.
    """

    if TYPE_CHECKING:
        _knn_graph: Graph[Any]
        _cluster_key: str
        _obs: pd.Series
        _shortest_path_lengths_df: pd.DataFrame
        _min_cluster_size: int
        _cluster_counts: pd.Series
        _n_neighbors: int
        _representation: str
        _n_components: Optional[int]
        _metric: Metric
        _method: Shortest_Path_Method

    def __init__(
        self,
        representation: Any = _UNSET,
        n_components: Optional[int] = None,
        n_neighbors: Optional[int] = None,
        metric: Optional[Metric] = None,
        *,
        use_rep: Any = _UNSET,
        **metric_kwargs: Any,
    ):
        self._deprecated_init_params: Dict[str, Any] = {}
        self._metric_kwargs: Dict[str, Any] = {}

        representation = _resolve_representation_argument(
            representation,
            use_rep,
            default=None,
            stacklevel=2,
        )

        if (
            n_neighbors is not None
            or representation is not None
            or n_components is not None
            or metric is not None
            or metric_kwargs
        ):
            warnings.warn(
                "Passing KNNSC configuration parameters to `__init__` is "
                "deprecated and will be removed in 2.0.0; pass them to "
                "`KNNSC.fit(...)` instead.",
                FutureWarning,
                stacklevel=2,
            )
            (
                resolved_representation,
                resolved_n_components,
                resolved_n_neighbors,
                resolved_metric,
                resolved_metric_kwargs,
            ) = _normalize_knnsc_configuration(
                representation=representation,
                n_components=n_components,
                n_neighbors=n_neighbors,
                metric=metric,
                metric_kwargs=metric_kwargs,
            )
            self._deprecated_init_params = {
                "representation": resolved_representation,
                "n_components": resolved_n_components,
                "n_neighbors": resolved_n_neighbors,
                "metric": resolved_metric,
                "metric_kwargs": resolved_metric_kwargs,
            }
            self._metric_kwargs = resolved_metric_kwargs

    def __repr__(self) -> str:
        """
        Return a compact representation of the estimator parameters.

        Returns
        -------
            str
            String representation of the KNNSC configuration.
        """

        if hasattr(self, "_cluster_key"):
            fit_params = self.params_
            return (
                f"{self.__class__.__name__}("
                f"cluster_key={fit_params['cluster_key']}, "
                f"representation={fit_params['representation']}, "
                f"n_components={fit_params['n_components']}, "
                f"n_neighbors={fit_params['n_neighbors']}, "
                f"metric={fit_params['metric']}, "
                f"metric_kwargs={fit_params['metric_kwargs']}, "
                f"min_cluster_size={fit_params['min_cluster_size']}, "
                f"method={fit_params['method']})"
            )

        if not self._deprecated_init_params:
            return f"{self.__class__.__name__}()"

        return (
            f"{self.__class__.__name__}("
            f"representation={self._deprecated_init_params['representation']}, "
            f"n_components={self._deprecated_init_params['n_components']}, "
            f"n_neighbors={self._deprecated_init_params['n_neighbors']}, "
            f"metric={self._deprecated_init_params['metric']}, "
            f"metric_kwargs={self._deprecated_init_params['metric_kwargs']})"
        )

    @property
    def knn_graph(self) -> Graph[Any]:
        """
        Fitted k-nearest-neighbor graph.
        """

        return cast(Graph, self._require_fitted_attribute("_knn_graph"))

    @property
    def cluster_key(self) -> str:
        """
        Observation column used as cluster labels during fitting.
        """

        return cast(str, self._require_fitted_attribute("_cluster_key"))

    @property
    def min_cluster_size(self) -> int:
        """
        Minimum cluster size used during fitting.
        """

        return cast(int, self._require_fitted_attribute("_min_cluster_size"))

    @property
    def cluster_counts(self) -> pd.Series:
        """
        Cluster sizes observed during fitting.
        """

        return cast(pd.Series, self._require_fitted_attribute("_cluster_counts"))

    @property
    def obs(self) -> pd.Series:
        """
        Fitted cluster labels after removing ineligible categories.
        """

        return cast(pd.Series, self._require_fitted_attribute("_obs"))

    @property
    def shortest_path_lengths_df(self) -> pd.DataFrame:
        """
        Fitted shortest-path lengths between cells and barycenters.
        """

        return cast(
            pd.DataFrame,
            self._require_fitted_attribute("_shortest_path_lengths_df"),
        )

    @property
    def n_neighbors(self) -> int:
        """
        Number of nearest neighbors used during fitting.
        """

        return cast(int, self._require_fitted_attribute("_n_neighbors"))

    @property
    def representation(self) -> str:
        """
        Representation key used during fitting.
        """

        return cast(str, self._require_fitted_attribute("_representation"))

    @property
    def use_rep(self) -> str:
        """
        Deprecated alias for `representation`.
        """

        _warn_deprecated_argument("use_rep", "representation", stacklevel=2)
        return self.representation

    @property
    def n_components(self) -> Optional[int]:
        """
        Number of representation dimensions used during fitting.
        """

        return cast(Optional[int], self._require_fitted_attribute("_n_components"))

    @property
    def metric(self) -> Metric:
        """
        Distance metric used during fitting.
        """

        return cast(Metric, self._require_fitted_attribute("_metric"))

    @property
    def metric_kwargs(self) -> Dict[str, Any]:
        """
        Additional distance metric parameters used during fitting.
        """

        if hasattr(self, "_metric_kwargs"):
            return self._metric_kwargs
        return cast(
            Dict[str, Any],
            self._deprecated_init_params.get("metric_kwargs", {}),
        )

    @property
    def method(self) -> Shortest_Path_Method:
        """
        Shortest-path algorithm used during fitting.
        """

        return cast(Shortest_Path_Method, self._require_fitted_attribute("_method"))

    @property
    def params_(self) -> Dict[str, Any]:
        """
        Fitted KNNSC configuration.
        """

        return {
            "cluster_key": self.cluster_key,
            "representation": self.representation,
            "n_components": self.n_components,
            "n_neighbors": self.n_neighbors,
            "metric": self.metric,
            "metric_kwargs": self.metric_kwargs,
            "min_cluster_size": self.min_cluster_size,
            "method": self.method,
        }

    @property
    def metric_kwds(self) -> Dict[str, Any]:
        """
        Deprecated alias for `metric_kwargs`.
        """

        _warn_deprecated_argument("metric_kwds", "metric_kwargs", stacklevel=2)
        return self.metric_kwargs

    @metric_kwds.setter
    def metric_kwds(self, value: Dict[str, Any]) -> None:
        _warn_deprecated_argument("metric_kwds", "metric_kwargs", stacklevel=2)
        self._metric_kwargs = value
        if self._deprecated_init_params:
            self._deprecated_init_params["metric_kwargs"] = value

    @property
    def kneighbors_graph(self) -> Graph[Any]:
        """
        Deprecated alias for `knn_graph`.
        """

        _warn_deprecated(
            "`kneighbors_graph`",
            replacement="`knn_graph`",
            stacklevel=2,
        )
        return self.knn_graph

    @kneighbors_graph.setter
    def kneighbors_graph(self, value: Graph[Any]) -> None:
        _warn_deprecated(
            "`kneighbors_graph`",
            replacement="`knn_graph`",
            stacklevel=2,
        )
        self._knn_graph = value

    def fit(
        self,
        adata: AnnData,
        cluster_key: Optional[str] = None,
        representation: Any = _UNSET,
        n_components: Optional[int] = None,
        n_neighbors: Optional[int] = None,
        metric: Optional[Metric] = None,
        metric_kwargs: Optional[Dict[str, Any]] = None,
        n_jobs: int = 1,
        method: Shortest_Path_Method = "dijkstra",
        *,
        obs: Optional[str] = None,
        use_rep: Any = _UNSET,
        min_cluster_size: int = 1,
    ) -> None:
        """
        Fit the KNNSC estimator using an embedding space.

        The method builds the k-nearest-neighbor graph, adds cluster
        barycenters, connects graph components when necessary, then computes
        `shortest_path_lengths_df`. All cells are kept in the graph, but only
        clusters containing at least `min_cluster_size` cells are used as
        barycenter sources for shortest-path computations.

        Parameters
        ----------
        adata: AnnData
            Unimodal annotated data matrix.
        cluster_key: str
            Observation column in `adata.obs` defining clusters.
        representation: str, optional (default: "X_pca")
            Representation key in `adata.obsm`.
        n_components: int, optional
            Number of dimensions to use. If `None`, use all dimensions.
        n_neighbors: int
            Number of nearest neighbors.
        metric: Metric (default: 'euclidean')
            Metric used when calculating pairwise distances between observations.
        metric_kwargs: dict, optional
            Additional keyword arguments passed to the distance function.
        n_jobs: int (default: 1)
            Number of allocated processors.
        method: 'dijkstra' | 'bellman-ford' (default: 'dijkstra')
            Algorithm used to compute the shortest path lengths.
        obs: str, optional
            Deprecated alias for `cluster_key`.
        use_rep: str, optional
            Deprecated alias for `representation`.
        min_cluster_size: int (default: 1)
            Minimum number of cells required for a non-empty cluster to be used
            as a barycenter source. Smaller clusters remain in the graph but
            are ignored as candidate subclusters. A value of 0 disables
            size-based filtering but still ignores empty categories.

        Notes
        -----
        The method updates the KNNSC object in place and adds the
        `knn_graph` and `shortest_path_lengths_df` attributes.
        """

        from sklearn.metrics import pairwise_distances

        if obs is not None:
            _warn_deprecated_argument("obs", "cluster_key", stacklevel=2)
            if cluster_key is not None:
                raise TypeError(
                    "received both 'cluster_key' and deprecated 'obs'; "
                    "please use only 'cluster_key'"
                )
            cluster_key = obs

        if cluster_key is None:
            raise TypeError("missing required argument: 'cluster_key'")

        representation = _resolve_representation_argument(
            representation,
            use_rep,
            default=None,
            stacklevel=2,
        )

        init_params = self._deprecated_init_params
        if n_neighbors is None and "n_neighbors" in init_params:
            n_neighbors = cast(int, init_params["n_neighbors"])
        if representation is None:
            representation = cast(Optional[str], init_params.get("representation"))
        if n_components is None:
            n_components = cast(Optional[int], init_params.get("n_components"))
        if metric is None:
            metric = cast(Optional[Metric], init_params.get("metric"))
        if metric_kwargs is None:
            metric_kwargs = cast(
                Optional[Dict[str, Any]], init_params.get("metric_kwargs")
            )

        (
            resolved_representation,
            resolved_n_components,
            resolved_n_neighbors,
            resolved_metric,
            resolved_metric_kwargs,
        ) = _normalize_knnsc_configuration(
            representation=representation,
            n_components=n_components,
            n_neighbors=n_neighbors,
            metric=metric,
            metric_kwargs=metric_kwargs,
        )

        min_cluster_size = _as_non_negative_integer(
            min_cluster_size,
            "min_cluster_size",
        )

        raw_obs = cast(pd.Series, adata.obs[cluster_key])
        if not hasattr(raw_obs, "cat"):
            raise AttributeError(
                f"adata.obs[{cluster_key!r}] object has no attribute 'cat'"
            )

        cluster_counts = cast(
            pd.Series,
            raw_obs.value_counts().reindex(raw_obs.cat.categories).fillna(0),
        )
        ineligible_categories = [
            category
            for category in raw_obs.cat.categories
            if int(cluster_counts.loc[category]) == 0
            or int(cluster_counts.loc[category]) < min_cluster_size
        ]
        candidate_obs = raw_obs.cat.remove_categories(ineligible_categories)
        if len(candidate_obs.cat.categories) == 0:
            raise ValueError(
                "no non-empty clusters contain at least "
                f"min_cluster_size={min_cluster_size} cells"
            )

        representation_mtx = get_representation(
            adata,
            obsm=resolved_representation,
            n_components=resolved_n_components,
        )

        _knn_graph = knn_graph(
            adata,
            representation=resolved_representation,
            n_components=resolved_n_components,
            n_neighbors=resolved_n_neighbors,
            metric=resolved_metric,
            create_using=nx.Graph,
            index_or_name="name",
            n_jobs=n_jobs,
            **resolved_metric_kwargs,
        )

        _barycenters = {
            category: np.nanmean(
                cast(Any, representation_mtx)[raw_obs == category],
                axis=0,
            )
            for category in candidate_obs.cat.categories
        }
        for key, value in _barycenters.items():
            barycenter_coordinate = value.reshape(1, -1)
            distances = pairwise_distances(
                representation_mtx,
                barycenter_coordinate,
                metric=resolved_metric,
                n_jobs=n_jobs,
                **resolved_metric_kwargs,
            ).reshape(1, -1)
            _knn_graph.add_node(key)
            knn_indices = np.argpartition(distances, kth=resolved_n_neighbors, axis=1)[
                :, :resolved_n_neighbors
            ].reshape(-1)
            distances = list(distances[0, knn_indices].reshape(-1))
            for obs_name, distance in zip(adata.obs.index.take(knn_indices), distances):
                _knn_graph.add_edge(key, obs_name, distance=distance)

        if nx.number_connected_components(_knn_graph) > 1:
            warnings.warn(
                "'knn_graph' not weakly connected: add edges for "
                "joining connected components"
            )
            scc = list(nx.connected_components(_knn_graph))
            scc = [list(cc - set(_barycenters.keys())) for cc in scc]
            for paired_scc in combinations(scc, 2):
                adata_any = cast(Any, adata)
                x_representation_mtx = get_representation(
                    adata_any[paired_scc[0], :],
                    obsm=resolved_representation,
                    n_components=resolved_n_components,
                )
                y_representation_mtx = get_representation(
                    adata_any[paired_scc[1], :],
                    obsm=resolved_representation,
                    n_components=resolved_n_components,
                )
                dists = pairwise_distances(
                    x_representation_mtx,
                    y_representation_mtx,
                    metric=resolved_metric,
                    n_jobs=n_jobs,
                    **resolved_metric_kwargs,
                )
                i, j = np.unravel_index(np.argmin(dists), shape=dists.shape, order="C")
                _knn_graph.add_edge(
                    paired_scc[0][i], paired_scc[1][j], distance=dists[i, j]
                )

        self._knn_graph = _knn_graph
        self._cluster_key = cluster_key
        self._min_cluster_size = min_cluster_size
        self._cluster_counts = cluster_counts
        self._obs = cast(pd.Series, candidate_obs)
        self._representation = resolved_representation
        self._n_components = resolved_n_components
        self._n_neighbors = resolved_n_neighbors
        self._metric = resolved_metric
        self._metric_kwargs = resolved_metric_kwargs
        self._method = method
        self._compute_shortest_path_lengths(method=method, n_jobs=n_jobs)
        return None

    def compute_shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:
        """
        Deprecated. Recompute pairwise shortest path lengths.

        Use `KNNSC.fit(...)` instead. This method is kept temporarily for
        compatibility and updates `shortest_path_lengths_df` in place.
        """

        warnings.warn(
            "`compute_shortest_path_lengths()` is deprecated and no longer "
            "required after `fit()`. Calling it explicitly recomputes "
            "`shortest_path_lengths_df`.",
            FutureWarning,
            stacklevel=2,
        )
        return self._compute_shortest_path_lengths(method=method, n_jobs=n_jobs)

    def shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:
        """
        Deprecated alias for `compute_shortest_path_lengths`.
        """

        _warn_deprecated(
            "`shortest_path_lengths`",
            replacement="`KNNSC.fit(...)`",
            stacklevel=2,
        )
        return self._compute_shortest_path_lengths(method=method, n_jobs=n_jobs)

    def select_peripheral_cells(
        self,
        subcluster_size: int = 30,
        key: Optional[str] = "knnsc",
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
        key: str, optional (default: 'knnsc')
            Pandas Series name. If `None`, leave the returned Series unnamed.
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
        clusters_iterable = list(clusters_iterable)
        self._validate_candidate_clusters(clusters_iterable, "clusters")

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
            cast(Any, subclusters_series.loc)[_obs] = cluster

        return subclusters_series

    def find_furthest_cells_to_other_barycenters(
        self,
        size: int = 30,
        key: Optional[str] = "knnsc",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Deprecated alias for `select_peripheral_cells`.
        """

        _warn_deprecated(
            "`find_furthest_cells_to_other_barycenters`",
            replacement="`select_peripheral_cells`",
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
        key: Optional[str] = "knnsc",
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
        key: str, optional (default: 'knnsc')
            Pandas Series name. If `None`, leave the returned Series unnamed.
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
        clusters_iterable = list(clusters_iterable)
        self._validate_candidate_clusters(clusters_iterable, "clusters")

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
            cast(Any, subclusters_series.loc)[_obs] = cluster

        return subclusters_series

    def find_closest_cells_to_self_barycenter(
        self,
        size: int = 30,
        key: Optional[str] = "knnsc",
        clusters: Optional[Iterable[str]] = None,
    ) -> pd.Series:
        """
        Deprecated alias for `select_central_cells`.
        """

        _warn_deprecated(
            "`find_closest_cells_to_self_barycenter`",
            replacement="`select_central_cells`",
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
        key: str = "knnsc",
        peripheral_clusters: Optional[Sequence[str]] = None,
        central_clusters: Optional[Sequence[str]] = None,
    ) -> pd.Series:
        """
        Find cluster-related cell manifolds with the KNNSC algorithm.

        Call `fit()` before calling this method.

        Parameters
        ----------
        subcluster_size: int (default: 30)
            Number of cells in each macrostate.
            If `subcluster_size` is greater than the cluster size, the corresponding
            subcluster contains the full cluster.
        key: str (default: 'knnsc')
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
        self._validate_candidate_clusters(peripheral_cluster_set, "peripheral_clusters")
        self._validate_candidate_clusters(central_cluster_set, "central_clusters")

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
        if hasattr(subclusters, "cat"):
            subclusters = subclusters.cat.remove_unused_categories()
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

        _warn_deprecated("`knnbs`", replacement="`predict`", stacklevel=2)
        return self.predict(
            subcluster_size=size,
            key=key,
            peripheral_clusters=subclusters_maximizing_distances,
            central_clusters=subclusters_minimizing_distances,
        )

    def _require_fitted_attribute(self, name: str) -> Any:

        if not hasattr(self, name):
            raise AttributeError("KNNSC has not been fitted yet")
        return getattr(self, name)

    def _compute_shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:

        def shortest_path_lengths_from(
            source: Hashable,
        ):
            if method == "dijkstra":
                lengths = nx.single_source_dijkstra_path_length(
                    self.knn_graph,
                    source=source,
                    weight="distance",
                )
            elif method == "bellman-ford":
                lengths = nx.single_source_bellman_ford_path_length(
                    self.knn_graph,
                    source=source,
                    weight="distance",
                )
            else:
                raise ValueError(f"unsupported shortest path method: {method!r}")

            distances = pd.Series(lengths, dtype=float).reindex(self.obs.index)
            missing_mask = distances.isna().to_numpy(dtype=bool)
            if np.any(missing_mask):
                missing = np.asarray(distances.index, dtype=object)[missing_mask][0]
                raise nx.NetworkXNoPath(
                    f"node {missing!r} is not reachable from {source!r}"
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

        self._shortest_path_lengths_df = pd.DataFrame.from_dict(
            data=shortest_path_lengths_dict, orient="columns"
        )
        return None

    def _validate_candidate_clusters(
        self,
        clusters: Iterable[str],
        argument: str,
    ) -> None:

        candidate_clusters = set(self.obs.cat.categories)
        invalid_clusters = [
            cluster for cluster in clusters if cluster not in candidate_clusters
        ]
        if not invalid_clusters:
            return None

        below_minimum = []
        empty = []
        missing = []
        cluster_counts = cast(pd.Series, self.cluster_counts)
        for cluster in invalid_clusters:
            if cluster in cluster_counts.index:
                size = int(cluster_counts.loc[cluster])
                if size == 0:
                    empty.append(str(cluster))
                else:
                    below_minimum.append(
                        f"{cluster} "
                        f"(size={size}, "
                        f"min_cluster_size={self.min_cluster_size})"
                    )
            else:
                missing.append(str(cluster))

        details = []
        if empty:
            details.append("empty: " + ", ".join(empty))
        if below_minimum:
            details.append("below minimum size: " + ", ".join(below_minimum))
        if missing:
            details.append("not found: " + ", ".join(missing))

        raise ValueError(f"invalid cluster(s) in '{argument}': {'; '.join(details)}")


class Knnbs(KNNSC):
    """
    Deprecated alias for `KNNSC`.
    """

    def __init__(
        self,
        n_neighbors: int,
        representation: Any = _UNSET,
        n_components: Optional[int] = None,
        metric: Metric = "euclidean",
        *,
        use_rep: Any = _UNSET,
        **metric_kwargs: Any,
    ):
        _warn_deprecated(
            "`bt.omics.tl.Knnbs`",
            replacement="`bt.omics.tl.KNNSC`",
            stacklevel=2,
        )
        representation = cast(
            str,
            _resolve_representation_argument(
                representation,
                use_rep,
                default="X_pca",
                stacklevel=2,
            ),
        )
        super().__init__(
            n_neighbors=n_neighbors,
            representation=representation,
            n_components=n_components,
            metric=metric,
            **metric_kwargs,
        )


def _normalize_knnsc_configuration(
    representation: Optional[str],
    n_components: Optional[int],
    n_neighbors: Optional[int],
    metric: Optional[Metric],
    metric_kwargs: Optional[Dict[str, Any]],
) -> Tuple[str, Optional[int], int, Metric, Dict[str, Any]]:

    if n_components is not None:
        n_components = _as_positive_integer(n_components, "n_components")

    resolved_representation = "X_pca" if representation is None else representation
    resolved_representation = _as_string(
        resolved_representation,
        "representation",
    )

    if n_neighbors is None:
        raise TypeError("missing required argument: 'n_neighbors'")

    n_neighbors = _as_positive_integer(n_neighbors, "n_neighbors")

    resolved_metric: Metric = "euclidean" if metric is None else metric
    resolved_metric = _as_literal(
        resolved_metric,
        choices=get_args(Metric),
        name="metric",
    )

    if metric_kwargs is None:
        resolved_metric_kwargs: Dict[str, Any] = {}
    elif isinstance(metric_kwargs, dict):
        resolved_metric_kwargs = dict(metric_kwargs)
    else:
        raise TypeError(
            f"unsupported argument type for 'metric_kwargs': "
            f"expected {dict} but received {type(metric_kwargs)}"
        )

    return (
        resolved_representation,
        n_components,
        n_neighbors,
        resolved_metric,
        resolved_metric_kwargs,
    )
