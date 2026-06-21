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
    Tuple,
    Type,
    Union,
    cast,
    overload,
)

import networkx as nx
import numpy as np
import pandas as pd
from networkx import DiGraph, Graph
from scipy.sparse import (
    coo_matrix,
    csr_matrix,
    diags,
    issparse,
)

from ..._compat import Literal, get_args
from ..._typing import RandomStateSeed
from ..._validation import (
    _as_literal,
    _as_non_negative_integer,
    _as_non_negative_number,
    _as_positive_integer,
    _as_seed,
    _as_string,
)
from ..._warnings import _warn_deprecated, _warn_deprecated_argument
from .._metadata import _format_random_state
from .._typing import (
    AnnData,
    Matrix,
    Metric,
    ScData,
    Shortest_Path_Method,
    anndata_or_mudata_checker,
)
from ._utils import get_representation

IndexOrName = Literal["index", "name"]
NeighborBackend = Literal["exact", "pynndescent"]
_SMOOTH_K_TOLERANCE = 1e-5
_MIN_K_DIST_SCALE = 0.001
_NPY_FLOATMAX = np.float32(3.4028235e38)


def _kneighbors_distance_matrix(
    scdata: ScData,
    use_rep: Optional[str],
    n_components: Optional[int],
    n_neighbors: int,
    metric: Metric,
    n_jobs: int,
) -> csr_matrix:

    from sklearn import neighbors as sklearn_neighbors

    representation_mtx = cast(
        Matrix,
        get_representation(
            scdata,
            use_rep=use_rep,
            n_components=n_components,
        ),
    )
    matrix = sklearn_neighbors.kneighbors_graph(
        X=representation_mtx,
        n_neighbors=n_neighbors,
        mode="distance",
        metric=metric,
        n_jobs=n_jobs,
    )
    return cast(csr_matrix, matrix)


def _kneighbors_graph_matrices(
    scdata: ScData,
    use_rep: Optional[str],
    n_components: Optional[int],
    n_neighbors: int,
    backend: NeighborBackend,
    metric: Metric,
    seed: RandomStateSeed,
    n_jobs: int,
) -> Tuple[csr_matrix, csr_matrix]:

    representation_mtx = cast(
        Matrix,
        get_representation(
            scdata,
            use_rep=use_rep,
            n_components=n_components,
        ),
    )
    if backend == "exact":
        knn_indices, knn_distances = _exact_knn_arrays(
            representation_mtx=representation_mtx,
            n_neighbors=n_neighbors,
            metric=metric,
            n_jobs=n_jobs,
        )
    else:
        knn_indices, knn_distances = _pynndescent_knn_arrays(
            representation_mtx=representation_mtx,
            n_neighbors=n_neighbors,
            metric=metric,
            seed=seed,
            n_jobs=n_jobs,
        )

    distances = _sparse_distances_from_knn(
        knn_indices=knn_indices,
        knn_distances=knn_distances,
        n_obs=scdata.n_obs,
    )
    connectivities = _umap_connectivities_from_knn(
        knn_indices,
        knn_distances,
        n_obs=scdata.n_obs,
    )
    return distances, connectivities


def _get_pynndescent_transformer() -> Any:

    try:
        from pynndescent import PyNNDescentTransformer
    except ImportError as exc:
        raise ImportError(
            "backend='pynndescent' requires the optional dependency "
            "'pynndescent'. Install it with:\n\n"
            "pip install pynndescent\n\n"
            "or use backend='exact'."
        ) from exc

    return PyNNDescentTransformer


def _exact_knn_arrays(
    representation_mtx: Matrix,
    n_neighbors: int,
    metric: Metric,
    n_jobs: int,
) -> Tuple[np.ndarray, np.ndarray]:

    from sklearn import neighbors as sklearn_neighbors

    neighbors_model = sklearn_neighbors.NearestNeighbors(
        n_neighbors=n_neighbors,
        algorithm="brute",
        metric=metric,
        n_jobs=n_jobs,
    )
    sklearn_representation_mtx = cast(Any, representation_mtx)
    neighbors_model.fit(sklearn_representation_mtx)
    knn_distances, knn_indices = neighbors_model.kneighbors(sklearn_representation_mtx)
    return (
        knn_indices.astype(np.int32, copy=False),
        knn_distances.astype(np.float32, copy=False),
    )


def _pynndescent_knn_arrays(
    representation_mtx: Matrix,
    n_neighbors: int,
    metric: Metric,
    seed: RandomStateSeed,
    n_jobs: int,
) -> Tuple[np.ndarray, np.ndarray]:

    PyNNDescentTransformer = _get_pynndescent_transformer()
    transformer = PyNNDescentTransformer(
        n_neighbors=n_neighbors,
        metric=metric,
        random_state=_as_seed(seed),
        n_jobs=n_jobs,
    )
    sparse_distances = cast(csr_matrix, transformer.fit_transform(representation_mtx))
    return _knn_arrays_from_sparse_distances(
        sparse_distances=sparse_distances.tocsr(),
        n_neighbors=n_neighbors,
    )


def _knn_arrays_from_sparse_distances(
    sparse_distances: csr_matrix,
    n_neighbors: int,
) -> Tuple[np.ndarray, np.ndarray]:

    row_sizes = np.diff(sparse_distances.indptr)
    if np.any(row_sizes < n_neighbors):
        raise ValueError(
            "invalid neighbor graph: backend returned fewer neighbors than expected"
        )

    if np.all(row_sizes == row_sizes[0]):
        row_size = int(row_sizes[0])
        n_obs = cast(Tuple[int, int], sparse_distances.shape)[0]
        indices = sparse_distances.indices.reshape(n_obs, row_size)
        distances = sparse_distances.data.reshape(n_obs, row_size)
        order = np.argsort(distances, axis=1)
        indices = np.take_along_axis(indices, order, axis=1)
        distances = np.take_along_axis(distances, order, axis=1)

        rows = np.arange(indices.shape[0])
        if np.all((indices[:, 0] == rows) & (distances[:, 0] == 0.0)):
            return (
                indices[:, :n_neighbors].astype(np.int32, copy=False),
                distances[:, :n_neighbors].astype(np.float32, copy=False),
            )

    return _knn_arrays_from_sparse_distances_slow(
        sparse_distances=sparse_distances,
        n_neighbors=n_neighbors,
    )


def _knn_arrays_from_sparse_distances_slow(
    sparse_distances: csr_matrix,
    n_neighbors: int,
) -> Tuple[np.ndarray, np.ndarray]:

    n_obs = cast(Tuple[int, int], sparse_distances.shape)[0]
    knn_indices = np.empty((n_obs, n_neighbors), dtype=np.int32)
    knn_distances = np.empty((n_obs, n_neighbors), dtype=np.float32)
    for row in range(n_obs):
        start = sparse_distances.indptr[row]
        end = sparse_distances.indptr[row + 1]
        indices = sparse_distances.indices[start:end]
        distances = sparse_distances.data[start:end]

        order = np.argsort(distances)
        indices = indices[order]
        distances = distances[order]

        self_positions = np.flatnonzero((indices == row) & (distances == 0.0))
        if self_positions.size == 0:
            indices = np.concatenate(([row], indices))
            distances = np.concatenate(([0.0], distances))
        elif self_positions[0] != 0:
            self_position = int(self_positions[0])
            indices = np.concatenate(
                (
                    indices[self_position : self_position + 1],
                    indices[:self_position],
                    indices[self_position + 1 :],
                )
            )
            distances = np.concatenate(
                (
                    distances[self_position : self_position + 1],
                    distances[:self_position],
                    distances[self_position + 1 :],
                )
            )

        if indices.shape[0] < n_neighbors:
            raise ValueError(
                "invalid neighbor graph: backend returned fewer neighbors than expected"
            )

        knn_indices[row, :] = indices[:n_neighbors]
        knn_distances[row, :] = distances[:n_neighbors]

    return knn_indices, knn_distances


def _sparse_distances_from_knn(
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
    n_obs: int,
) -> csr_matrix:

    n_neighbors = knn_indices.shape[1]
    rows = np.repeat(np.arange(n_obs), n_neighbors)
    columns = knn_indices.ravel()
    data = knn_distances.ravel()
    non_self = rows != columns
    distances = csr_matrix(
        (data[non_self], (rows[non_self], columns[non_self])),
        shape=(n_obs, n_obs),
    )
    distances.eliminate_zeros()
    return distances


def _umap_smooth_knn_distances(
    distances: np.ndarray,
    n_neighbors: int,
) -> Tuple[np.ndarray, np.ndarray]:

    distances = distances.astype(np.float32, copy=False)
    target = np.float32(np.log2(float(n_neighbors)))
    n_obs = distances.shape[0]

    rhos = np.zeros(n_obs, dtype=np.float32)
    for i in range(n_obs):
        non_zero_distances = distances[i][distances[i] > 0.0]
        if non_zero_distances.shape[0] > 0:
            rhos[i] = non_zero_distances[0]

    lo = np.zeros(n_obs, dtype=np.float32)
    hi = np.full(n_obs, _NPY_FLOATMAX, dtype=np.float32)
    sigmas = np.ones(n_obs, dtype=np.float32)
    active = np.ones(n_obs, dtype=bool)
    shifted_distances = distances[:, 1:] - rhos[:, None]

    for _ in range(64):
        with np.errstate(over="ignore"):
            probabilities = np.where(
                shifted_distances > 0.0,
                np.exp(-(shifted_distances / sigmas[:, None])),
                1.0,
            )
        probability_sums = probabilities.sum(axis=1)
        close = np.abs(probability_sums - target) < _SMOOTH_K_TOLERANCE
        update = active & ~close
        if not np.any(update):
            break

        greater = update & (probability_sums > target)
        hi[greater] = sigmas[greater]
        sigmas[greater] = (lo[greater] + hi[greater]) / 2.0

        lower = update & ~greater
        lo[lower] = sigmas[lower]
        unset_hi = lower & (hi >= _NPY_FLOATMAX)
        sigmas[unset_hi] *= 2.0
        set_hi = lower & ~unset_hi
        sigmas[set_hi] = (lo[set_hi] + hi[set_hi]) / 2.0

        active &= ~close

    mean_distances = np.mean(distances)
    mean_row_distances = np.mean(distances, axis=1)
    min_sigmas = np.where(
        rhos > 0.0,
        _MIN_K_DIST_SCALE * mean_row_distances,
        _MIN_K_DIST_SCALE * mean_distances,
    )
    sigmas = np.maximum(sigmas, min_sigmas)

    return sigmas.astype(np.float32, copy=False), rhos


def _umap_connectivities_from_knn(
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
    n_obs: int,
) -> csr_matrix:

    knn_distances = knn_distances.astype(np.float32, copy=False)
    n_neighbors = knn_indices.shape[1]
    sigmas, rhos = _umap_smooth_knn_distances(knn_distances, n_neighbors)

    rows = np.repeat(np.arange(n_obs, dtype=np.int32), n_neighbors)
    columns = knn_indices.astype(np.int32, copy=False).ravel()
    shifted_distances = knn_distances - rhos[:, None]
    with np.errstate(over="ignore"):
        values = np.exp(-(shifted_distances / sigmas[:, None]))
    values = np.where(
        (shifted_distances <= 0.0) | (sigmas[:, None] == 0.0),
        1.0,
        values,
    ).astype(np.float32, copy=False)
    values[knn_indices == np.arange(n_obs)[:, None]] = 0.0

    valid = columns != -1
    graph = coo_matrix(
        (
            values.ravel()[valid],
            (rows[valid], columns[valid]),
        ),
        shape=(n_obs, n_obs),
    ).tocsr()

    transpose = graph.transpose()
    intersections = graph.multiply(transpose)
    graph = graph + transpose - intersections
    graph.eliminate_zeros()

    return cast(csr_matrix, graph.tocsr())


def _normalize_knnsc_configuration(
    use_rep: Optional[str],
    n_components: Optional[int],
    n_neighbors: Optional[int],
    metric: Optional[Metric],
    metric_kwargs: Optional[Dict[str, Any]],
) -> Tuple[str, Optional[int], int, Metric, Dict[str, Any]]:

    if n_components is not None:
        n_components = _as_positive_integer(n_components, "n_components")

    resolved_use_rep = "X_pca" if use_rep is None else use_rep
    resolved_use_rep = _as_string(resolved_use_rep, "use_rep")

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
        resolved_use_rep,
        n_components,
        n_neighbors,
        resolved_metric,
        resolved_metric_kwargs,
    )


@overload
def neighbors(
    scdata: ScData,
    n_neighbors: int = 15,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    backend: NeighborBackend = "exact",
    metric: Metric = "euclidean",
    key_added: str = "neighbors",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    copy: Literal[False] = False,
) -> None: ...


@overload
def neighbors(
    scdata: ScData,
    n_neighbors: int = 15,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    backend: NeighborBackend = "exact",
    metric: Metric = "euclidean",
    key_added: str = "neighbors",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    *,
    copy: Literal[True],
) -> ScData: ...


@overload
def neighbors(
    scdata: ScData,
    n_neighbors: int = 15,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    backend: NeighborBackend = "exact",
    metric: Metric = "euclidean",
    key_added: str = "neighbors",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    copy: bool = False,
) -> Optional[ScData]: ...


@anndata_or_mudata_checker
def neighbors(
    scdata: ScData,  # type: ignore
    n_neighbors: int = 15,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    backend: NeighborBackend = "exact",
    metric: Metric = "euclidean",
    key_added: str = "neighbors",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    copy: bool = False,
) -> Optional[ScData]:  # type: ignore
    """
    Compute a k-nearest-neighbor graph of observations.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    n_neighbors: int (default: 15)
        Number of nearest neighbors.
    representation: str, optional (default: 'X_pca')
        Representation key in `scdata.obsm`.
    n_pcs: int, optional
        Number of representation dimensions to use. If None, use all
        dimensions in `representation`.
    backend: {'exact', 'pynndescent'} (default: 'exact')
        Neighbor-search backend. If `"exact"`, use scikit-learn brute-force
        exact nearest neighbors. If `"pynndescent"`, use approximate
        PyNNDescent neighbors.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    key_added: str (default: 'neighbors')
        Key used to store neighborhood graph metadata in `scdata.uns`.
    distances_key: str, optional
        Key used to store distances in `scdata.obsp`. Defaults to
        `"distances"` when `key_added="neighbors"` and
        `f"{key_added}_distances"` otherwise.
    connectivities_key: str, optional
        Key used to store connectivities in `scdata.obsp`. Defaults to
        `"connectivities"` when `key_added="neighbors"` and
        `f"{key_added}_connectivities"` otherwise.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by approximate neighbor-search
        backends.
    n_jobs: int (default: 1)
        Number of allocated processors.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=50)
    >>> bt.sct.tl.neighbors(adata, representation="X_pca")

    Returns
    -------
    ScData or None
        If `copy=True`, returns a copy of `scdata` with neighborhood graph
        results added. Otherwise, updates `scdata` in place and returns None.

        Neighborhood graph results are stored in:

        - `scdata.uns[key_added]`: neighborhood graph metadata;
        - `scdata.obsp[distances_key]`: pairwise neighbor distances;
        - `scdata.obsp[connectivities_key]`: neighborhood graph connectivities.

    Notes
    -----
    PyNNDescent is approximate. Its recall can decrease substantially as the
    embedding dimensionality increases. For PCA representations with many
    components (e.g. 50 PCs or more), the exact backend may be preferable even
    when PyNNDescent is slightly faster.
    """

    n_neighbors = _as_positive_integer(n_neighbors, "n_neighbors")
    if n_neighbors < 2:
        raise ValueError(
            f"invalid argument value for 'n_neighbors': "
            f"expected value >= 2 but received {n_neighbors!r}"
        )
    if n_neighbors > scdata.n_obs:
        raise ValueError(
            f"invalid argument value for 'n_neighbors': "
            f"expected value <= number of observations "
            f"({scdata.n_obs}) but received {n_neighbors!r}"
        )

    if n_pcs is not None:
        n_pcs = _as_positive_integer(n_pcs, "n_pcs")

    backend = _as_literal(
        backend,
        choices=get_args(NeighborBackend),
        name="backend",
    )
    representation = (
        "X_pca"
        if representation is None
        else _as_string(representation, "representation")
    )
    metric = _as_literal(
        metric,
        choices=get_args(Metric),
        name="metric",
    )
    key_added = _as_string(key_added, "key_added")
    _as_seed(seed)
    if not isinstance(n_jobs, int):
        raise TypeError(
            f"unsupported argument type for 'n_jobs': "
            f"expected {int} but received {type(n_jobs)}"
        )

    if distances_key is None:
        distances_key = (
            "distances" if key_added == "neighbors" else f"{key_added}_distances"
        )
    else:
        distances_key = _as_string(distances_key, "distances_key")

    if connectivities_key is None:
        connectivities_key = (
            "connectivities"
            if key_added == "neighbors"
            else f"{key_added}_connectivities"
        )
    else:
        connectivities_key = _as_string(connectivities_key, "connectivities_key")

    scdata = scdata.copy() if copy else scdata
    resolved_n_pcs = (
        n_pcs
        if n_pcs is not None
        else int(cast(Any, scdata.obsm[representation]).shape[1])
    )
    distances, connectivities = _kneighbors_graph_matrices(
        scdata=scdata,
        use_rep=representation,
        n_components=n_pcs,
        n_neighbors=n_neighbors,
        backend=backend,
        metric=metric,
        seed=seed,
        n_jobs=n_jobs,
    )

    scdata.obsp[distances_key] = distances
    scdata.obsp[connectivities_key] = connectivities
    scdata.uns[key_added] = {
        "distances_key": distances_key,
        "connectivities_key": connectivities_key,
        "params": {
            "n_neighbors": n_neighbors,
            "n_pcs": resolved_n_pcs,
            "representation": representation,
            "backend": backend,
            "metric": metric,
            "connectivity_method": "umap",
            "seed": _format_random_state(seed),
        },
    }

    return scdata if copy else None


@anndata_or_mudata_checker
def knn_graph(
    scdata: ScData,  # type: ignore
    use_rep: Optional[str] = None,
    n_components: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    metric: Metric = "euclidean",
    create_using: Type[Any] = DiGraph,
    edge_attr: str = "distance",
    index_or_name: IndexOrName = "index",
    n_jobs: int = 1,
    **metric_kwargs: Any,
) -> Graph[Any]:
    """
    Compute a k-nearest-neighbor graph from an embedding space.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    use_rep: str, optional
        Representation key in `scdata.obsm`.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.
    n_neighbors: int
        Number of nearest neighbors.
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

    if n_neighbors is None:
        raise TypeError("missing required argument: 'n_neighbors'")

    if metric_kwargs:
        from sklearn import neighbors as sklearn_neighbors

        representation_mtx = get_representation(
            scdata,
            use_rep=use_rep,
            n_components=n_components,
        )
        weighted_adjacency_matrix = sklearn_neighbors.kneighbors_graph(
            X=representation_mtx,
            n_neighbors=n_neighbors,
            mode="distance",
            metric=metric,
            n_jobs=n_jobs,
            **metric_kwargs,
        )
    else:
        weighted_adjacency_matrix = _kneighbors_distance_matrix(
            scdata=scdata,
            use_rep=use_rep,
            n_components=n_components,
            n_neighbors=n_neighbors,
            metric=metric,
            n_jobs=n_jobs,
        )
    try:
        graph = nx.from_scipy_sparse_array(
            weighted_adjacency_matrix,
            create_using=create_using,
            edge_attribute=edge_attr,
        )
    except AttributeError:
        from_scipy_sparse_matrix = getattr(nx, "from_scipy_sparse_matrix")
        graph = from_scipy_sparse_matrix(
            weighted_adjacency_matrix,
            create_using=create_using,
            edge_attribute=edge_attr,
        )

    index_or_name = _as_literal(
        index_or_name,
        choices=("index", "name"),
        name="index_or_name",
    )

    if index_or_name == "index":
        return graph
    if index_or_name == "name":
        return nx.relabel_nodes(
            graph,
            dict(zip(list(range(len(scdata.obs.index))), scdata.obs.index)),
        )


def kneighbors_graph(*args: Any, **kwargs: Any) -> Graph[Any]:
    """
    Deprecated alias for `knn_graph`.
    """

    _warn_deprecated(
        "`bt.sct.tl.kneighbors_graph`",
        replacement="`bt.sct.tl.knn_graph`",
        stacklevel=2,
    )
    return knn_graph(*args, **kwargs)


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
        _use_rep: str
        _n_components: Optional[int]
        _metric: Metric
        _method: Shortest_Path_Method

    def __init__(
        self,
        use_rep: Optional[str] = None,
        n_components: Optional[int] = None,
        n_neighbors: Optional[int] = None,
        metric: Optional[Metric] = None,
        **metric_kwargs: Any,
    ):
        self._deprecated_init_params: Dict[str, Any] = {}
        self._metric_kwargs: Dict[str, Any] = {}

        if (
            n_neighbors is not None
            or use_rep is not None
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
                resolved_use_rep,
                resolved_n_components,
                resolved_n_neighbors,
                resolved_metric,
                resolved_metric_kwargs,
            ) = _normalize_knnsc_configuration(
                use_rep=use_rep,
                n_components=n_components,
                n_neighbors=n_neighbors,
                metric=metric,
                metric_kwargs=metric_kwargs,
            )
            self._deprecated_init_params = {
                "use_rep": resolved_use_rep,
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
                f"use_rep={fit_params['use_rep']}, "
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
            f"use_rep={self._deprecated_init_params['use_rep']}, "
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
    def use_rep(self) -> str:
        """
        Representation key used during fitting.
        """

        return cast(str, self._require_fitted_attribute("_use_rep"))

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
            "use_rep": self.use_rep,
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
        use_rep: Optional[str] = None,
        n_components: Optional[int] = None,
        n_neighbors: Optional[int] = None,
        metric: Optional[Metric] = None,
        metric_kwargs: Optional[Dict[str, Any]] = None,
        n_jobs: int = 1,
        method: Shortest_Path_Method = "dijkstra",
        *,
        obs: Optional[str] = None,
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
        use_rep: str, optional (default: "X_pca")
            Representation key in `adata.obsm`.
        n_components: int, optional
            Number of dimensions to use. If None, use all dimensions.
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

        init_params = self._deprecated_init_params
        if n_neighbors is None and "n_neighbors" in init_params:
            n_neighbors = cast(int, init_params["n_neighbors"])
        if use_rep is None:
            use_rep = cast(Optional[str], init_params.get("use_rep"))
        if n_components is None:
            n_components = cast(Optional[int], init_params.get("n_components"))
        if metric is None:
            metric = cast(Optional[Metric], init_params.get("metric"))
        if metric_kwargs is None:
            metric_kwargs = cast(
                Optional[Dict[str, Any]], init_params.get("metric_kwargs")
            )

        (
            resolved_use_rep,
            resolved_n_components,
            resolved_n_neighbors,
            resolved_metric,
            resolved_metric_kwargs,
        ) = _normalize_knnsc_configuration(
            use_rep=use_rep,
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
            adata, use_rep=resolved_use_rep, n_components=resolved_n_components
        )

        _knn_graph = knn_graph(
            adata,
            use_rep=resolved_use_rep,
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
                    use_rep=resolved_use_rep,
                    n_components=resolved_n_components,
                )
                y_representation_mtx = get_representation(
                    adata_any[paired_scc[1], :],
                    use_rep=resolved_use_rep,
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
        self._use_rep = resolved_use_rep
        self._n_components = resolved_n_components
        self._n_neighbors = resolved_n_neighbors
        self._metric = resolved_metric
        self._metric_kwargs = resolved_metric_kwargs
        self._method = method
        self._compute_shortest_path_lengths(method=method, n_jobs=n_jobs)
        return None

    def _compute_shortest_path_lengths(
        self, method: Shortest_Path_Method = "dijkstra", n_jobs: int = 1
    ) -> None:

        def shortest_path_lengths_from(
            source: str,
        ):
            distances = pd.Series(data=self.obs.index, index=self.obs.index)
            distances = distances.apply(
                lambda x: nx.shortest_path_length(
                    self.knn_graph,
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

        self._shortest_path_lengths_df = pd.DataFrame.from_dict(
            data=shortest_path_lengths_dict, orient="columns"
        )
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
            subclusters_series.loc[_obs] = cluster

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
            subclusters_series.loc[_obs] = cluster

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
        for cluster in invalid_clusters:
            if cluster in self.cluster_counts.index:
                size = int(self.cluster_counts[cluster])
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
        use_rep: str = "X_pca",
        n_components: Optional[int] = None,
        metric: Metric = "euclidean",
        **metric_kwargs: Any,
    ):
        _warn_deprecated(
            "`bt.sct.tl.Knnbs`",
            replacement="`bt.sct.tl.KNNSC`",
            stacklevel=2,
        )
        super().__init__(
            n_neighbors=n_neighbors,
            use_rep=use_rep,
            n_components=n_components,
            metric=metric,
            **metric_kwargs,
        )


@anndata_or_mudata_checker
def _shared_nearest_neighbors_graph(
    scdata: ScData, cluster_key: str, prune_snn: float  # type: ignore
) -> csr_matrix:

    n_neighbors = scdata.uns[cluster_key]["params"]["n_neighbors"] - 1
    prune_snn = _as_non_negative_number(prune_snn, "prune_snn")
    if prune_snn < 1:
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


@overload
def shared_neighbors(
    scdata: ScData,
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[Union[float, int]] = 1 / 15,
    normalize_connectivities: bool = True,
    metric: Metric = "euclidean",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def shared_neighbors(
    scdata: ScData,
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[Union[float, int]] = 1 / 15,
    normalize_connectivities: bool = True,
    metric: Metric = "euclidean",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    *,
    copy: Literal[True],
) -> ScData: ...


@overload
def shared_neighbors(
    scdata: ScData,
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[Union[float, int]] = 1 / 15,
    normalize_connectivities: bool = True,
    metric: Metric = "euclidean",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False,
) -> Optional[ScData]: ...


@anndata_or_mudata_checker
def shared_neighbors(
    scdata: ScData,  # type: ignore
    knn_key: str = "neighbors",
    snn_key: str = "shared_neighbors",
    prune_snn: Optional[Union[float, int]] = 1 / 15,
    normalize_connectivities: bool = True,
    metric: Metric = "euclidean",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False,
) -> Optional[ScData]:  # type: ignore
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
    normalize_connectivities: bool (default: True)
        If False, connectivities store the absolute number of shared neighbors.
        Otherwise, connectivities are normalized.
    metric: Metric (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    distances_key: str (optional, default: None)
        Key used to store distances in `scdata.obsp`.
    connectivities_key: str (optional, default: None)
        Key used to store connectivities in `scdata.obsp`.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=50)
    >>> bt.sct.tl.neighbors(adata, representation="X_pca")
    >>> bt.sct.tl.shared_neighbors(adata)

    Returns
    -------
    ScData or None
        If `copy=True`, returns a copy of `scdata` with shared-neighbor
        results added. Otherwise, updates `scdata` in place and returns None.

        Shared-neighbor results are stored in:

        - `scdata.uns[snn_key]`: shared-neighbor graph metadata;
        - `scdata.obsp[distances_key]`: pairwise distances;
        - `scdata.obsp[connectivities_key]`: shared-neighbor connectivities.

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

    knn_params = scdata.uns[knn_key]["params"]
    n_pcs = knn_params["n_pcs"]
    obsm = knn_params.get("representation", knn_params.get("use_rep"))
    obsm = _as_string(obsm, "representation")

    representation_mtx = scdata.obsm[obsm][:, 0:n_pcs]
    zeros_ones = snn_graph.toarray()
    zeros_ones[zeros_ones > 0] = 1

    distances_matrix = pairwise_distances(representation_mtx, metric=metric)
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
