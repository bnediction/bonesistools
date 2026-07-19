#!/usr/bin/env python

from __future__ import annotations

import math
from functools import wraps
from typing import (
    Any,
    Callable,
    Dict,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
    overload,
)

import networkx as nx
import numpy as np
from networkx import DiGraph, Graph
from scipy.sparse import (
    coo_matrix,
    csr_matrix,
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
    Matrix,
    Metric,
    ScData,
    SNNMetric,
    anndata_or_mudata_checker,
)
from ._utils import (
    _UNSET,
    _resolve_representation_argument,
    get_representation,
)

IndexOrName = Literal["index", "name"]
NeighborBackend = Literal["exact", "pynndescent"]
NeighborConnectivityMethod = Literal["fuzzy", "binary"]
_F = TypeVar("_F", bound=Callable[..., Any])
_SMOOTH_K_TOLERANCE = 1e-5
_MIN_K_DIST_SCALE = 0.001
_NPY_FLOATMAX = np.float32(3.4028235e38)
_SHARED_SCORE_CHUNK_SIZE = 500_000


@overload
def neighbors(
    scdata: ScData,
    n_neighbors: int = 15,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    backend: NeighborBackend = "exact",
    connectivity_method: NeighborConnectivityMethod = "fuzzy",
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
    connectivity_method: NeighborConnectivityMethod = "fuzzy",
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
    connectivity_method: NeighborConnectivityMethod = "fuzzy",
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
    connectivity_method: NeighborConnectivityMethod = "fuzzy",
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
        Number of representation dimensions to use. If `None`, use all
        dimensions in `representation`.
    backend: {'exact', 'pynndescent'} (default: 'exact')
        Neighbor-search backend. If `"exact"`, use scikit-learn brute-force
        exact nearest neighbors. If `"pynndescent"`, use approximate
        PyNNDescent neighbors.
    connectivity_method: {'fuzzy', 'binary'} (default: 'fuzzy')
        Method used to convert nearest neighbors into graph connectivities. If
        `"fuzzy"`, use UMAP fuzzy simplicial-set connectivities. If
        `"binary"`, use a binary symmetrized k-nearest-neighbor graph.
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
    >>> bt.omics.tl.pca(adata, n_components=50)
    >>> bt.omics.tl.neighbors(adata, representation="X_pca")

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
    connectivity_method = _as_literal(
        connectivity_method,
        choices=get_args(NeighborConnectivityMethod),
        name="connectivity_method",
    )
    if connectivity_method == "binary" and backend != "exact":
        raise ValueError(
            "invalid argument combination: "
            "connectivity_method='binary' requires backend='exact'"
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
        representation=representation,
        n_components=n_pcs,
        n_neighbors=n_neighbors,
        backend=backend,
        connectivity_method=connectivity_method,
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
            "connectivity_method": connectivity_method,
            "seed": _format_random_state(seed),
        },
    }

    return scdata if copy else None


@anndata_or_mudata_checker
def knn_graph(
    scdata: ScData,  # type: ignore
    representation: Any = _UNSET,
    n_components: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    metric: Metric = "euclidean",
    create_using: Type[Any] = DiGraph,
    edge_attr: str = "distance",
    index_or_name: IndexOrName = "index",
    n_jobs: int = 1,
    use_rep: Any = _UNSET,
    **metric_kwargs: Any,
) -> Graph[Any]:
    """
    Compute a k-nearest-neighbor graph from an embedding space.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    representation: str, optional
        Representation key in `scdata.obsm`.
    n_components: int, optional
        Number of dimensions to use. If `None`, use all dimensions.
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
    use_rep: str, optional
        Deprecated alias for `representation`.

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

    representation = _resolve_representation_argument(
        representation,
        use_rep,
        default=None,
        stacklevel=2,
    )

    if metric_kwargs:
        from sklearn import neighbors as sklearn_neighbors

        representation_mtx = get_representation(
            scdata,
            obsm=representation,
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
            representation=representation,
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
        "`bt.omics.tl.kneighbors_graph`",
        replacement="`bt.omics.tl.knn_graph`",
        stacklevel=2,
    )
    return knn_graph(*args, **kwargs)


def _kneighbors_graph_matrices(
    scdata: ScData,
    representation: Optional[str],
    n_components: Optional[int],
    n_neighbors: int,
    backend: NeighborBackend,
    connectivity_method: NeighborConnectivityMethod,
    metric: Metric,
    seed: RandomStateSeed,
    n_jobs: int,
) -> Tuple[csr_matrix, csr_matrix]:

    representation_mtx = cast(
        Matrix,
        get_representation(
            scdata,
            obsm=representation,
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
    if connectivity_method == "binary":
        from sklearn import neighbors as sklearn_neighbors

        connectivities = cast(
            csr_matrix,
            sklearn_neighbors.kneighbors_graph(
                X=representation_mtx,
                n_neighbors=n_neighbors,
                mode="connectivity",
                metric=metric,
                include_self=True,
                n_jobs=n_jobs,
            ),
        )
        connectivities = cast(csr_matrix, 0.5 * (connectivities + connectivities.T))
        return distances, connectivities

    connectivities = _umap_connectivities_from_knn(
        knn_indices,
        knn_distances,
        n_obs=scdata.n_obs,
    )
    return distances, connectivities


def _kneighbors_distance_matrix(
    scdata: ScData,
    representation: Optional[str],
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
            obsm=representation,
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
    return _sort_knn_arrays(
        knn_indices=knn_indices.astype(np.int32, copy=False),
        knn_distances=knn_distances.astype(np.float32, copy=False),
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


def _umap_connectivities_from_knn(
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
    n_obs: int,
) -> csr_matrix:

    knn_indices = knn_indices.astype(np.int32, copy=False)
    knn_distances = knn_distances.astype(np.float32, copy=False)
    knn_distances = _zero_self_distances(
        knn_indices=knn_indices,
        knn_distances=knn_distances,
    )
    n_neighbors = knn_indices.shape[1]
    sigmas, rhos = _umap_smooth_knn_distances(knn_distances, n_neighbors)

    rows = np.repeat(np.arange(n_obs, dtype=np.int32), n_neighbors)
    columns = knn_indices.ravel()
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
            return _sort_knn_arrays(
                knn_indices=indices[:, :n_neighbors].astype(np.int32, copy=False),
                knn_distances=distances[:, :n_neighbors].astype(
                    np.float32,
                    copy=False,
                ),
            )

    return _knn_arrays_from_sparse_distances_slow(
        sparse_distances=sparse_distances,
        n_neighbors=n_neighbors,
    )


def _sort_knn_arrays(
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:

    knn_distances = _zero_self_distances(
        knn_indices=knn_indices,
        knn_distances=knn_distances,
    )
    rows = np.arange(knn_indices.shape[0], dtype=knn_indices.dtype)[:, None]
    self_order = np.where(knn_indices == rows, 0, 1)
    order = np.lexsort(
        (
            knn_indices,
            knn_distances,
            self_order,
        ),
        axis=1,
    )
    return (
        np.take_along_axis(knn_indices, order, axis=1),
        np.take_along_axis(knn_distances, order, axis=1),
    )


def _zero_self_distances(
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
) -> np.ndarray:
    """Set numerical self-distances to exact zeros."""

    rows = np.arange(knn_indices.shape[0], dtype=knn_indices.dtype)[:, None]
    non_zero_self = (knn_indices == rows) & (knn_distances != 0.0)
    if not np.any(non_zero_self):
        return knn_distances

    knn_distances = knn_distances.copy()
    knn_distances[non_zero_self] = 0.0
    return knn_distances


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

    return _sort_knn_arrays(knn_indices=knn_indices, knn_distances=knn_distances)


def _support_legacy_shared_neighbors_arguments(function: _F) -> _F:
    """Accept deprecated shared-neighbor argument names."""

    deprecated_arguments = {
        "knn_key": "neighbors_key",
        "snn_key": "key_added",
    }

    @wraps(function)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        for old_name, new_name in deprecated_arguments.items():
            if old_name not in kwargs:
                continue
            _warn_deprecated_argument(old_name, new_name, stacklevel=4)
            if new_name in kwargs:
                raise TypeError(
                    f"invalid argument combination: use either '{old_name}' "
                    f"or '{new_name}', not both"
                )
            kwargs[new_name] = kwargs.pop(old_name)
        return function(*args, **kwargs)

    return cast(_F, wrapper)


@overload
def shared_neighbors(
    scdata: ScData,
    *,
    neighbors_key: str = "neighbors",
    key_added: str = "shared_neighbors",
    prune: Optional[Union[float, int]] = 1 / 15,
    metric: SNNMetric = "jaccard",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def shared_neighbors(
    scdata: ScData,
    *,
    neighbors_key: str = "neighbors",
    key_added: str = "shared_neighbors",
    prune: Optional[Union[float, int]] = 1 / 15,
    metric: SNNMetric = "jaccard",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: Literal[True],
) -> ScData: ...


@overload
def shared_neighbors(
    scdata: ScData,
    *,
    neighbors_key: str = "neighbors",
    key_added: str = "shared_neighbors",
    prune: Optional[Union[float, int]] = 1 / 15,
    metric: SNNMetric = "jaccard",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False,
) -> Optional[ScData]: ...


@_support_legacy_shared_neighbors_arguments
@anndata_or_mudata_checker
def shared_neighbors(
    scdata: ScData,  # type: ignore
    *,
    neighbors_key: str = "neighbors",
    key_added: str = "shared_neighbors",
    prune: Optional[Union[float, int]] = 1 / 15,
    metric: SNNMetric = "jaccard",
    distances_key: Optional[str] = None,
    connectivities_key: Optional[str] = None,
    copy: bool = False,
) -> Optional[ScData]:  # type: ignore
    """
    Compute a shared nearest neighbor (SNN) graph of observations.

    Shared-neighbor connectivities are computed between the neighborhoods
    stored by a previous call to `neighbors` or an equivalent tool.

    Examples
    --------
    >>> bt.omics.tl.pca(adata, n_components=50)
    >>> bt.omics.tl.neighbors(adata, representation="X_pca")
    >>> bt.omics.tl.shared_neighbors(adata)

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    neighbors_key: str (default: 'neighbors')
        Key in `scdata.uns` containing the source neighborhood graph metadata.
    key_added: str (default: 'shared_neighbors')
        Key used to store shared-neighbor graph metadata in `scdata.uns`.
    prune: float | int (default: 1/15)
        Prune edges according to the raw number of shared neighbors before
        computing the selected similarity metric. If zero, no pruning is
        performed. Otherwise, edges at or below the threshold are removed.
        The threshold can be relative (float between 0 and 1) or absolute
        (integer between 1 and k).
    metric: {'jaccard', 'weighted-jaccard', 'overlap', 'binary-cosine',
             'ranked-jaccard'}
        Shared-neighbor similarity metric:
        - `jaccard` (default): intersection divided by union;
        - `weighted-jaccard`: source-connectivity-weighted Jaccard;
        - `overlap`: intersection divided by the smaller neighborhood;
        - `binary-cosine`: cosine similarity of binary neighborhood vectors;
        - `ranked-jaccard`: rank-weighted Jaccard normalized to `[0, 1]`.
    distances_key: str (optional, default: None)
        Key used to store distances in `scdata.obsp`. If `None`, uses
        `f"{key_added}_distances"`.
    connectivities_key: str (optional, default: None)
        Key used to store connectivities in `scdata.obsp`. If `None`, uses
        `f"{key_added}_connectivities"`.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Notes
    -----
    Unlike `neighbors`, this function constructs a shared nearest neighbor
    (SNN) graph whose edge weights measure the similarity between local
    neighborhoods rather than pairwise proximity. All supported similarity
    metrics return connectivities normalized to the interval `[0, 1]`,
    allowing the resulting graph to be used directly with graph-based
    algorithms such as Leiden, Louvain, PAGA, diffusion maps and UMAP.

    Distances are defined as

        distance = 1 - connectivity

    for every retained edge.

    Returns
    -------
    ScData or None
        If `copy=True`, returns a copy of `scdata` with shared-neighbor
        results added. Otherwise, updates `scdata` in place and returns None.

        Shared-neighbor results are stored in:

        - `scdata.uns[key_added]`: shared-neighbor graph metadata;
        - `scdata.obsp[distances_key]`: pairwise distances;
        - `scdata.obsp[connectivities_key]`: shared-neighbor connectivities.

    Raises
    ------
    KeyError
        If `neighbors_key` is not found in `scdata.uns`.
    ValueError
        If `prune` is negative or greater than or equal to the number of
        neighbors, or if `metric` is invalid.
    """

    if neighbors_key not in scdata.uns:
        raise KeyError(
            "neighborhood graph not found in 'scdata': "
            "please run 'scanpy.pp.neighbors' or specify 'neighbors_key'"
        )
    if prune is None:
        prune = 0
    if distances_key is None:
        distances_key = f"{key_added}_distances"
    if connectivities_key is None:
        connectivities_key = f"{key_added}_connectivities"
    metric = _as_literal(
        metric,
        choices=get_args(SNNMetric),
        name="metric",
    )
    n_neighbors = _as_positive_integer(
        scdata.uns[neighbors_key]["params"]["n_neighbors"],
        "n_neighbors",
    )
    neighborhood_size = n_neighbors - 1
    prune_threshold = _resolve_shared_neighbor_pruning(
        prune,
        neighborhood_size=neighborhood_size,
    )

    scdata = scdata.copy() if copy else scdata

    source_distances = _get_neighbor_distances(
        scdata,
        neighbors_key=neighbors_key,
    )
    neighborhood_graph = _binary_neighborhood_graph(source_distances)
    shared_counts = _shared_neighbor_counts(
        neighborhood_graph,
        prune_threshold=prune_threshold,
    )
    degrees = np.diff(neighborhood_graph.indptr)
    connectivities_matrix = _shared_connectivities(
        scdata,
        neighbors_key=neighbors_key,
        source_distances=source_distances,
        neighborhood_graph=neighborhood_graph,
        shared_counts=shared_counts,
        degrees=degrees,
        metric=metric,
    )
    distances_matrix = connectivities_matrix.copy()
    distances_matrix.data = np.asarray(
        1.0 - distances_matrix.data,
        dtype=np.float32,
    )

    scdata.obsp[distances_key] = distances_matrix
    scdata.obsp[connectivities_key] = connectivities_matrix

    scdata.uns[key_added] = dict()
    scdata.uns[key_added]["distances_key"] = distances_key
    scdata.uns[key_added]["connectivities_key"] = connectivities_key
    scdata.uns[key_added]["params"] = {
        "knn_base": f"scdata.uns['{neighbors_key}']",
        "prune": prune_threshold,
        "metric": metric,
    }

    return scdata if copy else None


def _resolve_shared_neighbor_pruning(
    prune: Optional[Union[float, int]],
    *,
    neighborhood_size: int,
) -> int:
    """Resolve a relative or absolute shared-neighbor threshold."""

    resolved = _as_non_negative_number(
        0 if prune is None else prune,
        "prune",
    )
    threshold = (
        math.ceil(neighborhood_size * resolved)
        if resolved < 1
        else _as_non_negative_integer(resolved, "prune")
    )
    if threshold >= neighborhood_size:
        raise ValueError(
            f"invalid argument values for 'prune' and 'n_neighbors': "
            f"expected prune < n_neighbors but received "
            f"prune={prune!r} and n_neighbors={neighborhood_size!r}"
        )
    return threshold


def _get_neighbor_distances(
    scdata: ScData,
    *,
    neighbors_key: str,
) -> csr_matrix:
    """Retrieve and validate a CSR copy of KNN distances without self-loops."""

    distances_key = scdata.uns[neighbors_key]["distances_key"]
    source = scdata.obsp[distances_key]
    if issparse(source):
        distances = cast(csr_matrix, source.tocsr(copy=True))
    else:
        distances = csr_matrix(source)
    distances.sum_duplicates()
    _remove_csr_self_loops(distances)
    distances.sort_indices()
    if np.any(~np.isfinite(distances.data)) or np.any(distances.data < 0):
        raise ValueError("invalid KNN distance graph: expected finite distances >= 0")
    return distances


def _binary_neighborhood_graph(source_distances: csr_matrix) -> csr_matrix:
    """Encode non-zero KNN entries as binary neighborhood vectors."""

    return csr_matrix(
        (
            np.ones(source_distances.nnz, dtype=np.int32),
            source_distances.indices.copy(),
            source_distances.indptr.copy(),
        ),
        shape=source_distances.shape,
    )


def _shared_neighbor_counts(
    neighborhood_graph: csr_matrix,
    *,
    prune_threshold: int,
) -> csr_matrix:
    """Count shared neighbors for every pair with a non-empty intersection."""

    counts = cast(
        csr_matrix,
        neighborhood_graph @ neighborhood_graph.transpose(),
    )
    counts = counts.tocsr()
    if np.all(np.diff(neighborhood_graph.indptr) > 0):
        counts.setdiag(0)
        counts.eliminate_zeros()
    else:
        _remove_csr_self_loops(counts)
    if prune_threshold:
        counts.data[counts.data <= prune_threshold] = 0
        counts.eliminate_zeros()
    counts.sort_indices()
    return counts


def _shared_connectivities(
    scdata: ScData,
    *,
    neighbors_key: str,
    source_distances: csr_matrix,
    neighborhood_graph: csr_matrix,
    shared_counts: csr_matrix,
    degrees: np.ndarray,
    metric: SNNMetric,
) -> csr_matrix:
    """Normalize shared-neighbor scores with the selected metric."""

    if metric in ("jaccard", "overlap", "binary-cosine"):
        return _binary_shared_connectivities(
            shared_counts,
            degrees=degrees,
            metric=metric,
        )
    if metric == "weighted-jaccard":
        weights = _shared_neighborhood_weights(
            scdata,
            neighbors_key=neighbors_key,
            neighborhood_graph=neighborhood_graph,
        )
        return _weighted_jaccard_connectivities(
            weights,
            shared_counts=shared_counts,
        )
    return _ranked_jaccard_connectivities(
        source_distances,
        shared_counts=shared_counts,
        degrees=degrees,
    )


def _binary_shared_connectivities(
    shared_counts: csr_matrix,
    *,
    degrees: np.ndarray,
    metric: SNNMetric,
) -> csr_matrix:
    """Normalize binary-neighborhood intersections."""

    rows = _csr_rows(shared_counts)
    columns = shared_counts.indices
    intersections = shared_counts.data.astype(np.float64, copy=False)
    left_degrees = degrees[rows].astype(np.float64, copy=False)
    right_degrees = degrees[columns].astype(np.float64, copy=False)
    if metric == "jaccard":
        denominators = left_degrees + right_degrees - intersections
    elif metric == "overlap":
        denominators = np.minimum(left_degrees, right_degrees)
    else:
        denominators = np.sqrt(left_degrees * right_degrees)
    values = np.divide(
        intersections,
        denominators,
        out=np.zeros_like(intersections),
        where=denominators > 0,
    ).astype(np.float32, copy=False)
    return _connectivity_matrix(shared_counts, values)


def _shared_neighborhood_weights(
    scdata: ScData,
    *,
    neighbors_key: str,
    neighborhood_graph: csr_matrix,
) -> csr_matrix:
    """Restrict source KNN connectivities to directed neighborhoods."""

    connectivities_key = scdata.uns[neighbors_key].get("connectivities_key")
    if connectivities_key is None or connectivities_key not in scdata.obsp:
        raise KeyError(
            "KNN connectivities not found in 'scdata': "
            "metric='weighted-jaccard' requires a source connectivity graph"
        )
    source = scdata.obsp[connectivities_key]
    if issparse(source):
        weights = cast(csr_matrix, source.tocsr(copy=True))
    else:
        weights = csr_matrix(source)
    weights.sum_duplicates()
    _remove_csr_self_loops(weights)
    if np.any(~np.isfinite(weights.data)) or np.any(weights.data < 0):
        raise ValueError("invalid KNN connectivity graph: expected finite weights >= 0")
    weights = cast(csr_matrix, weights.multiply(neighborhood_graph).tocsr())
    weights.eliminate_zeros()
    weights.sort_indices()
    return weights.astype(np.float32, copy=False)


def _weighted_jaccard_connectivities(
    weights: csr_matrix,
    *,
    shared_counts: csr_matrix,
) -> csr_matrix:
    """Normalize minimum weighted intersections by weighted unions."""

    intersections = _pairwise_shared_scores(weights, method="minimum")
    intersections = _restrict_shared_scores(intersections, shared_counts)
    rows = _csr_rows(intersections)
    columns = intersections.indices
    weight_sums = np.asarray(weights.sum(axis=1)).ravel()
    denominators = weight_sums[rows] + weight_sums[columns] - intersections.data
    values = np.divide(
        intersections.data,
        denominators,
        out=np.zeros_like(intersections.data),
        where=denominators > 0,
    ).astype(np.float32, copy=False)
    return _connectivity_matrix(intersections, values)


def _ranked_jaccard_connectivities(
    source_distances: csr_matrix,
    *,
    shared_counts: csr_matrix,
    degrees: np.ndarray,
) -> csr_matrix:
    """Normalize rank-weighted shared-neighbor scores."""

    ranks = _shared_neighbor_rank_matrix(source_distances)
    scores = _pairwise_shared_scores(ranks, method="ranked")
    scores = _restrict_shared_scores(scores, shared_counts)
    max_degree = int(degrees.max(initial=0))
    if max_degree == 0:
        return csr_matrix(shared_counts.shape, dtype=np.float32)
    maximum_score = np.sum(1.0 / (2.0 * np.arange(1, max_degree + 1, dtype=np.float64)))
    values = np.asarray(scores.data / maximum_score, dtype=np.float32)
    return _connectivity_matrix(scores, values)


def _shared_neighbor_rank_matrix(source_distances: csr_matrix) -> csr_matrix:
    """Encode each neighbor by its one-based distance rank."""

    ranks = _binary_neighborhood_graph(source_distances).astype(np.float32)
    n_rows = cast(Tuple[int, int], source_distances.shape)[0]
    for row in range(n_rows):
        start = source_distances.indptr[row]
        end = source_distances.indptr[row + 1]
        if end == start:
            continue
        order = np.lexsort(
            (
                source_distances.indices[start:end],
                source_distances.data[start:end],
            )
        )
        row_ranks = np.empty(end - start, dtype=np.float32)
        row_ranks[order] = np.arange(1, end - start + 1, dtype=np.float32)
        ranks.data[start:end] = row_ranks
    return ranks


def _pairwise_shared_scores(
    values: csr_matrix,
    *,
    method: Literal["minimum", "ranked"],
) -> csr_matrix:
    """Accumulate weighted contributions over shared sparse columns."""

    columns = values.tocsc()
    shape = cast(Tuple[int, int], values.shape)
    scores = csr_matrix(shape, dtype=np.float32)
    row_chunks = []
    column_chunks = []
    value_chunks = []
    pair_indices: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    pending_scores = 0
    n_columns = cast(Tuple[int, int], columns.shape)[1]
    for feature in range(n_columns):
        start = columns.indptr[feature]
        end = columns.indptr[feature + 1]
        size = end - start
        if size < 2:
            continue
        pairs = pair_indices.get(size)
        if pairs is None:
            pairs = np.triu_indices(size, k=1)
            pair_indices[size] = pairs
        left, right = pairs
        rows = columns.indices[start:end]
        feature_values = columns.data[start:end]
        row_chunks.append(rows[left])
        column_chunks.append(rows[right])
        if method == "minimum":
            contributions = np.minimum(
                feature_values[left],
                feature_values[right],
            )
        else:
            contributions = 1.0 / (feature_values[left] + feature_values[right])
        value_chunks.append(np.asarray(contributions, dtype=np.float32))
        pending_scores += left.size
        if pending_scores >= _SHARED_SCORE_CHUNK_SIZE:
            scores = cast(
                csr_matrix,
                (
                    scores
                    + _symmetric_shared_score_chunk(
                        row_chunks,
                        column_chunks,
                        value_chunks,
                        shape=shape,
                    )
                ).tocsr(),
            )
            row_chunks.clear()
            column_chunks.clear()
            value_chunks.clear()
            pending_scores = 0
    if row_chunks:
        scores = cast(
            csr_matrix,
            (
                scores
                + _symmetric_shared_score_chunk(
                    row_chunks,
                    column_chunks,
                    value_chunks,
                    shape=shape,
                )
            ).tocsr(),
        )
    scores.sum_duplicates()
    scores.sort_indices()
    return scores


def _symmetric_shared_score_chunk(
    row_chunks: Sequence[np.ndarray],
    column_chunks: Sequence[np.ndarray],
    value_chunks: Sequence[np.ndarray],
    *,
    shape: Tuple[int, int],
) -> csr_matrix:
    """Combine sparse upper-triangle contributions and their transpose."""

    upper_rows = np.concatenate(row_chunks)
    upper_columns = np.concatenate(column_chunks)
    contributions = np.concatenate(value_chunks)
    return cast(
        csr_matrix,
        coo_matrix(
            (
                np.concatenate((contributions, contributions)),
                (
                    np.concatenate((upper_rows, upper_columns)),
                    np.concatenate((upper_columns, upper_rows)),
                ),
            ),
            shape=shape,
            dtype=np.float32,
        ).tocsr(),
    )


def _restrict_shared_scores(
    scores: csr_matrix,
    shared_counts: csr_matrix,
) -> csr_matrix:
    """Restrict weighted scores to retained shared-neighbor edges."""

    retained = shared_counts.copy()
    retained.data.fill(1)
    scores = cast(csr_matrix, scores.multiply(retained).tocsr())
    scores.eliminate_zeros()
    scores.sort_indices()
    return scores


def _connectivity_matrix(
    structure: csr_matrix,
    values: np.ndarray,
) -> csr_matrix:
    """Build a bounded connectivity matrix with an existing CSR structure."""

    values = np.asarray(values, dtype=np.float32)
    np.clip(values, 0.0, 1.0, out=values)
    connectivities = csr_matrix(
        (
            values,
            structure.indices.copy(),
            structure.indptr.copy(),
        ),
        shape=structure.shape,
    )
    connectivities.eliminate_zeros()
    return connectivities


def _csr_rows(matrix: csr_matrix) -> np.ndarray:
    """Return the row index associated with each CSR data entry."""

    n_rows = cast(Tuple[int, int], matrix.shape)[0]
    return np.repeat(
        np.arange(n_rows, dtype=matrix.indices.dtype),
        np.diff(matrix.indptr),
    )


def _remove_csr_self_loops(matrix: csr_matrix) -> None:
    """Remove existing diagonal entries without inserting new CSR elements."""

    if np.any(matrix.diagonal() != 0):
        diagonal = matrix.indices == _csr_rows(matrix)
        matrix.data[diagonal] = 0
    matrix.eliminate_zeros()
