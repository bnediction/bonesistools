#!/usr/bin/env python

from __future__ import annotations

import numbers
import random
from importlib import import_module
from typing import Any, Dict, Optional, Tuple, cast, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ..._compat import Literal
from ..._typing import AutoInteger, RandomStateSeed
from ..._validation import (
    _as_boolean,
    _as_literal,
    _as_positive_integer,
    _as_positive_number,
    _as_seed,
    _as_string,
)
from .._metadata import _format_random_state
from .._typing import anndata_checker
from ._utils import get_pairwise, get_representation


@overload
def kmeans(
    adata: AnnData,
    n_clusters: int = 8,
    *,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    n_init: int = 10,
    key_added: str = "kmeans",
    seed: RandomStateSeed = 0,
    copy: Literal[False] = False,
    **kwargs: Any,
) -> None: ...


@overload
def kmeans(
    adata: AnnData,
    n_clusters: int = 8,
    *,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    n_init: int = 10,
    key_added: str = "kmeans",
    seed: RandomStateSeed = 0,
    copy: Literal[True],
    **kwargs: Any,
) -> AnnData: ...


@overload
def kmeans(
    adata: AnnData,
    n_clusters: int = 8,
    *,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    n_init: int = 10,
    key_added: str = "kmeans",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]: ...


@anndata_checker
def kmeans(
    adata: AnnData,
    n_clusters: int = 8,
    *,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    n_init: int = 10,
    key_added: str = "kmeans",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]:
    """
    Cluster observations with k-means.

    K-means partitions observations into a fixed number of clusters by
    minimizing within-cluster squared Euclidean distances in an observation
    representation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    n_clusters: int (default: 8)
        Number of clusters to compute.
    representation: str, optional (default: 'X_pca')
        Representation key in `adata.obsm`.
    n_pcs: int, optional
        Number of representation dimensions to use. If None, use all
        dimensions from `representation`.
    n_init: int (default: 10)
        Number of independent clustering initializations. The solution with
        the lowest inertia is retained.
    key_added: str (default: 'kmeans')
        Column in `adata.obs` where cluster labels are stored.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by k-means.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    **kwargs: Any
        Additional keyword arguments passed to `sklearn.cluster.KMeans`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=50)
    >>> bt.sct.tl.kmeans(adata, n_clusters=8, representation="X_pca")

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with k-means clustering
        results added. Otherwise, updates `adata` in place and returns None.

        K-means clustering results are stored in:

        - `adata.obs[key_added]`: categorical cluster labels;
        - `adata.uns[key_added]`: clustering metadata.
    """

    from sklearn.cluster import KMeans

    KMeansClass = cast(Any, KMeans)

    n_clusters = _as_positive_integer(n_clusters, "n_clusters")
    if n_pcs is not None:
        n_pcs = _as_positive_integer(n_pcs, "n_pcs")

    n_init = _as_positive_integer(n_init, "n_init")
    key_added = _as_string(key_added, "key_added")
    resolved_random_state = _as_seed(seed)

    if n_clusters > adata.n_obs:
        raise ValueError(
            f"invalid argument value for 'n_clusters': "
            f"expected value smaller than or equal to number of observations "
            f"({adata.n_obs}) but received {n_clusters!r}"
        )

    adata = adata.copy() if copy else adata
    representation_mtx = get_representation(
        adata,
        obsm=representation,
        n_components=n_pcs,
    )
    if sparse.issparse(representation_mtx):
        representation_mtx = cast(Any, representation_mtx).tocsr()
    else:
        representation_mtx = np.asarray(representation_mtx)

    if representation_mtx.ndim != 2:
        raise ValueError(
            "invalid representation matrix: expected a two-dimensional matrix"
        )

    labels = cast(
        np.ndarray,
        KMeansClass(
            n_clusters=n_clusters,
            n_init=n_init,
            random_state=resolved_random_state,
            **kwargs,
        ).fit_predict(representation_mtx),
    )

    labels = pd.Categorical([str(cluster) for cluster in labels])
    adata.obs[key_added] = pd.Series(labels, index=adata.obs_names, name=key_added)
    adata.uns[key_added] = {
        "params": {
            "method": "kmeans",
            "n_clusters": n_clusters,
            "representation": representation,
            "n_pcs": n_pcs,
            "n_init": n_init,
            "seed": _format_random_state(seed),
            **kwargs,
        }
    }

    return adata if copy else None


@overload
def louvain(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    weighted: bool = True,
    key_added: str = "louvain",
    seed: RandomStateSeed = 0,
    copy: Literal[False] = False,
    **kwargs: Any,
) -> None: ...


@overload
def louvain(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    weighted: bool = True,
    key_added: str = "louvain",
    seed: RandomStateSeed = 0,
    copy: Literal[True],
    **kwargs: Any,
) -> AnnData: ...


@overload
def louvain(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    weighted: bool = True,
    key_added: str = "louvain",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]: ...


@anndata_checker
def louvain(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    weighted: bool = True,
    key_added: str = "louvain",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]:
    """
    Cluster observations with the Louvain algorithm.

    Louvain clustering detects communities in an undirected neighborhood graph
    by optimizing modularity with a multilevel heuristic.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    resolution: float (default: 1.0)
        Resolution parameter controlling cluster granularity. Higher values
        usually lead to more clusters.
    neighbors_key: str, optional (default: 'neighbors')
        Key in `adata.uns` describing the precomputed neighborhood graph.
        Louvain reads the graph from
        `adata.obsp[adata.uns[neighbors_key]["connectivities_key"]]`.
    obsp: str, optional
        Key in `adata.obsp` storing the adjacency matrix. If provided, do not
        specify `neighbors_key`.
    weighted: bool (default: True)
        Whether edge weights from the adjacency graph are used.
    key_added: str (default: 'louvain')
        Column in `adata.obs` where cluster labels are stored.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by Louvain.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    **kwargs: Any
        Additional keyword arguments passed to `igraph.Graph.community_multilevel`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=50)
    >>> bt.sct.tl.neighbors(adata, representation="X_pca")
    >>> bt.sct.tl.louvain(adata, resolution=1.0)

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with Louvain clustering
        results added. Otherwise, updates `adata` in place and returns None.

        Louvain clustering results are stored in:

        - `adata.obs[key_added]`: categorical cluster labels;
        - `adata.uns[key_added]`: clustering metadata.

    References
    ----------
    Blondel et al. (2008). Fast unfolding of communities in large networks.
    Journal of Statistical Mechanics: Theory and Experiment, 2008(10), P10008.
    """

    igraph = import_module("igraph")

    resolution = _as_positive_number(resolution, "resolution")
    weighted = _as_boolean(weighted, "weighted")
    key_added = _as_string(key_added, "key_added")
    louvain_seed = _seed_to_integer(seed)

    adjacency_key, adjacency = _clustering_adjacency(
        adata=adata,
        neighbors_key=neighbors_key,
        obsp=obsp,
    )

    graph = _undirected_igraph_from_adjacency(igraph, adjacency)
    if weighted:
        weights = "weight"
    else:
        weights = None

    if louvain_seed is None:
        partition = graph.community_multilevel(
            weights=weights,
            resolution=resolution,
            **kwargs,
        )
    else:
        igraph.set_random_number_generator(random.Random(louvain_seed))
        try:
            partition = graph.community_multilevel(
                weights=weights,
                resolution=resolution,
                **kwargs,
            )
        finally:
            igraph.set_random_number_generator(None)

    adata = adata.copy() if copy else adata
    labels = pd.Categorical([str(cluster) for cluster in partition.membership])
    adata.obs[key_added] = pd.Series(labels, index=adata.obs_names, name=key_added)
    adata.uns[key_added] = {
        "params": {
            "method": "louvain",
            "resolution": resolution,
            "neighbors_key": neighbors_key,
            "obsp": adjacency_key if obsp is not None else None,
            "weighted": weighted,
            "seed": _format_random_state(seed),
            **kwargs,
        }
    }

    return adata if copy else None


@overload
def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    directed: bool = True,
    weighted: bool = True,
    n_iterations: AutoInteger = "auto",
    key_added: str = "leiden",
    seed: RandomStateSeed = 0,
    copy: Literal[False] = False,
    **kwargs: Any,
) -> None: ...


@overload
def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    directed: bool = True,
    weighted: bool = True,
    n_iterations: AutoInteger = "auto",
    key_added: str = "leiden",
    seed: RandomStateSeed = 0,
    copy: Literal[True],
    **kwargs: Any,
) -> AnnData: ...


@overload
def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    directed: bool = True,
    weighted: bool = True,
    n_iterations: AutoInteger = "auto",
    key_added: str = "leiden",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]: ...


@anndata_checker
def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    directed: bool = True,
    weighted: bool = True,
    n_iterations: AutoInteger = "auto",
    key_added: str = "leiden",
    seed: RandomStateSeed = 0,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[AnnData]:
    """
    Cluster observations with the Leiden algorithm.

    Leiden clustering detects communities in a neighborhood graph and is an
    improved version of the Louvain algorithm.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    resolution: float (default: 1.0)
        Resolution parameter controlling cluster granularity. Higher values
        usually lead to more clusters.
    neighbors_key: str, optional (default: 'neighbors')
        Key in `adata.uns` describing the precomputed neighborhood graph.
        Leiden reads the graph from
        `adata.obsp[adata.uns[neighbors_key]["connectivities_key"]]`.
    obsp: str, optional
        Key in `adata.obsp` storing the adjacency matrix. If provided, do not
        specify `neighbors_key`.
    directed: bool (default: True)
        Whether to treat the graph as directed.
    weighted: bool (default: True)
        Whether edge weights from the adjacency graph are used.
    n_iterations: int or "auto" (default: "auto")
        Number of Leiden iterations. If `"auto"`, optimize until convergence.
    key_added: str (default: 'leiden')
        Column in `adata.obs` where cluster labels are stored.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by Leiden.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    **kwargs: Any
        Additional keyword arguments passed to `leidenalg.find_partition`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=50)
    >>> bt.sct.tl.neighbors(adata, representation="X_pca")
    >>> bt.sct.tl.leiden(adata, resolution=1.0)

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with Leiden clustering
        results added. Otherwise, updates `adata` in place and returns None.

        Leiden clustering results are stored in:

        - `adata.obs[key_added]`: categorical cluster labels;
        - `adata.uns[key_added]`: clustering metadata.

    References
    ----------
    Traag et al. (2019). From Louvain to Leiden: guaranteeing well-connected
    communities. Scientific Reports, 9, 5233.
    """

    igraph = import_module("igraph")
    leidenalg = import_module("leidenalg")

    resolution = _as_positive_number(resolution, "resolution")
    weighted = _as_boolean(weighted, "weighted")
    directed = _as_boolean(directed, "directed")
    key_added = _as_string(key_added, "key_added")

    if isinstance(n_iterations, str):
        n_iterations = _as_literal(
            n_iterations,
            choices=("auto",),
            name="n_iterations",
        )
        leiden_iterations = -1
    else:
        if not isinstance(n_iterations, int) or isinstance(n_iterations, bool):
            raise TypeError(
                f"unsupported argument type for 'n_iterations': "
                f"expected {int} or 'auto' but received {type(n_iterations)}"
            )
        if n_iterations <= 0:
            raise ValueError(
                f"invalid argument value for 'n_iterations': "
                f"expected 'auto' or positive value but received {n_iterations!r}"
            )
        leiden_iterations = n_iterations

    leiden_seed = _seed_to_integer(seed)

    adjacency_key, adjacency = _clustering_adjacency(
        adata=adata,
        neighbors_key=neighbors_key,
        obsp=obsp,
    )

    graph = _igraph_from_adjacency(igraph, adjacency, directed=directed)

    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight" if weighted else None,
        n_iterations=leiden_iterations,
        seed=leiden_seed,
        resolution_parameter=resolution,
        **kwargs,
    )

    adata = adata.copy() if copy else adata
    labels = pd.Categorical([str(cluster) for cluster in partition.membership])
    adata.obs[key_added] = pd.Series(labels, index=adata.obs_names, name=key_added)
    adata.uns[key_added] = {
        "params": {
            "method": "leiden",
            "resolution": resolution,
            "neighbors_key": neighbors_key,
            "obsp": adjacency_key if obsp is not None else None,
            "directed": directed,
            "weighted": weighted,
            "n_iterations": n_iterations,
            "seed": _format_random_state(seed),
        }
    }

    return adata if copy else None


def _seed_to_integer(seed: RandomStateSeed) -> Optional[int]:

    _as_seed(seed)
    if seed is np.random:
        return int(np.random.randint(0, np.iinfo(np.int32).max))
    if isinstance(seed, numbers.Integral):
        return int(seed)
    if isinstance(seed, np.random.RandomState):
        return int(seed.randint(0, np.iinfo(np.int32).max))
    return None


def _clustering_adjacency(
    adata: AnnData,
    neighbors_key: Optional[str],
    obsp: Optional[str],
) -> Tuple[str, sparse.csr_matrix]:

    if neighbors_key is not None and obsp is not None:
        raise ValueError(
            "invalid argument combination: "
            "'neighbors_key' and 'obsp' cannot be both specified"
        )

    if obsp is not None:
        adjacency_key = _as_string(obsp, "obsp")
    else:
        neighbors_key = "neighbors" if neighbors_key is None else neighbors_key
        neighbors_key = _as_string(neighbors_key, "neighbors_key")
        if neighbors_key not in adata.uns:
            raise KeyError(
                f"key {neighbors_key!r} not found in adata.uns: "
                f"please run `bt.sct.tl.neighbors(..., key_added={neighbors_key!r})`"
            )
        neighbors = cast(Dict[str, Any], adata.uns[neighbors_key])
        adjacency_key = _as_string(
            neighbors["connectivities_key"],
            "connectivities_key",
        )

    if adjacency_key not in adata.obsp:
        raise KeyError(f"key {adjacency_key!r} not found in adata.obsp")

    pairwise_mtx = get_pairwise(adata, adjacency_key, axis="obs")
    adjacency = pairwise_mtx
    if sparse.issparse(adjacency):
        adjacency = cast(Any, adjacency).tocsr(copy=True)
        return adjacency_key, cast(sparse.csr_matrix, adjacency)
    return adjacency_key, sparse.csr_matrix(adjacency)


def _igraph_from_adjacency(
    igraph: Any,
    adjacency: sparse.csr_matrix,
    directed: bool,
) -> Any:

    n_vertices = cast(Tuple[int, int], adjacency.shape)[0]
    coo_adjacency = cast(Any, adjacency).tocoo(copy=True)
    graph = igraph.Graph(
        n=n_vertices,
        edges=list(
            zip(
                coo_adjacency.row.tolist(),
                coo_adjacency.col.tolist(),
            )
        ),
        directed=directed,
    )
    graph.es["weight"] = coo_adjacency.data.astype(float).tolist()

    return graph


def _undirected_igraph_from_adjacency(
    igraph: Any,
    adjacency: sparse.csr_matrix,
) -> Any:

    n_vertices = cast(Tuple[int, int], adjacency.shape)[0]
    upper = sparse.triu(adjacency, k=1, format="coo")
    lower = cast(Any, sparse.tril(adjacency, k=-1, format="coo")).T
    undirected = cast(Any, upper + lower).tocoo()
    graph = igraph.Graph(
        n=n_vertices,
        edges=list(
            zip(
                undirected.row.tolist(),
                undirected.col.tolist(),
            )
        ),
        directed=False,
    )
    graph.es["weight"] = undirected.data.astype(float).tolist()

    return graph
