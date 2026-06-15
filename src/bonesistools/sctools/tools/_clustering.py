#!/usr/bin/env python

from __future__ import annotations

import numbers
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
    _as_positive_number,
    _as_seed,
    _as_string,
)
from .._metadata import _format_random_state
from .._typing import anndata_checker
from ._utils import get_pairwise


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

    _as_seed(seed)
    leiden_seed = None
    if seed is np.random:
        leiden_seed = int(np.random.randint(0, np.iinfo(np.int32).max))
    elif isinstance(seed, numbers.Integral):
        leiden_seed = int(seed)
    elif isinstance(seed, np.random.RandomState):
        leiden_seed = int(seed.randint(0, np.iinfo(np.int32).max))

    adjacency_key, adjacency = _leiden_adjacency(
        adata=adata,
        neighbors_key=neighbors_key,
        obsp=obsp,
    )

    adjacency = adjacency.tocoo(copy=True)
    graph = igraph.Graph(
        n=adata.n_obs,
        edges=list(zip(adjacency.row.tolist(), adjacency.col.tolist())),
        directed=directed,
    )
    weights = adjacency.data.astype(float).tolist()
    if weighted:
        graph.es["weight"] = weights

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


def _leiden_adjacency(
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
