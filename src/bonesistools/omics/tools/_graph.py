#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, cast, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.sparse.csgraph import minimum_spanning_tree

from ..._compat import Literal
from ..._validation import _as_string
from .._typing import anndata_checker
from ._utils import get_pairwise

_DENSE_PAGA_GROUP_PAIRS_LIMIT = 1_000_000


@overload
def paga(
    adata: AnnData,
    groupby: str,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    key_added: str = "paga",
    copy: Literal[False] = False,
) -> None: ...


@overload
def paga(
    adata: AnnData,
    groupby: str,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    key_added: str = "paga",
    copy: Literal[True],
) -> AnnData: ...


@overload
def paga(
    adata: AnnData,
    groupby: str,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    key_added: str = "paga",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def paga(
    adata: AnnData,
    groupby: str,
    *,
    neighbors_key: Optional[str] = "neighbors",
    obsp: Optional[str] = None,
    key_added: str = "paga",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Compute partition-based graph abstraction connectivities.

    PAGA quantifies the coarse-grained connectivity between groups of
    observations by comparing observed inter-group edges in the neighborhood
    graph to the number expected under a random null model.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation column in `adata.obs` defining groups.
    neighbors_key: str, optional (default: 'neighbors')
        Key in `adata.uns` describing the precomputed neighborhood graph.
        PAGA reads distances from
        `adata.obsp[adata.uns[neighbors_key]["distances_key"]]`.
    obsp: str, optional
        Key in `adata.obsp` storing the adjacency matrix. If provided, do not
        specify `neighbors_key`.
    key_added: str (default: 'paga')
        Key in `adata.uns` where PAGA results are stored.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with PAGA results added.
        Otherwise, updates `adata` in place and returns None.

        PAGA results are stored in:

        - `adata.uns[key_added]["connectivities"]`: abstracted graph
          connectivities;
        - `adata.uns[key_added]["connectivities_tree"]`: maximum-spanning
          tree-like subgraph;
        - `adata.uns[key_added]["params"]`: PAGA metadata.

    Examples
    --------
    >>> bt.omics.tl.pca(adata, n_components=50)
    >>> bt.omics.tl.neighbors(adata, representation="X_pca")
    >>> bt.omics.tl.leiden(adata, resolution=1.0)
    >>> bt.omics.tl.paga(adata, groupby="leiden")

    References
    ----------
    Wolf et al. (2019). PAGA: graph abstraction reconciles clustering with
    trajectory inference through a topology preserving map of single cells.
    Genome Biology, 20(1), 59.
    """

    groupby = _as_string(groupby, "groupby")
    key_added = _as_string(key_added, "key_added")
    adjacency_key, adjacency = _paga_adjacency(
        adata=adata,
        neighbors_key=neighbors_key,
        obsp=obsp,
    )

    labels = pd.Categorical(adata.obs[groupby]).remove_unused_categories()
    if np.any(labels.codes < 0):
        raise ValueError(
            f"invalid observation labels in adata.obs[{groupby!r}]: "
            "PAGA requires non-missing group labels"
        )

    connectivities, connectivities_tree = _paga_connectivities(
        adjacency=adjacency,
        codes=np.asarray(labels.codes, dtype=np.int64),
        n_groups=len(labels.categories),
    )

    adata = adata.copy() if copy else adata
    adata.uns[key_added] = {
        "connectivities": connectivities,
        "connectivities_tree": connectivities_tree,
        "params": {
            "method": "paga",
            "groupby": groupby,
            "neighbors_key": neighbors_key,
            "obsp": adjacency_key if obsp is not None else None,
        },
    }

    return adata if copy else None


def _paga_adjacency(
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
                f"please run `bt.omics.tl.neighbors(..., key_added={neighbors_key!r})`"
            )
        neighbors = cast(Dict[str, Any], adata.uns[neighbors_key])
        adjacency_key = _as_string(neighbors["distances_key"], "distances_key")

    if adjacency_key not in adata.obsp:
        raise KeyError(f"key {adjacency_key!r} not found in adata.obsp")

    pairwise_mtx = get_pairwise(adata, adjacency_key, axis="obs")
    if sparse.issparse(pairwise_mtx):
        adjacency = cast(Any, pairwise_mtx).tocsr(copy=False)
    else:
        adjacency = sparse.csr_matrix(pairwise_mtx)

    expected_shape = (adata.n_obs, adata.n_obs)
    if adjacency.shape != expected_shape:
        raise ValueError(
            "invalid adjacency matrix shape: "
            f"expected {expected_shape} but received {adjacency.shape}"
        )

    if adjacency.nnz > 0 and np.any(adjacency.data == 0):
        adjacency = adjacency.copy()
        adjacency.eliminate_zeros()

    return adjacency_key, cast(sparse.csr_matrix, adjacency)


def _paga_connectivities(
    adjacency: sparse.csr_matrix,
    codes: np.ndarray,
    n_groups: int,
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix]:

    if codes.shape[0] < 2:
        raise ValueError(
            "invalid number of observations: PAGA requires at least two observations"
        )

    graph = adjacency.tocoo(copy=False)
    source_codes = codes[graph.row]
    target_codes = codes[graph.col]
    inter_cluster = source_codes != target_codes

    inner_edges = np.bincount(
        source_codes[~inter_cluster],
        minlength=n_groups,
    ).astype(float)
    group_sizes = np.bincount(codes, minlength=n_groups).astype(float)
    group_pair_count = n_groups * n_groups

    if group_pair_count <= _DENSE_PAGA_GROUP_PAIRS_LIMIT:
        inter_source_codes = source_codes[inter_cluster]
        inter_target_codes = target_codes[inter_cluster]
        inter_edges = np.bincount(
            inter_source_codes * n_groups + inter_target_codes,
            minlength=group_pair_count,
        ).reshape(n_groups, n_groups)
        inter_edges = inter_edges.astype(float, copy=False)

        outgoing_edges = inter_edges.sum(axis=1)
        observed_edges = inter_edges + inter_edges.T
        rows, columns = np.nonzero(observed_edges)
        observed = observed_edges[rows, columns]
    else:
        inter_edges = sparse.csr_matrix(
            (
                np.ones(int(inter_cluster.sum()), dtype=float),
                (source_codes[inter_cluster], target_codes[inter_cluster]),
            ),
            shape=(n_groups, n_groups),
        )
        outgoing_edges = np.asarray(inter_edges.sum(axis=1)).ravel()
        observed_edges = cast(sparse.csr_matrix, inter_edges + inter_edges.T).tocoo()
        rows = observed_edges.row
        columns = observed_edges.col
        observed = observed_edges.data.astype(float, copy=False)

    total_obs = int(group_sizes.sum())
    total_edges = inner_edges + outgoing_edges
    expected = (
        total_edges[rows] * group_sizes[columns]
        + total_edges[columns] * group_sizes[rows]
    ) / (total_obs - 1)
    values = np.ones_like(observed, dtype=float)
    nonzero_expected = expected != 0
    values[nonzero_expected] = observed[nonzero_expected] / expected[nonzero_expected]
    values = np.minimum(values, 1.0)

    connectivities = sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(n_groups, n_groups),
    )
    connectivities.eliminate_zeros()

    return connectivities, _paga_connectivities_tree(connectivities)


def _paga_connectivities_tree(connectivities: sparse.csr_matrix) -> sparse.csr_matrix:

    if connectivities.nnz == 0:
        return sparse.csr_matrix(connectivities.shape, dtype=float)

    inverse = connectivities.copy()
    inverse.data = 1.0 / inverse.data
    tree = minimum_spanning_tree(inverse).tocoo()
    tree_values = np.asarray(connectivities[tree.row, tree.col]).ravel()
    return sparse.csr_matrix(
        (tree_values, (tree.row, tree.col)),
        shape=connectivities.shape,
    )
