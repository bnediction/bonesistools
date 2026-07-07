#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Optional, Set, Tuple

from ._influence_graph import InfluenceGraph

SignedEdge = Tuple[Any, Any, int]


def similarity(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    metric: str = "jaccard",
    domain: Optional[InfluenceGraph] = None,
) -> float:
    """
    Compute similarity between two influence graphs.

    The comparison is performed on signed edges represented as
    `(source, target, sign)`. Two edges with the same source and target but
    different signs are therefore treated as different edges.

    Parameters
    ----------
    ig1, ig2: InfluenceGraph
        Influence graphs to compare.
    metric: {"jaccard", "hamming"} (default: 'jaccard')
        Similarity metric.
    domain: InfluenceGraph, optional
        Signed-edge universe used when `metric="hamming"`. If `None`, use all
        signed directed edges between nodes in `union(V1, V2)`, including
        signed self-edges. If provided, the domain, corresponding to a reference
        graph, must be a strict signed-edge superset of both input graphs.

    Examples
    --------
    >>> bt.bpy.ig.similarity(ig1, ig2)

    >>> dorothea = bt.dbs.omnipath.dorothea()
    >>> bt.bpy.ig.similarity(ig1, ig2, metric="hamming", domain=dorothea)

    Returns
    -------
    float
        Signed-edge similarity.
    """

    edges1 = _signed_edge_set(ig1)
    edges2 = _signed_edge_set(ig2)

    if metric == "jaccard":
        if domain is not None:
            raise ValueError("`domain` is only supported with metric='hamming'")

        intersection = edges1 & edges2
        union = edges1 | edges2

        if not union:
            return 1.0

        return len(intersection) / len(union)

    if metric == "hamming":
        domain_size = _hamming_domain_size(ig1, ig2, edges1, edges2, domain)

        if domain_size == 0:
            return 1.0

        symmetric_difference = edges1 ^ edges2

        return 1.0 - len(symmetric_difference) / domain_size

    else:
        raise ValueError(f"unsupported similarity metric: {metric!r}")


def distance(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    metric: str = "jaccard",
    domain: Optional[InfluenceGraph] = None,
) -> float:
    """
    Compute distance between two influence graphs.

    The comparison is performed on signed edges represented as
    `(source, target, sign)`. Two edges with the same source and target but
    different signs are therefore treated as different edges.

    Parameters
    ----------
    ig1, ig2: InfluenceGraph
        Influence graphs to compare.
    metric: {"jaccard", "hamming"} (default: 'jaccard')
        Distance metric.
    domain: InfluenceGraph, optional
        Signed-edge universe used when `metric="hamming"`. If `None`, use all
        signed directed edges between nodes in `union(V1, V2)`, including
        signed self-edges. If provided, the domain, corresponding to a reference
        graph, must be a strict signed-edge superset of both input graphs.

    Examples
    --------
    >>> bt.bpy.ig.distance(ig1, ig2)

    >>> dorothea = bt.dbs.omnipath.dorothea()
    >>> bt.bpy.ig.distance(ig1, ig2, metric="hamming", domain=dorothea)

    Returns
    -------
    float
        Signed-edge distance.
    """

    if metric == "jaccard":
        return 1.0 - similarity(ig1, ig2, metric=metric, domain=domain)

    if metric == "hamming":
        edges1 = _signed_edge_set(ig1)
        edges2 = _signed_edge_set(ig2)

        domain_size = _hamming_domain_size(ig1, ig2, edges1, edges2, domain)

        if domain_size == 0:
            return 0.0

        symmetric_difference = edges1 ^ edges2

        return len(symmetric_difference) / domain_size

    else:
        raise ValueError(f"unsupported distance metric: {metric!r}")


def _signed_edge_set(ig: InfluenceGraph) -> Set[SignedEdge]:

    return {
        _signed_edge_identity(source, target, sign)
        for source, target, sign in ig.edges(data="sign")
    }


def _hamming_domain_size(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    edges1: Set[SignedEdge],
    edges2: Set[SignedEdge],
    domain: Optional[InfluenceGraph],
) -> int:

    if domain is None:
        return _implicit_hamming_domain_size(ig1, ig2)

    return len(_validate_hamming_domain(edges1, edges2, domain))


def _signed_edge_identity(source: Any, target: Any, sign: Any) -> SignedEdge:

    return (source, target, InfluenceGraph._normalize_sign(sign))


def _implicit_hamming_domain_size(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
) -> int:

    nodes = set(ig1.nodes()) | set(ig2.nodes())

    return 2 * len(nodes) ** 2


def _validate_hamming_domain(
    edges1: Set[SignedEdge],
    edges2: Set[SignedEdge],
    domain: InfluenceGraph,
) -> Set[SignedEdge]:

    domain_edges = _signed_edge_set(domain)

    if not edges1 < domain_edges or not edges2 < domain_edges:
        raise ValueError(
            "`domain` must be a strict signed-edge superset of both input graphs"
        )

    return domain_edges
