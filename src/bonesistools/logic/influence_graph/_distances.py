#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Set, Tuple, Union

from ..._compat import Literal
from ._influence_graph import InfluenceGraph

SignedEdge = Tuple[Any, Any, int]
InfluenceGraphUniverse = Union[
    Literal["union", "complete", "complete_loopless"],
    InfluenceGraph,
]


def similarity(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    universe: InfluenceGraphUniverse = "union",
) -> float:
    """
    Compute similarity between two influence graphs.

    Influence graphs are compared through their signed edges
    `(source, target, sign)`. Similarity is the proportion of agreeing signed
    edges over the selected edge universe. Opposite signs between the same
    pair of nodes represent distinct edges; all other edge attributes are
    ignored.

    Examples
    --------
    >>> import bonesistools as bt
    >>> ig1 = bt.logic.ig.InfluenceGraph()
    >>> ig1.add_edge("A", "B", sign=1)
    >>> ig2 = ig1.copy()
    >>> ig2.add_edge("B", "A", sign=-1)
    >>> bt.logic.ig.similarity(ig1, ig2, universe="union")
    0.5
    >>> bt.logic.ig.similarity(ig1, ig2, universe="complete")
    0.875
    >>> prior = ig2.copy()
    >>> prior.add_edge("A", "A", sign=1)
    >>> bt.logic.ig.similarity(ig1, ig2, universe=prior)
    0.6666666666666666

    Parameters
    ----------
    ig1, ig2: InfluenceGraph
        Influence graphs to compare.
    universe: {"union", "complete", "complete_loopless"} or InfluenceGraph
        Signed-edge universe used for the comparison. Default is `union`.

        - `union` uses signed edges observed in at least one graph. This gives
          classical signed-edge Jaccard similarity.
        - `complete` uses every positive and negative directed edge between
          nodes appearing in either graph, including self-loops. This gives
          normalized Hamming similarity.
        - `complete_loopless` uses the same complete signed-edge universe but
          excludes self-loops. This gives loopless normalized Hamming
          similarity.
        - An `InfluenceGraph` uses its signed edges as an explicit universe. It
          must contain every signed edge from both input graphs.

    Returns
    -------
    float
        Signed-edge similarity in the interval `[0, 1]`.

    Notes
    -----
    Choosing `universe="union"` yields classical Jaccard similarity. Choosing
    a complete or explicit universe yields normalized Hamming similarity over
    that signed-edge universe, where both shared presences and shared absences
    count as agreements.

    Using the same explicit universe is recommended for pairwise similarity
    matrices when graph pairs may contain different node sets.

    Raises
    ------
    ValueError
        If `complete_loopless` is selected for an input containing a self-loop,
        an explicit universe omits an input edge, or `universe` is unsupported.
    """

    edges1 = _signed_edge_set(ig1)
    edges2 = _signed_edge_set(ig2)
    difference_size = len(edges1 ^ edges2)
    universe_size = _edge_universe_size(
        ig1,
        ig2,
        edges1,
        edges2,
        universe,
    )

    if universe_size == 0:
        return 1.0

    agreement_size = universe_size - difference_size

    return agreement_size / universe_size


def distance(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    universe: InfluenceGraphUniverse = "union",
) -> float:
    """
    Compute distance between two influence graphs.

    Influence graphs are compared through their signed edges
    `(source, target, sign)`. Distance is the proportion of disagreeing signed
    edges over the selected edge universe. Opposite signs between the same
    pair of nodes represent distinct edges; all other edge attributes are
    ignored.

    Examples
    --------
    >>> import bonesistools as bt
    >>> ig1 = bt.logic.ig.InfluenceGraph()
    >>> ig1.add_edge("A", "B", sign=1)
    >>> ig2 = ig1.copy()
    >>> ig2.add_edge("B", "A", sign=-1)
    >>> bt.logic.ig.distance(ig1, ig2, universe="union")
    0.5
    >>> bt.logic.ig.distance(ig1, ig2, universe="complete")
    0.125
    >>> prior = ig2.copy()
    >>> prior.add_edge("A", "A", sign=1)
    >>> bt.logic.ig.distance(ig1, ig2, universe=prior)
    0.3333333333333333

    Parameters
    ----------
    ig1, ig2: InfluenceGraph
        Influence graphs to compare.
    universe: {"union", "complete", "complete_loopless"} or InfluenceGraph
        Signed-edge universe used for the comparison. Default is `union`.

        - `union` uses signed edges observed in at least one graph. This gives
          classical signed-edge Jaccard distance.
        - `complete` uses every positive and negative directed edge between
          nodes appearing in either graph, including self-loops. This gives
          normalized Hamming distance.
        - `complete_loopless` uses the same complete signed-edge universe but
          excludes self-loops. This gives loopless normalized Hamming distance.
        - An `InfluenceGraph` uses its signed edges as an explicit universe. It
          must contain every signed edge from both input graphs.

    Returns
    -------
    float
        Signed-edge distance in the interval `[0, 1]`.

    Notes
    -----
    Choosing `universe="union"` yields classical Jaccard distance. Choosing a
    complete or explicit universe yields normalized Hamming distance over that
    signed-edge universe. Shared absences enlarge the universe without
    contributing a disagreement, so Hamming distance can be smaller for sparse
    graphs.

    Using the same explicit universe is recommended for pairwise distance
    matrices when graph pairs may contain different node sets.

    Raises
    ------
    ValueError
        If `complete_loopless` is selected for an input containing a self-loop,
        an explicit universe omits an input edge, or `universe` is unsupported.
    """

    edges1 = _signed_edge_set(ig1)
    edges2 = _signed_edge_set(ig2)
    difference_size = len(edges1 ^ edges2)
    universe_size = _edge_universe_size(
        ig1,
        ig2,
        edges1,
        edges2,
        universe,
    )

    if universe_size == 0:
        return 0.0

    return difference_size / universe_size


def _signed_edge_set(ig: InfluenceGraph) -> Set[SignedEdge]:

    return {
        (source, target, InfluenceGraph._normalize_sign(sign))
        for source, target, sign in ig.edges(data="sign")
    }


def _edge_universe_size(
    ig1: InfluenceGraph,
    ig2: InfluenceGraph,
    edges1: Set[SignedEdge],
    edges2: Set[SignedEdge],
    universe: InfluenceGraphUniverse,
) -> int:

    input_edges = edges1 | edges2

    if universe == "union":
        return len(input_edges)

    nodes = set(ig1.nodes()) | set(ig2.nodes())

    if universe == "complete":
        return 2 * len(nodes) ** 2

    if universe == "complete_loopless":
        if any(source == target for source, target, _ in input_edges):
            raise ValueError(
                "`universe='complete_loopless'` cannot contain input self-loops"
            )

        return 2 * len(nodes) * (len(nodes) - 1)

    if isinstance(universe, InfluenceGraph):
        universe_edges = _signed_edge_set(universe)

        if not input_edges <= universe_edges:
            raise ValueError(
                "`universe` must contain every signed edge from both input graphs"
            )

        return len(universe_edges)

    raise ValueError(f"unsupported edge universe: {universe!r}")
