#!/usr/bin/env python

from __future__ import annotations

import networkx as nx


def walks_from(
    graph: nx.Graph,
    source: str,
    max_depth: int = 3,
) -> list[list[str]]:
    """
    Return directed walks starting from a source node up to a maximum depth.

    Walks are enumerated through a recursive depth-first traversal strategy with
    bounded depth. Unlike NetworkX functions such as `nx.all_simple_paths`,
    nodes may be revisited multiple times within the same walk. This allows
    feedback traversals and cyclic propagation patterns to be represented.

    Consequently, the returned walks are not restricted to simple paths and may
    contain repeated nodes whenever the graph contains directed cycles.

    Examples
    --------
    Consider the following influence graph:

        A → B → C
            ↑   ↓
            └───┘

    >>> walks_from(G, "A", max_depth=4)
    [
        ['A', 'B'],
        ['A', 'B', 'C'],
        ['A', 'B', 'C', 'B'],
        ['A', 'B', 'C', 'B', 'C'],
    ]

    In contrast:
    >>> list(nx.all_simple_paths(G, "A", "C"))
    [['A', 'B', 'C']]
    because simple paths forbid repeated nodes.

    Parameters
    ----------
    graph: nx.DiGraph
        Directed graph.
    source: str
        Source node.
    max_depth: int (default: 3)
        Maximum number of traversed edges.

    Returns
    -------
    list of list of str
        Directed walks starting from `source`.
    """

    walks = []

    def explore(
        node: str,
        path: list[str],
        remaining_depth: int,
    ) -> None:

        if remaining_depth == 0:
            return

        for successor in getattr(graph, "successors")(node):
            new_path = path + [successor]
            walks.append(new_path)
            explore(successor, new_path, remaining_depth - 1)

    explore(source, [source], max_depth)

    return walks
