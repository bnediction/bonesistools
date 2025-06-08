#!/usr/bin/env python

from collections import namedtuple

from typing import Any, Sequence, Optional
from numbers import Number
from ._typing import Graph

import networkx as nx
from ._algorithms import depth_first_extraction

def get_edge_sign(
    grn: Graph,
    source: Any,
    target: Any
):
    """
    Get sign between source and target nodes.

    Parameters
    ----------
    grn: nx.Graph
        NetworkX Graph object.
    source: Any
        Starting node for edge.
    target: Any
        Ending node for edge.

    Returns
    -------
    Return sign of a source node upon a target node.
    """

    edge_data = grn.get_edge_data(source, target)
    signs = {value["sign"] for value in edge_data.values()}
    for sign in signs:
        if sign not in [-1, 1]:
            raise ValueError(f"invalid value for edge attribute 'sign': expected 0 or 1 but received {sign}")
    if len(signs) == 1:
        return list(signs)[0]
    else:
        return 0

def get_path_sign(grn: Graph, *nodes) -> int:
    """
    Compute whether the source node has a positive or negative effect upon the target node within a path.

    Parameters
    ----------
    grn: nx.Graph
        NetworkX Graph object.
    *nodes
        nodes in a specific order depicting an existing path.

    Returns
    -------
    Return 1 if effect is positive and -1 if effect is negative.
    """

    path_sign = 1
    u = nodes[0]
    for v in nodes[1:]:
        _sign = get_edge_sign(grn, u, v)
        if _sign == -1:
            path_sign = _sign * path_sign
        elif _sign == 0:
            return 0
        elif _sign == 1:
            pass
        else:
            raise ValueError(f"invalid value for edge sign between nodes {u} and {v}: expected -1, 0 or 1 but got {_sign}")
        u = v
    return path_sign

def path_to_string(grn: nx.Graph, *nodes) -> str:
    """
    Get a human-readable string describing a path.

    Parameters
    ----------
    grn: nx.Graph
        NetworkX Graph object.
    *nodes
        nodes in a specific order depicting an existing path.

    Returns
    -------
    Return a string with nodes separated by an arrow depicting the positive, negative or bi-signed effect.
    """

    u = nodes[0]
    string = str(u)
    for v in nodes[1:]:
        sign = get_edge_sign(grn, u, v)
        if sign == -1:
            string += f" -| {v}"
        elif sign == 1:
            string += f" -> {v}"
        else:
            string += f" -- {v}"
        u = v
    return string

def statistics(
    grn: Graph,
    weights: Sequence[Number],
    radius: int = 3,
    genes: Optional[Sequence[str]] = None,
    allow_null_path_sign: Optional[bool] = None
) -> dict:
    """
    Compute statistics (score, path number and maximum score) upon all existing paths within a radius for all gene pairwise in a gene set.

    Parameters
    ----------
    grn: nx.Graph
        NetworkX Graph object.
    weights: Sequence[Number]
        Specify path length weighting in the computation of scores.
    radius: int (default: 3)
        Specify the maximum search depth when sampling paths.
    genes: Sequence[str] (optional, default: None)
        Set of genes being considered.

    Returns
    -------
    Return a two-keywords dictionary as follows:
    - 1st key: source
    - 2nd key: target
    - value: statistics
    """

    interactions = dict()
    genes = genes if genes is not None else set(grn.nodes)
    for source in genes:
        targets = genes.difference([source])
        paths_from_source = depth_first_extraction(graph=grn, source=source, limit_depth=radius)
        _interactions_from_source = {target: [0, 0, 0] for target in targets}
        for _path in paths_from_source:
            _target = _path[-1]
            if _target in targets:
                _path_sign = get_path_sign(grn, *_path)
                if allow_null_path_sign is True or _path_sign != 0:
                    _score, _maxscore, _path_number = _interactions_from_source[_target]
                    _weight = weights[len(_path)-2]
                    _score += _path_sign * _weight
                    _maxscore += _weight
                    _path_number += 1
                    _interactions_from_source[_target] = [_score, _maxscore, _path_number]
        for gene, value in _interactions_from_source.items():
            _score, _maxscore, _path_number = value
            _interactions_from_source[gene] = namedtuple("Tuple", ["score", "maxscore", "path_number"])(_score, _maxscore, _path_number)
        interactions[source] = _interactions_from_source
        del _interactions_from_source
    return interactions
