#!/usr/bin/env python

from typing import Any, Sequence, Optional
from numbers import Number
from collections import namedtuple

import networkx as nx
from ._algorithms import depth_first_extraction

def get_edge_sign(graph: nx.Graph, source: Any, target: Any):
    """
    Get sign between source and target genes.

    Parameters
    ----------
    graph
        NetworkX graph
    source
        node
    target
        node

    Returns
    -------
    Return sign of a source gene upon a target gene. 
    """

    edge_data = graph.get_edge_data(source, target)
    signs = {value["sign"] for value in edge_data.values()}
    for sign in signs:
        if sign not in [-1, 1]:
            raise ValueError("edge attribute `sign` is not equal to -1 or 1")
    if len(signs) == 1:
        return list(signs)[0]
    else:
        return 0

def get_path_sign(graph: nx.Graph, *nodes) -> int:
    """
    Compute whether the source gene has a positive or negative effect upon the target gene within a path.

    Parameters
    ----------
    graph
        NetworkX graph
    *nodes
        nodes in a specific order depicting an existing path

    Returns
    -------
    Return 1 if effect is positive and -1 if effect is negative.
    """

    path_sign = 1
    u = nodes[0]
    for v in nodes[1:]:
        _sign = get_edge_sign(graph, u, v)
        if _sign == -1:
            path_sign = _sign * path_sign
        elif _sign == 0:
            return 0
        elif _sign == 1:
            pass
        else:
            raise ValueError("value of `sign` between {u} and {v} genes is not equal to -1, 0 or 1: `sign` = {_sign}")
        u = v
    return path_sign

def path_to_string(graph: nx.Graph, *nodes) -> str:
    """
    Get a human-readable string describing a path.

    Parameters
    ----------
    graph
        NetworkX graph
    *nodes
        nodes in a specific order depicting an existing path

    Returns
    -------
    Return a string with nodes separated by an arrow depicting the positive, negative or bi-signed effect.
    """

    u = nodes[0]
    string = str(u)
    for v in nodes[1:]:
        sign = get_edge_sign(graph, u, v)
        if sign == -1:
            string += f" -| {v}"
        elif sign == 1:
            string += f" -> {v}"
        else:
            string += f" -- {v}"
        u = v
    return string

def scoring(
    graph: nx.Graph,
    weights: Sequence[Number],
    radius: int = 3,
    gene_set: Optional[Sequence[str]] = None,
    allow_null_path_sign: Optional[bool] = None
) -> dict:
    """
    Compute statistics (score, path number and maximum score) upon all existing paths within a radius for all gene pairwise in a gene set.

    Parameters
    ----------
    graph
        NetworkX graph
    weights
        specify path length weighting in the computation of scores.
    radius
        specify the maximum search depth when sampling paths.
    gene_set
        set of genes

    Returns
    -------
    Return a two-keyword dictionary where one considers:
    - first key: source
    - second key: target
    - value: namedtuple with score, path_number and maxscore values
    """

    interactions = dict()
    gene_set = gene_set if gene_set is not None else set(graph.nodes)
    for source in gene_set:
        targets = gene_set.difference([source])
        paths_from_source = depth_first_extraction(graph=graph, source=source, limit_depth=radius)
        _interactions_from_source = {target: [0, 0, 0] for target in targets}
        for _path in paths_from_source:
            _target = _path[-1]
            if _target in targets:
                _path_sign = get_path_sign(graph, *_path)
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
