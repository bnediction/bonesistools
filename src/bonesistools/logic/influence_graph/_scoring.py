#!/usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

import networkx as nx

from ..._validation import _as_positive_integer, _as_probability


@dataclass(frozen=True)
class _InteractionScore:
    """
    Signed influence score accumulated from directed walks.

    Parameters
    ----------
    score: float
        Weighted signed sum of contributing walks.
    total_weight: float
        Sum of absolute walk weights.
    path_number: int
        Number of signed walks contributing to the interaction.
    """

    score: float
    total_weight: float
    path_number: int

    @property
    def normalized_score(self) -> float:
        """
        Return the score normalized by total walk weight.

        Returns
        -------
        float
            Value in [-1, 1], or 0.0 if no walk contributes.
        """

        if self.total_weight == 0:
            return 0.0

        return self.score / self.total_weight


_EMPTY_INTERACTION_SCORE = _InteractionScore(
    score=0.0,
    total_weight=0.0,
    path_number=0,
)


def interaction_scores_from_walks(
    graph: nx.Graph[Any],
    genes: Optional[Iterable[str]] = None,
    max_depth: int = 3,
    weights: Optional[Sequence[float]] = None,
) -> Dict[str, Dict[str, _InteractionScore]]:
    """
    Estimate signed interaction scores from bounded directed walks.

    Walks are aggregated from each source gene up to `max_depth`. A walk
    contributes the product of its edge signs multiplied by a length-dependent
    weight. Parallel signed edges in a MultiDiGraph are counted independently.

    Parameters
    ----------
    graph: nx.Graph
        Directed signed graph.
    genes: Iterable[str], optional
        Genes used as sources and targets. If `None`, use all graph nodes.
    max_depth: int (default: 3)
        Maximum number of traversed edges.
    weights: Sequence[float], optional
        Walk weights indexed by path length minus one. If `None`, use
        `1 / path_length**2`.

    Returns
    -------
    dict
        Nested mapping `source -> target -> score object`. Each score object
        exposes `score`, `total_weight`, `path_number` and `normalized_score`.

    Raises
    ------
    ValueError
        If `max_depth` is not positive, `weights` is too short or an edge sign
        is invalid.
    """

    max_depth = _as_positive_integer(max_depth, "max_depth")

    weights = _default_weights(max_depth) if weights is None else weights
    all_nodes_are_sources = genes is None
    genes = list(graph.nodes if genes is None else genes)
    gene_set = set(genes)

    scores: Dict[str, Dict[str, _InteractionScore]] = {gene: {} for gene in genes}
    accumulators: Dict[str, Dict[str, Tuple[float, float, int]]] = {
        gene: {} for gene in genes
    }
    adjacency: Dict[str, Tuple[Tuple[str, int, int], ...]] = (
        {node: _signed_successor_counts(graph, node) for node in graph}
        if all_nodes_are_sources
        else {}
    )

    for source in genes:
        if source not in graph:
            continue

        frontier: Dict[str, Tuple[int, int]] = {source: (1, 0)}

        for path_length in range(1, max_depth + 1):
            next_frontier: Dict[str, Tuple[int, int]] = {}

            for node, signed_counts in frontier.items():
                positive_count, negative_count = signed_counts

                if positive_count == 0 and negative_count == 0:
                    continue

                if node not in adjacency:
                    adjacency[node] = _signed_successor_counts(graph, node)

                for successor, positive_edges, negative_edges in adjacency[node]:
                    previous_positive, previous_negative = next_frontier.get(
                        successor,
                        (0, 0),
                    )
                    next_frontier[successor] = (
                        previous_positive
                        + positive_count * positive_edges
                        + negative_count * negative_edges,
                        previous_negative
                        + negative_count * positive_edges
                        + positive_count * negative_edges,
                    )

            if not next_frontier:
                break

            weight = _path_weight(path_length, weights)

            for target, signed_counts in next_frontier.items():
                if target == source:
                    continue

                if target not in gene_set:
                    continue

                positive_count, negative_count = signed_counts
                path_number = positive_count + negative_count

                if path_number == 0:
                    continue

                score, total_weight, total_paths = accumulators[source].get(
                    target,
                    (0.0, 0.0, 0),
                )
                accumulators[source][target] = (
                    score + (positive_count - negative_count) * weight,
                    total_weight + path_number * weight,
                    total_paths + path_number,
                )

            frontier = next_frontier

    for source, targets in accumulators.items():
        for target, values in targets.items():
            scores[source][target] = _InteractionScore(
                score=values[0],
                total_weight=values[1],
                path_number=values[2],
            )

    return scores


def infer_signed_interactions(
    scores: Mapping[str, Mapping[str, _InteractionScore]],
    genes: Optional[Iterable[str]] = None,
    minimum_path_number: int = 1,
    threshold: float = 0.75,
    allow_bidirectional: bool = False,
) -> List[Tuple[str, str, Dict[str, int]]]:
    """
    Infer signed directed interactions from normalized interaction scores.

    Parameters
    ----------
    scores: Mapping
        Nested interaction scores as returned by
        `interaction_scores_from_walks()`. Score objects must expose
        `path_number` and `normalized_score`.
    genes: Iterable[str], optional
        Genes to compare. If `None`, use keys from `scores`.
    minimum_path_number: int (default: 1)
        Minimum number of contributing walks required for an interaction.
    threshold: float (default: 0.75)
        Minimum absolute normalized score.
    allow_bidirectional: bool (default: False)
        If `True`, allow both directions of a gene pair to be retained.

    Returns
    -------
    list of tuple
        Interactions as `(source, target, {"sign": sign})`.

    Raises
    ------
    ValueError
        If `minimum_path_number` is less than 1 or `threshold` is outside
        [0, 1].
    """

    minimum_path_number = _as_positive_integer(
        minimum_path_number,
        "minimum_path_number",
    )

    threshold = _as_probability(threshold, "threshold")

    genes = list(scores if genes is None else genes)
    interactions = []
    gene_positions = {gene: index for index, gene in enumerate(genes)}
    if len(gene_positions) == len(genes):
        scored_pairs = set()
        for source in genes:
            source_index = gene_positions[source]
            for target in scores.get(source, {}):
                if target == source or target not in gene_positions:
                    continue
                target_index = gene_positions[target]
                scored_pairs.add(
                    (min(source_index, target_index), max(source_index, target_index))
                )
        ordered_pairs = sorted(scored_pairs)
    else:
        ordered_pairs = [
            (source_index, target_index)
            for source_index, source in enumerate(genes)
            for target_index, target in enumerate(
                genes[source_index + 1 :],
                start=source_index + 1,
            )
            if target in scores.get(source, {}) or source in scores.get(target, {})
        ]

    for source_index, target_index in ordered_pairs:
        source = genes[source_index]
        target = genes[target_index]
        forward = scores.get(source, {}).get(target, _EMPTY_INTERACTION_SCORE)
        backward = scores.get(target, {}).get(source, _EMPTY_INTERACTION_SCORE)

        candidates = [
            (source, target, forward, backward),
            (target, source, backward, forward),
        ]

        for candidate_source, candidate_target, score, opposite in candidates:
            normalized_score = score.normalized_score
            if (
                score.path_number < minimum_path_number
                or normalized_score == 0
                or abs(normalized_score) < threshold
            ):
                continue

            if not allow_bidirectional and abs(normalized_score) <= abs(
                opposite.normalized_score
            ):
                continue

            sign = 1 if normalized_score > 0 else -1
            interactions.append((candidate_source, candidate_target, {"sign": sign}))

    return interactions


def infer_signed_interactions_from_walks(
    graph: nx.Graph[Any],
    genes: Optional[Iterable[str]] = None,
    max_depth: int = 3,
    weights: Optional[Sequence[float]] = None,
    minimum_path_number: int = 1,
    threshold: float = 0.75,
    allow_bidirectional: bool = False,
) -> List[Tuple[str, str, Dict[str, int]]]:
    """
    Infer signed interactions directly from bounded directed walks.

    This convenience function combines `interaction_scores_from_walks()` and
    `infer_signed_interactions()`. Use `interaction_scores_from_walks()`
    directly when intermediate score attributes should be inspected.

    Parameters
    ----------
    graph: nx.Graph
        Directed signed graph.
    genes: Iterable[str], optional
        Genes used as sources and targets. If `None`, use all graph nodes.
    max_depth: int (default: 3)
        Maximum number of traversed edges.
    weights: Sequence[float], optional
        Walk weights indexed by path length minus one. If `None`, use
        `1 / path_length**2`.
    minimum_path_number: int (default: 1)
        Minimum number of contributing walks required for an interaction.
    threshold: float (default: 0.75)
        Minimum absolute normalized score.
    allow_bidirectional: bool (default: False)
        If `True`, allow both directions of a gene pair to be retained.

    Returns
    -------
    list of tuple
        Interactions as `(source, target, {"sign": sign})`.

    Raises
    ------
    ValueError
        If scoring or signed-interaction inference receives invalid
        parameters.
    """

    scores = interaction_scores_from_walks(
        graph=graph,
        genes=genes,
        max_depth=max_depth,
        weights=weights,
    )

    return infer_signed_interactions(
        scores=scores,
        genes=genes,
        minimum_path_number=minimum_path_number,
        threshold=threshold,
        allow_bidirectional=allow_bidirectional,
    )


def _default_weights(max_depth: int) -> List[float]:

    return [1 / (length**2) for length in range(1, max_depth + 1)]


def _edge_signs(
    graph: nx.Graph[Any],
    source: str,
    target: str,
) -> List[int]:

    edge_data = graph.get_edge_data(source, target)

    if edge_data is None:
        return []

    if graph.is_multigraph():
        data_iterable = edge_data.values()
    else:
        data_iterable = [edge_data]

    signs = []

    for data in data_iterable:
        sign = data.get("sign")

        if sign not in {-1, 1}:
            raise ValueError(
                "invalid interaction sign: expected -1 or 1 for edge "
                f"{source!r} -> {target!r} but received {sign!r}"
            )

        signs.append(sign)

    return signs


def _signed_successor_counts(
    graph: nx.Graph[Any],
    source: str,
) -> Tuple[Tuple[str, int, int], ...]:
    """Return signed parallel-edge counts for each direct successor."""

    successors = []
    for target in getattr(graph, "successors")(source):
        signs = _edge_signs(graph, source, target)
        positive = sum(sign == 1 for sign in signs)
        successors.append((target, positive, len(signs) - positive))

    return tuple(successors)


def _path_weight(
    path_length: int,
    weights: Sequence[float],
) -> float:

    try:
        return weights[path_length - 1]

    except IndexError:
        raise ValueError(
            "invalid argument value for 'weights': "
            f"expected at least {path_length} values but received {len(weights)}"
        )


def _walk_signs(
    graph: nx.Graph[Any],
    walk: Sequence[str],
) -> List[int]:

    edge_signs = [
        _edge_signs(graph, source, target) for source, target in zip(walk, walk[1:])
    ]

    if any(len(signs) == 0 for signs in edge_signs):
        return []

    signs = []

    for sign_combination in product(*edge_signs):
        path_sign = 1

        for sign in sign_combination:
            path_sign *= sign

        signs.append(path_sign)

    return signs
