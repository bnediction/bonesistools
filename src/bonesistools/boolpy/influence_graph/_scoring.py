#!/usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import (
    Any,
    Iterable,
    Mapping,
    Sequence,
)

import networkx as nx

from ._algorithms import walks_from


@dataclass(frozen=True)
class InteractionScore:
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


def _default_weights(max_depth: int) -> list[float]:
    return [1 / (length**2) for length in range(1, max_depth + 1)]


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


def _edge_signs(
    graph: nx.Graph[Any],
    source: str,
    target: str,
) -> list[int]:
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


def _walk_signs(
    graph: nx.Graph[Any],
    walk: Sequence[str],
) -> list[int]:
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


def interaction_scores_from_walks(
    graph: nx.Graph[Any],
    genes: Iterable[str] | None = None,
    max_depth: int = 3,
    weights: Sequence[float] | None = None,
) -> dict[str, dict[str, InteractionScore]]:
    """
    Estimate signed interaction scores from bounded directed walks.

    Walks are enumerated from each source gene up to `max_depth`. A walk
    contributes the product of its edge signs multiplied by a length-dependent
    weight. Parallel signed edges in a MultiDiGraph are expanded into all sign
    combinations.

    Parameters
    ----------
    graph: nx.Graph
        Directed signed graph.
    genes: Iterable[str], optional
        Genes used as sources and targets. If None, use all graph nodes.
    max_depth: int (default: 3)
        Maximum number of traversed edges.
    weights: Sequence[float], optional
        Walk weights indexed by path length minus one. If None, use
        `1 / path_length**2`.

    Returns
    -------
    dict
        Nested mapping `source -> target -> InteractionScore`.

    Raises
    ------
    ValueError
        If `max_depth` is not positive, `weights` is too short or an edge sign
        is invalid.
    """

    if max_depth < 1:
        raise ValueError(
            "invalid argument value for 'max_depth': "
            f"expected a positive integer but received {max_depth!r}"
        )

    weights = _default_weights(max_depth) if weights is None else weights
    genes = list(graph.nodes if genes is None else genes)
    gene_set = set(genes)

    scores: dict[str, dict[str, InteractionScore]] = {gene: {} for gene in genes}
    accumulators: dict[str, dict[str, dict[str, float]]] = {gene: {} for gene in genes}

    for source in genes:
        if source not in graph:
            continue

        for walk in walks_from(graph, source=source, max_depth=max_depth):
            target = walk[-1]

            if target == source:
                continue

            if target not in gene_set:
                continue

            path_length = len(walk) - 1
            weight = _path_weight(path_length, weights)

            for path_sign in _walk_signs(graph, walk):
                target_scores = accumulators[source].setdefault(
                    target,
                    {"score": 0.0, "total_weight": 0.0, "path_number": 0},
                )
                target_scores["score"] += path_sign * weight
                target_scores["total_weight"] += weight
                target_scores["path_number"] += 1

    for source, targets in accumulators.items():
        for target, values in targets.items():
            scores[source][target] = InteractionScore(
                score=values["score"],
                total_weight=values["total_weight"],
                path_number=int(values["path_number"]),
            )

    return scores


def _empty_score() -> InteractionScore:
    return InteractionScore(score=0.0, total_weight=0.0, path_number=0)


def _passes_threshold(
    score: InteractionScore,
    minimum_path_number: int,
    threshold: float,
) -> bool:
    return (
        score.path_number >= minimum_path_number
        and score.normalized_score != 0
        and abs(score.normalized_score) >= threshold
    )


def infer_signed_interactions(
    scores: Mapping[str, Mapping[str, InteractionScore]],
    genes: Iterable[str] | None = None,
    minimum_path_number: int = 1,
    threshold: float = 0.75,
    allow_bidirectional: bool = False,
) -> list[tuple[str, str, dict[str, int]]]:
    """
    Infer signed directed interactions from normalized interaction scores.

    Parameters
    ----------
    scores: Mapping[str, Mapping[str, InteractionScore]]
        Nested interaction scores.
    genes: Iterable[str], optional
        Genes to compare. If None, use keys from `scores`.
    minimum_path_number: int (default: 1)
        Minimum number of contributing walks required for an interaction.
    threshold: float (default: 0.75)
        Minimum absolute normalized score.
    allow_bidirectional: bool (default: False)
        If True, allow both directions of a gene pair to be retained.

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

    if minimum_path_number < 1:
        raise ValueError(
            "invalid argument value for 'minimum_path_number': "
            f"expected a positive integer but received {minimum_path_number!r}"
        )

    if not 0 <= threshold <= 1:
        raise ValueError(
            "invalid argument value for 'threshold': "
            f"expected a value in [0, 1] but received {threshold!r}"
        )

    genes = list(scores if genes is None else genes)
    interactions = []

    for idx, source in enumerate(genes):
        for target in genes[idx + 1 :]:
            forward = scores.get(source, {}).get(target, _empty_score())
            backward = scores.get(target, {}).get(source, _empty_score())

            candidates = [
                (source, target, forward, backward),
                (target, source, backward, forward),
            ]

            for candidate_source, candidate_target, score, opposite in candidates:
                if not _passes_threshold(score, minimum_path_number, threshold):
                    continue

                if not allow_bidirectional and abs(score.normalized_score) <= abs(
                    opposite.normalized_score
                ):
                    continue

                sign = 1 if score.normalized_score > 0 else -1
                interactions.append(
                    (candidate_source, candidate_target, {"sign": sign})
                )

    return interactions


def infer_signed_interactions_from_walks(
    graph: nx.Graph[Any],
    genes: Iterable[str] | None = None,
    max_depth: int = 3,
    weights: Sequence[float] | None = None,
    minimum_path_number: int = 1,
    threshold: float = 0.75,
    allow_bidirectional: bool = False,
) -> list[tuple[str, str, dict[str, int]]]:
    """
    Infer signed interactions directly from bounded directed walks.

    This convenience function combines `interaction_scores_from_walks()` and
    `infer_signed_interactions()`. Use `interaction_scores_from_walks()`
    directly when intermediate `InteractionScore` objects should be inspected.

    Parameters
    ----------
    graph: nx.Graph
        Directed signed graph.
    genes: Iterable[str], optional
        Genes used as sources and targets. If None, use all graph nodes.
    max_depth: int (default: 3)
        Maximum number of traversed edges.
    weights: Sequence[float], optional
        Walk weights indexed by path length minus one. If None, use
        `1 / path_length**2`.
    minimum_path_number: int (default: 1)
        Minimum number of contributing walks required for an interaction.
    threshold: float (default: 0.75)
        Minimum absolute normalized score.
    allow_bidirectional: bool (default: False)
        If True, allow both directions of a gene pair to be retained.

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
