#!/usr/bin/env python

from typing import Iterable, Tuple

import pytest

import bonesistools as bt
from bonesistools.logic.influence_graph._influence_graph import InfluenceSign


def _influence_graph(edges: Iterable[Tuple[str, str, InfluenceSign]]):

    graph = bt.logic.ig.InfluenceGraph()

    for source, target, sign in edges:
        graph.add_edge(source, target, sign=sign)

    return graph


def test_identical_signed_graphs_have_full_jaccard_similarity():

    graph = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )

    assert bt.logic.ig.similarity(graph, graph) == 1.0
    assert bt.logic.ig.distance(graph, graph) == 0.0


def test_disjoint_signed_graphs_have_zero_jaccard_similarity():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("C", "D", -1)])

    assert bt.logic.ig.similarity(graph1, graph2) == 0.0
    assert bt.logic.ig.distance(graph1, graph2) == 1.0


def test_signed_jaccard_similarity_counts_partial_overlap():

    graph1 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("C", "D", 1),
        ]
    )

    assert bt.logic.ig.similarity(graph1, graph2) == pytest.approx(1 / 3)
    assert bt.logic.ig.distance(graph1, graph2) == pytest.approx(2 / 3)


def test_signed_hamming_uses_complete_node_union_domain_by_default():

    graph1 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("C", "D", 1),
        ]
    )

    assert bt.logic.ig.similarity(graph1, graph2, metric="hamming") == pytest.approx(
        15 / 16
    )
    assert bt.logic.ig.distance(
        graph1, graph2, metric="hamming"
    ) == pytest.approx(1 / 16)


def test_signed_hamming_uses_explicit_domain():

    graph1 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("C", "D", 1),
        ]
    )
    dorothea_like_domain = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
            ("C", "D", 1),
            ("D", "A", -1),
        ]
    )

    assert bt.logic.ig.similarity(
        graph1,
        graph2,
        metric="hamming",
        domain=dorothea_like_domain,
    ) == pytest.approx(1 / 2)
    assert bt.logic.ig.distance(
        graph1,
        graph2,
        metric="hamming",
        domain=dorothea_like_domain,
    ) == pytest.approx(1 / 2)


def test_signed_jaccard_similarity_distinguishes_opposite_signs():

    positive = _influence_graph([("A", "B", 1)])
    negative = _influence_graph([("A", "B", -1)])

    assert bt.logic.ig.similarity(positive, negative) == 0.0
    assert bt.logic.ig.distance(positive, negative) == 1.0


def test_signed_hamming_distinguishes_opposite_signs():

    positive = _influence_graph([("A", "B", 1)])
    negative = _influence_graph([("A", "B", -1)])
    domain = _influence_graph(
        [
            ("A", "B", 1),
            ("A", "B", -1),
        ]
    )

    assert (
        bt.logic.ig.similarity(
            positive,
            negative,
            metric="hamming",
            domain=domain,
        )
        == 0.0
    )
    assert (
        bt.logic.ig.distance(
            positive,
            negative,
            metric="hamming",
            domain=domain,
        )
        == 1.0
    )


def test_signed_jaccard_similarity_ignores_edge_attributes():

    graph1 = bt.logic.ig.InfluenceGraph()
    graph1.add_edge("A", "B", sign=1, evidence="curated")

    graph2 = bt.logic.ig.InfluenceGraph()
    graph2.add_edge("A", "B", sign=1, confidence=0.7)

    assert bt.logic.ig.similarity(graph1, graph2) == 1.0
    assert bt.logic.ig.distance(graph1, graph2) == 0.0


def test_empty_graphs_have_full_jaccard_similarity():

    graph1 = bt.logic.ig.InfluenceGraph()
    graph2 = bt.logic.ig.InfluenceGraph()

    assert bt.logic.ig.similarity(graph1, graph2) == 1.0
    assert bt.logic.ig.distance(graph1, graph2) == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, metric="hamming") == 1.0
    assert bt.logic.ig.distance(graph1, graph2, metric="hamming") == 0.0


def test_hamming_rejects_non_strict_or_incomplete_explicit_domain():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("C", "D", 1)])

    with pytest.raises(ValueError, match="strict signed-edge superset"):
        bt.logic.ig.distance(
            graph1,
            graph2,
            metric="hamming",
            domain=graph1,
        )

    incomplete_domain = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    with pytest.raises(ValueError, match="strict signed-edge superset"):
        bt.logic.ig.distance(
            graph1,
            graph2,
            metric="hamming",
            domain=incomplete_domain,
        )


def test_unsupported_similarity_and_distance_metrics_raise_clear_errors():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("A", "B", 1)])

    with pytest.raises(ValueError, match="unsupported similarity metric"):
        bt.logic.ig.similarity(graph1, graph2, metric="cosine")

    with pytest.raises(ValueError, match="unsupported distance metric"):
        bt.logic.ig.distance(graph1, graph2, metric="cosine")
