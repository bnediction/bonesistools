#!/usr/bin/env python

from typing import Any, Iterable, Tuple, cast

import pytest

import bonesistools as bt
from bonesistools.logic.influence_graph._influence_graph import InfluenceSign


def _influence_graph(
    edges: Iterable[Tuple[str, str, InfluenceSign]] = (),
    *,
    nodes: Iterable[str] = (),
):

    graph = bt.logic.ig.InfluenceGraph()
    graph.add_nodes_from(nodes)

    for source, target, sign in edges:
        graph.add_edge(source, target, sign=sign)

    return graph


def test_union_universe_identical_graphs_have_distance_zero():

    graph = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )

    assert bt.logic.ig.distance(graph, graph) == 0.0
    assert bt.logic.ig.similarity(graph, graph) == 1.0


def test_union_universe_values_avoid_complement_rounding():

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

    assert bt.logic.ig.similarity(graph1, graph2) == 1 / 3
    assert bt.logic.ig.distance(graph1, graph2) == 2 / 3


def test_union_universe_disjoint_graphs_have_distance_one():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("C", "D", -1)])

    assert bt.logic.ig.distance(graph1, graph2) == 1.0
    assert bt.logic.ig.similarity(graph1, graph2) == 0.0


def test_union_universe_empty_graphs_are_identical():

    graph1 = _influence_graph(nodes=["A"])
    graph2 = _influence_graph(nodes=["B"])

    assert bt.logic.ig.distance(graph1, graph2) == 0.0
    assert bt.logic.ig.similarity(graph1, graph2) == 1.0


def test_union_universe_distinguishes_opposite_signs():

    positive = _influence_graph([("A", "B", 1)])
    negative = _influence_graph([("A", "B", -1)])

    assert bt.logic.ig.distance(positive, negative) == 1.0
    assert bt.logic.ig.similarity(positive, negative) == 0.0


def test_union_universe_treats_self_loops_as_observed_edges():

    graph1 = _influence_graph([("A", "A", 1)])
    graph2 = _influence_graph(nodes=["A"])

    assert bt.logic.ig.distance(graph1, graph2) == 1.0
    assert bt.logic.ig.similarity(graph1, graph2) == 0.0


def test_signed_edge_comparison_ignores_other_attributes():

    graph1 = bt.logic.ig.InfluenceGraph()
    graph1.add_edge("A", "B", sign=1, evidence="curated")

    graph2 = bt.logic.ig.InfluenceGraph()
    graph2.add_edge("A", "B", sign=1, confidence=0.7)

    assert bt.logic.ig.distance(graph1, graph2) == 0.0
    assert bt.logic.ig.similarity(graph1, graph2) == 1.0


def test_complete_universe_has_two_signed_positions_per_node_pair():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "A", -1),
        ]
    )

    assert bt.logic.ig.distance(graph1, graph2, universe="complete") == 1 / 8
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 7 / 8


def test_complete_universe_distinguishes_positive_and_negative_positions():

    positive = _influence_graph([("A", "B", 1)])
    negative = _influence_graph([("A", "B", -1)])

    assert bt.logic.ig.distance(positive, negative, universe="complete") == 1 / 4
    assert bt.logic.ig.similarity(positive, negative, universe="complete") == 3 / 4


def test_complete_universe_counts_joint_absences_as_agreements():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "A", -1),
        ]
    )

    assert bt.logic.ig.similarity(graph1, graph2, universe="union") == 1 / 2
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 7 / 8


def test_complete_universe_includes_signed_self_loops():

    graph1 = _influence_graph([("A", "A", 1)])
    graph2 = _influence_graph(nodes=["A"])

    assert bt.logic.ig.distance(graph1, graph2, universe="complete") == 1 / 2
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 1 / 2


def test_complete_universe_empty_graphs_with_nodes_are_identical():

    graph1 = _influence_graph(nodes=["A", "B"])
    graph2 = _influence_graph(nodes=["A", "B"])

    assert bt.logic.ig.distance(graph1, graph2, universe="complete") == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 1.0


def test_complete_universe_graphs_without_nodes_are_identical():

    graph1 = _influence_graph()
    graph2 = _influence_graph()

    assert bt.logic.ig.distance(graph1, graph2, universe="complete") == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 1.0


def test_complete_universe_values_avoid_complement_rounding():

    nodes = ("A", "B", "C")
    complete_edges = [
        (source, target, sign)
        for source in nodes
        for target in nodes
        for sign in (1, -1)
    ]
    graph1 = _influence_graph(complete_edges[:6], nodes=nodes)
    graph2 = _influence_graph(complete_edges[6:12], nodes=nodes)

    assert bt.logic.ig.similarity(graph1, graph2, universe="complete") == 1 / 3
    assert bt.logic.ig.distance(graph1, graph2, universe="complete") == 2 / 3


def test_complete_loopless_universe_excludes_signed_self_loops():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "A", -1),
        ]
    )

    assert bt.logic.ig.distance(graph1, graph2, universe="complete_loopless") == 1 / 4
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete_loopless") == 3 / 4


@pytest.mark.parametrize(
    "comparison",
    [bt.logic.ig.distance, bt.logic.ig.similarity],
)
def test_complete_loopless_universe_rejects_input_self_loops(comparison):

    graph1 = _influence_graph([("A", "A", 1)])
    graph2 = _influence_graph(nodes=["A"])

    with pytest.raises(ValueError, match="cannot contain input self-loops"):
        comparison(graph1, graph2, universe="complete_loopless")


def test_one_node_complete_loopless_universe_is_empty():

    graph1 = _influence_graph(nodes=["A"])
    graph2 = _influence_graph(nodes=["A"])

    assert bt.logic.ig.distance(graph1, graph2, universe="complete_loopless") == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, universe="complete_loopless") == 1.0


def test_explicit_universe_accepts_exact_input_edge_set():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("A", "B", 1)])
    universe = graph1.copy()

    assert bt.logic.ig.distance(graph1, graph2, universe=universe) == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, universe=universe) == 1.0


def test_explicit_universe_accepts_non_strict_superset():

    graph1 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    graph2 = _influence_graph([("A", "B", 1)])
    universe = graph1.copy()

    assert bt.logic.ig.distance(graph1, graph2, universe=universe) == 1 / 2
    assert bt.logic.ig.similarity(graph1, graph2, universe=universe) == 1 / 2


@pytest.mark.parametrize(
    "comparison",
    [bt.logic.ig.distance, bt.logic.ig.similarity],
)
@pytest.mark.parametrize("sign", [1, -1])
def test_explicit_universe_rejects_missing_signed_edge(comparison, sign):

    graph1 = _influence_graph([("A", "B", sign)])
    graph2 = _influence_graph(nodes=["A", "B"])
    universe = _influence_graph()

    with pytest.raises(ValueError, match="must contain every signed edge"):
        comparison(graph1, graph2, universe=universe)


def test_empty_explicit_universe_accepts_empty_graphs():

    graph1 = _influence_graph(nodes=["A"])
    graph2 = _influence_graph(nodes=["B"])
    universe = _influence_graph()

    assert bt.logic.ig.distance(graph1, graph2, universe=universe) == 0.0
    assert bt.logic.ig.similarity(graph1, graph2, universe=universe) == 1.0


def test_explicit_universe_distinguishes_opposite_signs():

    positive = _influence_graph([("A", "B", 1)])
    negative = _influence_graph([("A", "B", -1)])
    universe = _influence_graph(
        [
            ("A", "B", 1),
            ("A", "B", -1),
        ]
    )

    assert bt.logic.ig.distance(positive, negative, universe=universe) == 1.0
    assert bt.logic.ig.similarity(positive, negative, universe=universe) == 0.0


def test_distance_and_similarity_are_normalized_for_every_universe():

    graph1 = _influence_graph(
        [
            ("A", "B", 1),
            ("B", "C", -1),
        ]
    )
    graph2 = _influence_graph(
        [
            ("A", "B", 1),
            ("C", "B", 1),
        ]
    )
    explicit_universe = _influence_graph(
        [
            ("A", "A", 1),
            ("A", "B", 1),
            ("B", "C", -1),
            ("C", "B", 1),
        ]
    )

    for universe in (
        "union",
        "complete",
        "complete_loopless",
        explicit_universe,
    ):
        distance = bt.logic.ig.distance(graph1, graph2, universe=universe)
        similarity = bt.logic.ig.similarity(graph1, graph2, universe=universe)

        assert 0.0 <= distance <= 1.0
        assert 0.0 <= similarity <= 1.0
        assert distance + similarity == pytest.approx(1.0)


@pytest.mark.parametrize(
    "comparison",
    [bt.logic.ig.distance, bt.logic.ig.similarity],
)
def test_unsupported_edge_universe_raises_clear_error(comparison):

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("A", "B", 1)])

    with pytest.raises(ValueError, match="unsupported edge universe"):
        comparison(graph1, graph2, universe=cast(Any, "observed"))


def test_removed_metric_and_domain_arguments_are_rejected():

    graph1 = _influence_graph([("A", "B", 1)])
    graph2 = _influence_graph([("A", "B", 1)])
    prior = graph1.copy()
    distance = cast(Any, bt.logic.ig.distance)
    similarity = cast(Any, bt.logic.ig.similarity)

    with pytest.raises(TypeError):
        distance(graph1, graph2, metric="jaccard")

    with pytest.raises(TypeError):
        similarity(graph1, graph2, domain=prior)
