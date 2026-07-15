#!/usr/bin/env python

from typing import Any, cast

import networkx as nx
import pytest

from bonesistools.logic.influence_graph import (
    AggregatedInfluenceGraph,
    InfluenceGraph,
)


def _edge_counts(graph):
    return sorted(
        (source, target, data["sign"], data["count"])
        for source, target, data in graph.edges(data=True)
    )


def _single_edge_graph(total=4):
    graph = AggregatedInfluenceGraph(total=total)
    graph.add_edge("A", "B", sign=1, count=3)

    return graph


def _collapse_graph():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edges_from(
        [
            ("TF", "g1", {"sign": 1, "count": 3}),
            ("TF", "g2", {"sign": 1, "count": 3}),
            ("g1", "out", {"sign": -1, "count": 2}),
            ("g2", "out", {"sign": -1, "count": 2}),
        ]
    )

    return graph


def test_doc_example_constructor_edge_frequency_and_autoregulations():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=3)
    graph.add_edge("B", "C", sign=-1, count=1)

    assert graph.edge_frequency("A", "B") == 0.75
    assert graph.autoregulations() == []


def test_aggregated_influence_graph_is_public():
    from bonesistools import logic

    assert logic.ig.AggregatedInfluenceGraph is AggregatedInfluenceGraph


def test_doc_example_copy():
    graph = _single_edge_graph()

    copied = graph.copy()

    assert copied.total == 4
    assert copied.edge_count("A", "B") == 3

    with pytest.raises(NotImplementedError, match="view copies"):
        graph.copy(as_view=True)


def test_doc_example_add_edge_allows_distinct_signed_edges():
    graph = AggregatedInfluenceGraph(total=4)

    assert graph.add_edge("A", "B", sign=1, count=3) is None
    assert graph.add_edge("A", "B", sign=-1, count=1) is None
    assert _edge_counts(graph) == [("A", "B", -1, 1), ("A", "B", 1, 3)]

    key_style = AggregatedInfluenceGraph(total=2)
    key_style.add_edge("A", "B", 1, count=1)

    assert _edge_counts(key_style) == [("A", "B", 1, 1)]


def test_doc_example_add_edges_from_and_update():
    graph = AggregatedInfluenceGraph(total=3)
    graph.add_edges_from(
        [
            ("A", "B", {"sign": 1, "count": 2}),
            ("B", "C", {"sign": -1, "count": 1}),
        ]
    )

    assert _edge_counts(graph) == [("A", "B", 1, 2), ("B", "C", -1, 1)]

    with pytest.raises(ValueError, match="invalid edge count"):
        graph.add_edges_from([("C", "D", {"sign": 1, "count": 4})])

    assert "D" not in graph

    graph = AggregatedInfluenceGraph(total=2)
    graph.add_edge("A", "B", sign=1, count=1)
    graph.update(edges=[("B", "C", {"sign": -1, "count": 2})])

    assert graph.edge_count("B", "C") == 2


def test_doc_example_total_accessors_and_edge_frequency():
    graph = _single_edge_graph()

    assert graph.total == 4

    graph.total = 6

    assert graph.total == 6
    assert graph.edge_frequency("A", "B") == 0.5

    with pytest.raises(ValueError, match="invalid edge count"):
        graph.total = 2


def test_doc_example_edge_count_with_ambiguous_signed_edges():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=3)

    assert graph.edge_count("A", "B") == 3

    graph.add_edge("A", "B", sign=-1, count=1)

    assert graph.edge_count("A", "B", sign=-1) == 1


def test_doc_example_autoregulations_and_frequency_bins():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "A", sign=1, count=3)
    graph.add_edge("B", "B", sign=-1, count=1)

    assert graph.autoregulations() == [("A", 1, 0.75), ("B", -1, 0.25)]
    assert graph.autoregulations(sign=1) == [("A", 1, 0.75)]
    assert graph.frequency_bin("A", "A", bins=(0.0, 0.5, 0.75, 1.0)) == (
        0.5,
        0.75,
    )

    with pytest.raises(ValueError, match="cover"):
        graph.frequency_bin("A", "A", bins=(0.0, 0.5))


def test_deprecated_influence_graph_wrappers_are_disabled():
    graph = _collapse_graph()

    with pytest.raises(NotImplementedError, match="family_collapsed_graph"):
        graph.family_compressed_graph(preserve_feedback=False)

    with pytest.raises(NotImplementedError, match="collapsed_graph"):
        graph.compressed_graph()


def test_doc_example_validate_counts_rejects_invalid_count():
    invalid = nx.MultiDiGraph()
    invalid.add_edge("A", "B", sign=1, count=4)

    with pytest.raises(ValueError, match="invalid edge count"):
        AggregatedInfluenceGraph(invalid, total=3)


def test_validation_accepts_integer_float_count_and_rejects_fractional_count():
    missing_sign = nx.MultiDiGraph()
    missing_sign.add_edge("A", "B", count=1)

    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        AggregatedInfluenceGraph(missing_sign, total=2)

    graph = AggregatedInfluenceGraph(total=2)

    graph.add_edge("A", "B", sign=1, count=cast(Any, 1.0))
    assert cast(Any, graph["A"]["B"])[0]["count"] == 1

    with pytest.raises(ValueError, match="expected integer"):
        graph.add_edge("A", "B", sign=1, count=cast(Any, 1.5))


def test_invariant_validation_detects_low_level_graph_mutations():
    graph = AggregatedInfluenceGraph(total=2)
    graph.add_edge("A", "B", sign=1, count=1)

    nx.MultiDiGraph.add_edge(graph, "B", "C", sign=1, count=3)

    with pytest.raises(ValueError, match="invalid edge count"):
        graph.validate_counts()

    graph = AggregatedInfluenceGraph(total=2)
    nx.MultiDiGraph.add_edge(graph, "A", "B", sign=1)

    with pytest.raises(KeyError, match="missing required edge attribute 'count'"):
        graph.validate_counts()

    graph = AggregatedInfluenceGraph(total=2)
    nx.MultiDiGraph.add_edge(graph, "A", "B", count=1)

    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        graph._validate_graph()

    graph = AggregatedInfluenceGraph(total=2)
    nx.MultiDiGraph.add_edge(graph, "A", "B", sign=1, count=1)
    nx.MultiDiGraph.add_edge(graph, "A", "B", sign=1, count=1)

    with pytest.raises(ValueError, match="duplicated edge sign"):
        graph._validate_graph()


def test_from_influence_graphs_counts_signed_edges_separately():
    ig1 = InfluenceGraph()
    ig1.add_edge("A", "B", sign=1)
    ig1.add_edge("A", "B", sign=-1)
    ig1.add_node("isolated")

    ig2 = InfluenceGraph()
    ig2.add_edge("A", "B", sign=1)
    ig2.add_edge("B", "C", sign=-1)

    aggregated = AggregatedInfluenceGraph.from_influence_graphs(ig1, ig2)

    assert aggregated.total == 2
    assert sorted(aggregated.nodes()) == ["A", "B", "C", "isolated"]
    assert _edge_counts(aggregated) == [
        ("A", "B", -1, 1),
        ("A", "B", 1, 2),
        ("B", "C", -1, 1),
    ]
    assert aggregated.edge_count("A", "B", sign=1) == 2
    assert aggregated.edge_count("A", "B", sign=-1) == 1
    assert aggregated.edge_frequency("A", "B", sign=1) == 1.0

    with pytest.raises(ValueError, match="ambiguous aggregated edge"):
        aggregated.edge_count("A", "B")


def test_from_influence_graphs_rejects_empty_input():
    with pytest.raises(ValueError, match="at least one influence graph"):
        AggregatedInfluenceGraph.from_influence_graphs()


def test_total_is_changed_only_through_validated_setter():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=3)

    assert graph.total == 4

    graph.total = 10

    assert graph.total == 10

    with pytest.raises(TypeError, match="unsupported argument type"):
        graph.total = True

    with pytest.raises(ValueError, match="positive value"):
        graph.total = 0

    with pytest.raises(ValueError, match="invalid edge count"):
        graph.total = 2

    graph.total = 6

    assert graph.total == 6
    assert graph.total == 6
    assert graph.edge_frequency("A", "B") == 0.5


def test_constructor_and_mutations_validate_counts_transactionally():
    missing_count = nx.MultiDiGraph()
    missing_count.add_edge("A", "B", sign=1)

    with pytest.raises(KeyError, match="missing required edge attribute 'count'"):
        AggregatedInfluenceGraph(missing_count, total=2)

    graph = AggregatedInfluenceGraph(total=2)
    graph.add_edge("A", "B", sign=1, count=1)

    with pytest.raises(TypeError, match="missing required argument: 'count'"):
        graph.add_edge("B", "C", sign=1)

    with pytest.raises(ValueError, match="invalid edge count"):
        graph.add_edges_from([("B", "C", {"sign": 1, "count": 3})])

    assert _edge_counts(graph) == [("A", "B", 1, 1)]
    assert "C" not in graph
