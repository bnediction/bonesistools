#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, cast

import networkx as nx
import pytest

from bonesistools.boolpy.boolean_network import BooleanNetwork, BooleanNetworkEnsemble
from bonesistools.boolpy.influence_graph import (
    AggregatedInfluenceGraph,
    InfluenceGraph,
)


def _edge_counts(graph):
    return sorted(
        (source, target, data["sign"], data["count"])
        for source, target, data in graph.edges(data=True)
    )


def _rendered_edges(graph):
    return sorted(
        (source, target, attrs["label"], attrs["arrowhead"])
        for source, target, attrs in graph.edges
    )


def _pydot_get_string(obj, method):
    return cast(str, getattr(cast(Any, obj), method)()).strip('"')


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


def _feedback_graph():
    graph = AggregatedInfluenceGraph(total=3)
    graph.add_edges_from(
        [
            ("A", "B", {"sign": 1, "count": 2}),
            ("B", "A", {"sign": 1, "count": 2}),
            ("B", "C", {"sign": -1, "count": 1}),
        ]
    )

    return graph


def _three_boolean_networks():
    return (
        BooleanNetwork({"A": "B", "B": 0, "C": 1}),
        BooleanNetwork({"A": "B", "B": "C", "C": 0}),
        BooleanNetwork({"A": "!B", "B": "C", "C": 1}),
    )


def test_doc_example_constructor_edge_frequency_and_autoregulations():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=3)
    graph.add_edge("B", "C", sign=-1, count=1)

    assert graph.edge_frequency("A", "B") == 0.75
    assert graph.autoregulations() == []


def test_aggregated_influence_graph_is_public():
    from bonesistools import bpy

    assert bpy.ig.AggregatedInfluenceGraph is AggregatedInfluenceGraph


def test_doc_example_copy():
    graph = _single_edge_graph()

    copied = graph.copy()

    assert copied.total == 4
    assert copied.edge_count("A", "B") == 3

    with pytest.raises(NotImplementedError, match="view copies"):
        graph.copy(as_view=True)


def test_doc_example_add_edge_allows_distinct_signed_edges():
    graph = AggregatedInfluenceGraph(total=4)

    assert graph.add_edge("A", "B", sign=1, count=3) == 0
    assert graph.add_edge("A", "B", sign=-1, count=1) == 1
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

    with pytest.raises(ValueError, match="not covered by bins"):
        graph.frequency_bin("A", "A", bins=(0.0, 0.5))


def test_doc_example_structural_families():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edges_from(
        [
            ("TF", "g1", {"sign": 1, "count": 3}),
            ("TF", "g2", {"sign": 1, "count": 3}),
            ("TF", "g3", {"sign": 1, "count": 1}),
        ]
    )

    families = graph.structural_families(
        include_successors=False,
        bins=(0.0, 0.5, 1.0),
    )

    assert set(map(frozenset, families.values())) == {frozenset({"g1", "g2"})}


def test_doc_example_family_and_plain_collapsed_graphs():
    graph = _collapse_graph()

    family_graph = graph.family_collapsed_graph(exclude_feedback_nodes=False)

    assert isinstance(family_graph, InfluenceGraph)
    assert not isinstance(family_graph, AggregatedInfluenceGraph)
    assert sorted(family_graph.nodes()) == ["TF", "g1|g2", "out"]
    assert family_graph.nodes["g1|g2"]["members"] == {"g1", "g2"}
    assert cast(Any, family_graph["TF"]["g1|g2"])[0]["frequency"] == 0.75


def test_family_collapsed_graph_uses_pipe_separator_by_default():
    graph = _collapse_graph()

    collapsed = graph.family_collapsed_graph(exclude_feedback_nodes=False)

    assert "g1|g2" in collapsed
    assert "g1/g2" not in collapsed


def test_family_collapsed_graph_accepts_custom_separator():
    graph = _collapse_graph()

    collapsed = graph.family_collapsed_graph(
        exclude_feedback_nodes=False,
        sep="/",
    )

    assert "g1/g2" in collapsed
    assert "g1|g2" not in collapsed


def test_doc_example_feedback_induced_graph():
    graph = _feedback_graph()

    feedback = graph.feedback_induced_graph()

    assert sorted(feedback.nodes()) == ["A", "B"]
    assert feedback.total == 3


def test_doc_example_collapsed_graph():
    graph = _feedback_graph()

    collapsed = graph.collapsed_graph()

    assert isinstance(collapsed, InfluenceGraph)
    assert not isinstance(collapsed, AggregatedInfluenceGraph)
    assert sorted(collapsed.nodes()) == ["A", "B"]


def test_deprecated_influence_graph_wrappers_are_disabled():
    graph = _collapse_graph()

    with pytest.raises(NotImplementedError, match="family_collapsed_graph"):
        graph.family_compressed_graph(exclude_feedback_nodes=False)

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


def test_structural_families_can_protect_feedback_nodes():
    graph = AggregatedInfluenceGraph(total=3)
    graph.add_edges_from(
        [
            ("A", "B", {"sign": 1, "count": 2}),
            ("B", "A", {"sign": 1, "count": 2}),
            ("TF", "g1", {"sign": 1, "count": 2}),
            ("TF", "g2", {"sign": 1, "count": 2}),
        ]
    )

    families = graph.structural_families(
        include_successors=False,
        protect_feedback_nodes=True,
    )

    assert set(map(frozenset, families.values())) == {frozenset({"g1", "g2"})}
    assert all("A" not in family and "B" not in family for family in families.values())


def test_family_collapsed_graph_ignores_singleton_families():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("TF", "g1", sign=1, count=3)

    collapsed = graph.family_collapsed_graph(exclude_feedback_nodes=False)

    assert sorted(collapsed.nodes()) == ["TF", "g1"]
    assert sorted(collapsed.edges()) == [("TF", "g1")]


def test_family_collapse_does_not_merge_nodes_with_internal_edges():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edges_from(
        [
            ("TF", "g1", {"sign": 1, "count": 3}),
            ("TF", "g2", {"sign": 1, "count": 3}),
            ("g1", "g2", {"sign": 1, "count": 2}),
            ("g2", "g1", {"sign": 1, "count": 2}),
        ]
    )

    collapsed = graph.family_collapsed_graph(
        include_successors=False,
        exclude_feedback_nodes=False,
    )

    assert sorted(collapsed.nodes()) == ["TF", "g1", "g2"]
    assert sorted(collapsed.edges()) == [
        ("TF", "g1"),
        ("TF", "g2"),
        ("g1", "g2"),
        ("g2", "g1"),
    ]


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


def test_copy_and_collapse_preserve_expected_types():
    graph = _collapse_graph()

    copied = graph.copy()
    collapsed = graph.family_collapsed_graph(exclude_feedback_nodes=False)

    assert isinstance(copied, AggregatedInfluenceGraph)
    assert copied.total == 4
    assert isinstance(collapsed, InfluenceGraph)
    assert sorted(collapsed.nodes()) == ["TF", "g1|g2", "out"]
    assert cast(Any, collapsed["TF"]["g1|g2"])[0]["frequency"] == 0.75
    assert cast(Any, collapsed["g1|g2"]["out"])[0]["frequency"] == 0.5


def test_from_boolean_networks_accepts_varargs_and_ensemble():
    bn1, bn2, bn3 = _three_boolean_networks()

    aggregated = AggregatedInfluenceGraph.from_boolean_networks(
        bn1,
        bn2,
        bn3,
    )

    assert aggregated.total == 3
    assert aggregated.edge_count("B", "A", sign=1) == 2
    assert aggregated.edge_count("B", "A", sign=-1) == 1
    assert aggregated.edge_count("C", "B") == 2
    assert aggregated.nodes["A"]["function_count"] == 2
    assert aggregated.nodes["A"]["function_stability"] == pytest.approx(2 / 3)

    ensemble = BooleanNetworkEnsemble(bns=[bn1, bn2, bn3])
    from_ensemble = AggregatedInfluenceGraph.from_boolean_networks(ensemble)

    assert from_ensemble.total == 3
    assert from_ensemble.edge_count("B", "A", sign=1) == 2
    assert from_ensemble.edge_count("B", "A", sign=-1) == 1
    assert from_ensemble.edge_count("C", "B") == 2
    assert from_ensemble.nodes["A"]["function_count"] == 2
    assert from_ensemble.nodes["A"]["function_stability"] == pytest.approx(2 / 3)


def test_to_graphviz_supports_collapse_modes(fake_graphviz):
    graph = _collapse_graph()

    exact = graph.to_graphviz(rankdir="LR")
    family = graph.to_graphviz(collapse="family")

    assert isinstance(exact, fake_graphviz)
    assert exact.graph_attr["rankdir"] == "LR"
    assert _rendered_edges(exact) == [
        ("TF", "g1", "3", "normal"),
        ("TF", "g2", "3", "normal"),
        ("g1", "out", "2", "tee"),
        ("g2", "out", "2", "tee"),
    ]
    assert _rendered_edges(family) == [
        ("TF", "g1|g2", "3", "normal"),
        ("g1|g2", "out", "2", "tee"),
    ]

    feedback_graph = AggregatedInfluenceGraph(total=4)
    feedback_graph.add_edges_from(
        [
            ("A", "B", {"sign": 1, "count": 3}),
            ("B", "A", {"sign": 1, "count": 3}),
            ("B", "C", {"sign": -1, "count": 1}),
        ]
    )

    feedback = feedback_graph.to_graphviz(collapse="feedback")
    both = feedback_graph.to_graphviz(collapse="both")

    assert sorted((source, target) for source, target, _ in feedback.edges) == [
        ("A", "B"),
        ("B", "A"),
    ]
    assert _rendered_edges(feedback) == [
        ("A", "B", "3", "normal"),
        ("B", "A", "3", "normal"),
    ]
    assert sorted((source, target) for source, target, _ in both.edges) == [
        ("A", "B"),
        ("B", "A"),
    ]
    assert _rendered_edges(both) == [
        ("A", "B", "3", "normal"),
        ("B", "A", "3", "normal"),
    ]

    with pytest.raises(ValueError, match="unsupported collapse"):
        graph.to_graphviz(collapse=cast(Any, "bad"))

    with pytest.raises(ValueError, match="min_frequency"):
        graph.to_graphviz(min_frequency=2)


def test_to_graphviz_can_filter_and_disable_frequency_style(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(
        collapse="family",
        min_frequency=0.6,
        edge_label=None,
        edge_style=False,
    )

    assert isinstance(rendered, fake_graphviz)
    assert sorted((source, target) for source, target, _ in rendered.edges) == [
        ("TF", "g1|g2"),
    ]
    assert all("label" not in attrs for _, _, attrs in rendered.edges)
    assert all("style" not in attrs for _, _, attrs in rendered.edges)


def test_to_graphviz_can_force_frequency_edge_labels_on_exact_graph(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(edge_label="frequency")

    assert isinstance(rendered, fake_graphviz)
    assert _rendered_edges(rendered) == [
        ("TF", "g1", "0.75", "normal"),
        ("TF", "g2", "0.75", "normal"),
        ("g1", "out", "0.5", "tee"),
        ("g2", "out", "0.5", "tee"),
    ]


def test_to_graphviz_can_disable_edge_labels_with_edge_label_none(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(edge_label=None)

    assert isinstance(rendered, fake_graphviz)
    assert all("label" not in attrs for _, _, attrs in rendered.edges)


def test_to_graphviz_deprecates_show_edge_labels_without_disabling_labels(
    fake_graphviz,
):
    graph = _single_edge_graph()

    with pytest.warns(FutureWarning, match="show_edge_labels"):
        rendered = graph.to_graphviz(show_edge_labels=False)

    assert isinstance(rendered, fake_graphviz)
    assert _rendered_edges(rendered) == [("A", "B", "3", "normal")]


def test_to_graphviz_can_force_count_edge_labels_on_collapsed_graph(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(collapse="family", edge_label="count")

    assert isinstance(rendered, fake_graphviz)
    assert _rendered_edges(rendered) == [
        ("TF", "g1|g2", "3", "normal"),
        ("g1|g2", "out", "2", "tee"),
    ]


def test_to_graphviz_rejects_unknown_edge_label():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="unsupported edge_label"):
        graph.to_graphviz(edge_label=cast(Any, "bad"))


def test_to_graphviz_can_remove_isolated_nodes(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=4)
    graph.add_edge("C", "D", sign=-1, count=1)

    rendered = graph.to_graphviz(
        min_frequency=0.5,
        remove_isolated_nodes=True,
    )

    assert isinstance(rendered, fake_graphviz)
    assert sorted(node for node, _ in rendered.nodes) == ["A", "B"]
    assert sorted((source, target) for source, target, _ in rendered.edges) == [
        ("A", "B"),
    ]


def test_to_graphviz_can_style_nodes_from_function_count(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("A", function_count=1)
    graph.add_node("B", function_count=2)
    graph.add_edge("A", "B", sign=1, count=4)

    rendered = graph.to_graphviz(node_style="count")

    assert isinstance(rendered, fake_graphviz)

    node_attrs = {node: attrs for node, attrs in rendered.nodes}

    assert node_attrs["A"]["fillcolor"] == "darkgoldenrod2"
    assert node_attrs["A"]["style"] == "rounded,filled,bold"
    assert node_attrs["B"]["fillcolor"] == "lightgoldenrod1"
    assert node_attrs["B"]["style"] == "rounded,filled"


def test_to_graphviz_node_style_warns_for_partially_missing_attributes(
    fake_graphviz,
):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("A", function_stability=1.0)
    graph.add_node("B")
    graph.add_edge("A", "B", sign=1, count=4)

    with pytest.warns(UserWarning, match="B"):
        rendered = graph.to_graphviz(node_style="stability")

    node_attrs = {node: attrs for node, attrs in rendered.nodes}

    assert node_attrs["A"]["fillcolor"] == "darkgoldenrod2"
    assert "fillcolor" not in node_attrs["B"]


def test_to_graphviz_node_style_warning_truncates_missing_nodes(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("styled", function_stability=1.0)

    for index in range(10):
        graph.add_node(f"missing_{index}")

    graph.add_edge("styled", "missing_0", sign=1, count=4)

    with pytest.warns(UserWarning) as warnings:
        graph.to_graphviz(node_style="stability")

    message = str(warnings[0].message)

    assert "missing_0" in message
    assert "missing_7" in message
    assert "missing_8" not in message
    assert "and 2 more" in message


def test_to_graphviz_node_style_warns_without_listing_when_all_missing(
    fake_graphviz,
):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=4)

    with pytest.warns(UserWarning) as warnings:
        rendered = graph.to_graphviz(node_style="stability")

    message = str(warnings[0].message)

    assert isinstance(rendered, fake_graphviz)
    assert "no node was styled" in message
    assert "A" not in message
    assert "B" not in message


def test_to_graphviz_node_style_can_be_inferred_from_complete_attributes(
    fake_graphviz,
):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("A", function_stability=1.0)
    graph.add_node("B", function_stability=0.5)
    graph.add_edge("A", "B", sign=1, count=4)

    rendered = graph.to_graphviz(node_style=True)

    node_attrs = {node: attrs for node, attrs in rendered.nodes}

    assert node_attrs["A"]["fillcolor"] == "darkgoldenrod2"
    assert node_attrs["B"]["fillcolor"] == "cornsilk"


def test_to_graphviz_rejects_unknown_node_style():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="unsupported node_style"):
        graph.to_graphviz(node_style=cast(Any, "bad"))


def test_to_pydot_supports_collapse_modes():
    pytest.importorskip("pydot")

    graph = _collapse_graph()

    exact = graph.to_pydot(rankdir="LR")
    family = graph.to_pydot(collapse="family")

    assert cast(Any, exact).get_rankdir() == "LR"
    assert sorted(
        (
            _pydot_get_string(edge, "get_source"),
            _pydot_get_string(edge, "get_destination"),
            _pydot_get_string(edge, "get_label"),
            _pydot_get_string(edge, "get_arrowhead"),
        )
        for edge in exact.get_edges()
    ) == [
        ("TF", "g1", "3", "normal"),
        ("TF", "g2", "3", "normal"),
        ("g1", "out", "2", "tee"),
        ("g2", "out", "2", "tee"),
    ]
    assert sorted(
        (
            _pydot_get_string(edge, "get_source"),
            _pydot_get_string(edge, "get_destination"),
            _pydot_get_string(edge, "get_label"),
            _pydot_get_string(edge, "get_arrowhead"),
        )
        for edge in family.get_edges()
    ) == [
        ("TF", "g1|g2", "3", "normal"),
        ("g1|g2", "out", "2", "tee"),
    ]


def test_show_supports_collapse_modes(monkeypatch):
    calls = {}

    class FakeDot:
        def create_svg(self):
            return b'<svg width="100pt" height="200pt"></svg>'

    def fake_to_pydot(
        self,
        collapse=None,
        bins=(0.0, 0.25, 0.5, 0.75, 1.0),
        protect_feedback_nodes=True,
        include_singleton_selfloops=True,
        min_frequency=0.0,
        remove_isolated_nodes=False,
        node_style=None,
        show_edge_labels="not passed",
        edge_label="count",
        edge_style=None,
        program="dot",
        **kwargs,
    ):
        calls["to_pydot"] = {
            "collapse": collapse,
            "bins": bins,
            "protect_feedback_nodes": protect_feedback_nodes,
            "include_singleton_selfloops": include_singleton_selfloops,
            "min_frequency": min_frequency,
            "remove_isolated_nodes": remove_isolated_nodes,
            "node_style": node_style,
            "edge_label": edge_label,
            "edge_style": edge_style,
            "program": program,
            "kwargs": kwargs,
        }
        return FakeDot()

    class FakeSVG:
        def __init__(self, svg):
            self.svg = svg

    def fake_display(svg):
        calls["display"] = svg

    display_module = ModuleType("IPython.display")
    setattr(display_module, "SVG", FakeSVG)
    setattr(display_module, "display", fake_display)
    monkeypatch.setitem(sys.modules, "IPython.display", display_module)
    monkeypatch.setattr(AggregatedInfluenceGraph, "to_pydot", fake_to_pydot)

    graph = _collapse_graph()
    graph.show(
        collapse="family",
        bins=(0.0, 0.5, 1.0),
        protect_feedback_nodes=False,
        include_singleton_selfloops=False,
        min_frequency=0.5,
        remove_isolated_nodes=True,
        node_style="stability",
        edge_label="frequency",
        edge_style=False,
        program="neato",
        width="640px",
        rankdir="LR",
    )

    assert calls["to_pydot"] == {
        "collapse": "family",
        "bins": (0.0, 0.5, 1.0),
        "protect_feedback_nodes": False,
        "include_singleton_selfloops": False,
        "min_frequency": 0.5,
        "remove_isolated_nodes": True,
        "node_style": "stability",
        "edge_label": "frequency",
        "edge_style": False,
        "program": "neato",
        "kwargs": {"rankdir": "LR"},
    }
    assert isinstance(calls["display"], FakeSVG)
    assert calls["display"].svg == '<svg width="640px"></svg>'
