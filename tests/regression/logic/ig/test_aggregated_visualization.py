#!/usr/bin/env python

from typing import Any, cast

import pytest

from bonesistools.logic.influence_graph import AggregatedInfluenceGraph


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


def test_to_graphviz_can_filter_and_disable_frequency_style(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(
        collapse="family",
        min_frequency=0.6,
        edge_label=None,
        edge_style=None,
    )

    assert isinstance(rendered, fake_graphviz)
    assert sorted((source, target) for source, target, _ in rendered.edges) == [
        ("TF", "g1|g2"),
    ]
    assert all("label" not in attrs for _, _, attrs in rendered.edges)
    assert all("style" not in attrs for _, _, attrs in rendered.edges)


def test_to_graphviz_accepts_default_graph_node_and_edge_attributes(fake_graphviz):
    graph = _single_edge_graph()

    rendered = graph.to_graphviz(
        graph_attr={"rankdir": "LR", "bgcolor": None},
        node_attr={"fontsize": 20},
        edge_attr={"fontcolor": "gray"},
    )

    assert isinstance(rendered, fake_graphviz)
    assert rendered.graph_attr == {"rankdir": "LR"}
    assert rendered.node_attr == {"fontsize": "20"}
    assert rendered.edge_attr == {"fontcolor": "gray"}
    assert all("fontsize" not in attrs for _, attrs in rendered.nodes)
    assert all("fontcolor" not in attrs for _, _, attrs in rendered.edges)


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


def test_to_graphviz_can_label_edges_from_custom_attribute(fake_graphviz):
    graph = _single_edge_graph()
    cast(Any, graph["A"]["B"])[0]["support"] = "manual"

    rendered = graph.to_graphviz(edge_label="support")

    assert isinstance(rendered, fake_graphviz)
    assert _rendered_edges(rendered) == [("A", "B", "manual", "normal")]


def test_to_graphviz_can_disable_edge_labels_with_edge_label_none(fake_graphviz):
    graph = _collapse_graph()

    rendered = graph.to_graphviz(edge_label=None)

    assert isinstance(rendered, fake_graphviz)
    assert all("label" not in attrs for _, _, attrs in rendered.edges)


def test_to_graphviz_rejects_missing_edge_label_attribute():
    graph = _single_edge_graph()

    with pytest.raises(KeyError, match="missing edge attribute"):
        graph.to_graphviz(edge_label=cast(Any, "bad"))


def test_to_graphviz_rejects_non_string_edge_label():
    graph = _single_edge_graph()

    with pytest.raises(TypeError, match="edge_label"):
        graph.to_graphviz(edge_label=cast(Any, 1))


def test_to_graphviz_can_drop_isolates(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("A", "B", sign=1, count=4)
    graph.add_edge("C", "D", sign=-1, count=1)

    rendered = graph.to_graphviz(
        min_frequency=0.5,
        drop_isolates=True,
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


def test_to_graphviz_can_style_nodes_from_function_stability(
    fake_graphviz,
):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("A", function_stability=1.0)
    graph.add_node("B", function_stability=0.5)
    graph.add_edge("A", "B", sign=1, count=4)

    rendered = graph.to_graphviz(node_style="stability")

    node_attrs = {node: attrs for node, attrs in rendered.nodes}

    assert node_attrs["A"]["fillcolor"] == "darkgoldenrod2"
    assert node_attrs["B"]["fillcolor"] == "cornsilk"


def test_to_graphviz_callable_node_style_uses_named_attributes(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_node("A", function_stability=1.0, function_count=2)
    graph.add_node("B", function_stability=0.25, function_count=10)
    graph.add_edge("A", "B", sign=1, count=4)

    def style(function_stability, function_count):
        return {
            "fillcolor": "white" if function_stability >= 0.5 else "lightgray",
            "penwidth": str(function_count),
        }

    rendered = graph.to_graphviz(node_style=style)

    node_attrs = {node: attrs for node, attrs in rendered.nodes}

    assert isinstance(rendered, fake_graphviz)
    assert node_attrs["A"]["fillcolor"] == "white"
    assert node_attrs["A"]["penwidth"] == "2"
    assert node_attrs["B"]["fillcolor"] == "lightgray"
    assert node_attrs["B"]["penwidth"] == "10"


def test_to_graphviz_callable_node_style_rejects_kwargs():
    graph = _single_edge_graph()

    with pytest.raises(TypeError, match="explicitly named parameters"):
        graph.to_graphviz(node_style=lambda **data: {})


def test_to_graphviz_callable_node_style_rejects_missing_attributes():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="function_stability"):
        graph.to_graphviz(node_style=lambda function_stability: {})


def test_to_graphviz_callable_node_style_rejects_varargs():
    graph = _single_edge_graph()

    def style(*data):
        return {}

    with pytest.raises(TypeError, match="named parameters"):
        graph.to_graphviz(node_style=style)


def test_to_graphviz_rejects_unknown_node_style():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="unsupported node_style"):
        graph.to_graphviz(node_style=cast(Any, "bad"))


def test_to_graphviz_rejects_invalid_edge_style():
    graph = _single_edge_graph()

    with pytest.raises(TypeError, match="edge_style"):
        graph.to_graphviz(edge_style=cast(Any, False))


def test_to_graphviz_rejects_unknown_edge_style():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="unsupported edge_style"):
        graph.to_graphviz(edge_style=cast(Any, "bad"))


def test_to_graphviz_accepts_frequency_edge_style(fake_graphviz):
    graph = _single_edge_graph(total=4)

    rendered = graph.to_graphviz(edge_style="frequency")

    assert isinstance(rendered, fake_graphviz)
    assert rendered.edges[0][2]["style"] == "bold"
    assert rendered.edges[0][2]["penwidth"] == "2"


def test_to_graphviz_callable_edge_style_uses_named_attributes(fake_graphviz):
    graph = _single_edge_graph(total=4)

    rendered = graph.to_graphviz(
        edge_style=lambda frequency, count, sign: {
            "tooltip": f"{frequency:.2f}:{count}:{sign}"
        }
    )

    assert isinstance(rendered, fake_graphviz)
    assert rendered.edges[0][2]["tooltip"] == "0.75:3:1"


def test_to_graphviz_callable_edge_style_rejects_kwargs():
    graph = _single_edge_graph(total=4)

    with pytest.raises(TypeError, match="explicitly named parameters"):
        graph.to_graphviz(edge_style=lambda **data: {})


def test_to_graphviz_callable_edge_style_rejects_missing_attributes():
    graph = _single_edge_graph()

    with pytest.raises(ValueError, match="missing"):
        graph.to_graphviz(edge_style=lambda missing: {})


def test_to_graphviz_callable_edge_style_rejects_varargs():
    graph = _single_edge_graph()

    def style(*data):
        return {}

    with pytest.raises(TypeError, match="named parameters"):
        graph.to_graphviz(edge_style=style)


def test_to_pydot_accepts_default_graph_node_and_edge_attributes():
    pytest.importorskip("pydot")

    graph = _single_edge_graph()

    dot = graph.to_pydot(
        graph_attr={"rankdir": "LR", "bgcolor": None},
        node_attr={"fontsize": 20},
        edge_attr={"fontcolor": "gray"},
    )

    assert cast(Any, dot).get_rankdir() == "LR"
    assert cast(Any, dot).get_node_defaults() == [{"fontsize": "20"}]
    assert cast(Any, dot).get_edge_defaults() == [{"fontcolor": "gray"}]


def test_to_pydot_applies_node_and_edge_attributes_to_rendered_items():
    pytest.importorskip("pydot")

    graph = _single_edge_graph()
    graph.nodes["A"]["function_stability"] = 1.0
    graph.nodes["B"]["function_stability"] = 0.5

    dot = graph.to_pydot(
        node_attr={"fontsize": 100},
        node_style="stability",
        edge_attr={"fontsize": 100},
        edge_style="frequency",
    )

    node_attrs = {
        _pydot_get_string(node, "get_name"): cast(Any, node).get_attributes()
        for node in cast(Any, dot).get_nodes()
        if _pydot_get_string(node, "get_name") in {"A", "B"}
    }
    edge_attrs = cast(Any, dot).get_edges()[0].get_attributes()

    assert node_attrs["A"]["fontsize"] == "100"
    assert node_attrs["A"]["fillcolor"] == "darkgoldenrod2"
    assert node_attrs["B"]["fontsize"] == "100"
    assert node_attrs["B"]["fillcolor"] == "cornsilk"
    assert edge_attrs["fontsize"] == "100"
    assert edge_attrs["style"] == "bold"
