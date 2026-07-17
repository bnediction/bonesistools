#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, cast

import pytest

from bonesistools.logic.influence_graph import (
    AggregatedInfluenceGraph,
    InfluenceGraph,
)


def _rendered_edges(graph):
    return sorted(
        (source, target, attrs["label"], attrs["arrowhead"])
        for source, target, attrs in graph.edges
    )


def _rendered_nodes(graph):
    return {node: attrs for node, attrs in graph.nodes}


def _pydot_get_string(obj, method):
    return cast(str, getattr(cast(Any, obj), method)()).strip('"')


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


def test_structural_families_can_ignore_edge_frequencies():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edges_from(
        [
            ("TF", "g1", {"sign": 1, "count": 3}),
            ("TF", "g2", {"sign": 1, "count": 3}),
            ("TF", "g3", {"sign": 1, "count": 1}),
        ]
    )

    frequency_aware = graph.structural_families(
        include_successors=False,
        bins=(0.0, 0.5, 1.0),
    )
    signed_structure_only = graph.structural_families(
        include_successors=False,
        bins=None,
    )
    collapsed = graph.family_collapsed_graph(
        include_successors=False,
        preserve_feedback=False,
        bins=None,
    )

    assert set(map(frozenset, frequency_aware.values())) == {frozenset({"g1", "g2"})}
    assert set(map(frozenset, signed_structure_only.values())) == {
        frozenset({"g1", "g2", "g3"})
    }
    assert sorted(collapsed.nodes()) == ["TF", "g1|g2|g3"]
    assert collapsed.nodes["g1|g2|g3"]["members"] == {"g1", "g2", "g3"}


def test_structural_families_validate_frequency_bins():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("TF", "g1", sign=1, count=3)

    with pytest.raises(ValueError, match="at least two"):
        graph.structural_families(bins=(0.0,))

    with pytest.raises(ValueError, match="strictly increasing"):
        graph.structural_families(bins=(0.0, 0.5, 0.5, 1.0))

    with pytest.raises(ValueError, match="cover"):
        graph.structural_families(bins=(0.25, 1.0))


def test_doc_example_family_and_plain_collapsed_graphs():
    graph = _collapse_graph()

    family_graph = graph.family_collapsed_graph(preserve_feedback=False)

    assert isinstance(family_graph, InfluenceGraph)
    assert not isinstance(family_graph, AggregatedInfluenceGraph)
    assert sorted(family_graph.nodes()) == ["TF", "g1|g2", "out"]
    assert family_graph.nodes["g1|g2"]["members"] == {"g1", "g2"}
    assert cast(Any, family_graph["TF"]["g1|g2"])[0]["frequency"] == 0.75


def test_family_collapsed_graph_uses_pipe_separator_by_default():
    graph = _collapse_graph()

    collapsed = graph.family_collapsed_graph(preserve_feedback=False)

    assert "g1|g2" in collapsed
    assert "g1/g2" not in collapsed


def test_family_collapsed_graph_accepts_custom_separator():
    graph = _collapse_graph()

    collapsed = graph.family_collapsed_graph(
        preserve_feedback=False,
        sep="/",
    )

    assert "g1/g2" in collapsed
    assert "g1|g2" not in collapsed


def test_doc_example_feedback_induced_graph():
    graph = _feedback_graph()

    feedback = graph.feedback_induced_graph()

    assert list(feedback.nodes()) == ["A", "B"]
    assert sorted(feedback.nodes()) == ["A", "B"]
    assert feedback.total == 3


def test_doc_example_collapsed_graph():
    graph = _feedback_graph()

    collapsed = graph.collapsed_graph()

    assert isinstance(collapsed, InfluenceGraph)
    assert not isinstance(collapsed, AggregatedInfluenceGraph)
    assert sorted(collapsed.nodes()) == ["A", "B"]


def test_structural_families_can_preserve_feedback_nodes():
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
        preserve_feedback=True,
    )

    assert set(map(frozenset, families.values())) == {frozenset({"g1", "g2"})}
    assert all("A" not in family and "B" not in family for family in families.values())


def test_family_collapsed_graph_ignores_singleton_families():
    graph = AggregatedInfluenceGraph(total=4)
    graph.add_edge("TF", "g1", sign=1, count=3)

    collapsed = graph.family_collapsed_graph(preserve_feedback=False)

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
        preserve_feedback=False,
    )

    assert sorted(collapsed.nodes()) == ["TF", "g1", "g2"]
    assert sorted(collapsed.edges()) == [
        ("TF", "g1"),
        ("TF", "g2"),
        ("g1", "g2"),
        ("g2", "g1"),
    ]


def test_copy_and_collapse_preserve_expected_types():
    graph = _collapse_graph()

    copied = graph.copy()
    collapsed = graph.family_collapsed_graph(preserve_feedback=False)

    assert isinstance(copied, AggregatedInfluenceGraph)
    assert copied.total == 4
    assert isinstance(collapsed, InfluenceGraph)
    assert sorted(collapsed.nodes()) == ["TF", "g1|g2", "out"]
    assert cast(Any, collapsed["TF"]["g1|g2"])[0]["frequency"] == 0.75
    assert cast(Any, collapsed["g1|g2"]["out"])[0]["frequency"] == 0.5


def test_to_graphviz_supports_collapse_modes(fake_graphviz):
    graph = _collapse_graph()

    exact = cast(Any, graph.to_graphviz(graph_attr={"rankdir": "LR"}))
    family = cast(Any, graph.to_graphviz(collapse="family"))

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
    assert _rendered_nodes(family)["g1|g2"]["label"] == "g1|g2"
    assert _rendered_nodes(family)["g1|g2"]["margin"] == "0.08,0.04"
    assert _rendered_nodes(family)["g1|g2"]["shape"] == "box"
    assert _rendered_nodes(family)["g1|g2"]["style"] == "rounded"

    feedback_graph = AggregatedInfluenceGraph(total=4)
    feedback_graph.add_edges_from(
        [
            ("A", "B", {"sign": 1, "count": 3}),
            ("B", "A", {"sign": 1, "count": 3}),
            ("B", "C", {"sign": -1, "count": 1}),
        ]
    )

    feedback = cast(Any, feedback_graph.to_graphviz(collapse="feedback"))
    both = cast(Any, feedback_graph.to_graphviz(collapse="both"))

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


def test_to_graphviz_can_force_count_edge_labels_on_collapsed_graph(fake_graphviz):
    graph = _collapse_graph()

    rendered = cast(
        Any,
        graph.to_graphviz(collapse="family", edge_label="count"),
    )

    assert isinstance(rendered, fake_graphviz)
    assert _rendered_edges(rendered) == [
        ("TF", "g1|g2", "3", "normal"),
        ("g1|g2", "out", "2", "tee"),
    ]


def test_to_graphviz_styles_family_from_lowest_member_stability(fake_graphviz):
    graph = _collapse_graph()
    graph.nodes["TF"]["function_stability"] = 1.0
    graph.nodes["g1"]["function_stability"] = 1.0
    graph.nodes["g2"]["function_stability"] = 0.75
    graph.nodes["out"]["function_stability"] = 1.0

    rendered = cast(
        Any,
        graph.to_graphviz(collapse="family", node_style="stability"),
    )
    family_attrs = _rendered_nodes(rendered)["g1|g2"]

    assert family_attrs["function_stability"] == "0.75"
    assert family_attrs["fillcolor"] == "lightgoldenrod1"
    assert family_attrs["shape"] == "box"
    assert set(family_attrs["style"].split(",")) == {"rounded", "filled"}


def test_to_graphviz_styles_family_as_orange_when_all_members_are_orange(
    fake_graphviz,
):
    graph = _collapse_graph()
    graph.nodes["TF"]["function_stability"] = 1.0
    graph.nodes["g1"]["function_stability"] = 1.0
    graph.nodes["g2"]["function_stability"] = 1.0
    graph.nodes["out"]["function_stability"] = 1.0

    rendered = cast(
        Any,
        graph.to_graphviz(collapse="family", node_style="stability"),
    )
    family_attrs = _rendered_nodes(rendered)["g1|g2"]

    assert family_attrs["function_stability"] == "1.0"
    assert family_attrs["fillcolor"] == "darkgoldenrod2"
    assert family_attrs["shape"] == "box"
    assert set(family_attrs["style"].split(",")) == {"rounded", "filled", "bold"}


def test_to_graphviz_styles_family_from_highest_member_function_count(
    fake_graphviz,
):
    graph = _collapse_graph()
    graph.nodes["TF"]["function_count"] = 1
    graph.nodes["g1"]["function_count"] = 1
    graph.nodes["g2"]["function_count"] = 2
    graph.nodes["out"]["function_count"] = 1

    rendered = cast(
        Any,
        graph.to_graphviz(collapse="family", node_style="count"),
    )
    family_attrs = _rendered_nodes(rendered)["g1|g2"]

    assert family_attrs["function_count"] == "2"
    assert family_attrs["fillcolor"] == "lightgoldenrod1"
    assert family_attrs["shape"] == "box"
    assert set(family_attrs["style"].split(",")) == {"rounded", "filled"}


def test_to_pydot_supports_collapse_modes():
    pytest.importorskip("pydot")

    graph = _collapse_graph()

    exact = graph.to_pydot(graph_attr={"rankdir": "LR"})
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
    family_nodes = {
        _pydot_get_string(node, "get_name"): node for node in family.get_nodes()
    }
    family_node = cast(Any, family_nodes["g1|g2"])

    assert _pydot_get_string(family_node, "get_label") == "g1|g2"
    assert _pydot_get_string(family_node, "get_margin") == "0.08,0.04"
    assert _pydot_get_string(family_node, "get_shape") == "box"
    assert _pydot_get_string(family_node, "get_style") == "rounded"


def test_to_graphviz_wraps_long_family_labels(fake_graphviz):
    graph = AggregatedInfluenceGraph(total=4)
    family = (
        "Apex1",
        "Cbx3",
        "Exosc8",
        "Fasn",
        "Gnl3",
        "Hspa4",
        "Hspd1",
        "Hspe1",
        "Nap1l1",
        "Ncl",
        "Nme1",
        "Nop56",
        "Odc1",
        "Ppat",
        "Ppid",
        "Rcc1",
        "Srm",
        "Tfrc",
        "Thrap3",
    )
    for gene in family:
        graph.add_edge("TF", gene, sign=1, count=3)

    rendered = cast(Any, graph.to_graphviz(collapse="family"))
    family_name = "|".join(sorted(family))
    label = _rendered_nodes(rendered)[family_name]["label"]

    assert "\n" in label
    assert label.replace("\n", "|") == family_name
    assert all(len(line) <= 40 for line in label.splitlines())
    assert _rendered_nodes(rendered)[family_name]["margin"] == "0.08,0.04"
    assert _rendered_nodes(rendered)[family_name]["shape"] == "box"
    assert _rendered_nodes(rendered)[family_name]["style"] == "rounded"


def test_to_graphviz_family_attr_can_be_disabled(fake_graphviz):
    graph = _collapse_graph()

    rendered = cast(
        Any,
        graph.to_graphviz(collapse="family", family_attr=False),
    )
    family_attrs = _rendered_nodes(rendered)["g1|g2"]

    assert "height" not in family_attrs
    assert "label" not in family_attrs
    assert "margin" not in family_attrs
    assert "shape" not in family_attrs
    assert "style" not in family_attrs


def test_to_graphviz_family_attr_updates_default_family_attributes(fake_graphviz):
    graph = _collapse_graph()

    rendered = cast(
        Any,
        graph.to_graphviz(
            collapse="family",
            family_attr={
                "height": "0.1",
                "margin": "0.02",
                "shape": "oval",
                "style": "dashed",
            },
        ),
    )
    family_attrs = _rendered_nodes(rendered)["g1|g2"]

    assert family_attrs["label"] == "g1|g2"
    assert family_attrs["height"] == "0.1"
    assert family_attrs["margin"] == "0.02"
    assert family_attrs["shape"] == "oval"
    assert family_attrs["style"] == "dashed"


def test_show_supports_collapse_modes(monkeypatch):
    calls = {}

    class FakeDot:
        def create_svg(self):
            return b'<svg width="100pt" height="200pt"></svg>'

    def fake_to_pydot(
        self,
        collapse=None,
        bins=(0.0, 0.25, 0.5, 0.75, 1.0),
        preserve_feedback=True,
        include_selfloops=True,
        min_frequency=0.0,
        drop_isolates=False,
        node_style=None,
        family_attr=True,
        edge_label="count",
        edge_style=None,
        program="dot",
        graph_attr=None,
        node_attr=None,
        edge_attr=None,
        **kwargs,
    ):
        calls["to_pydot"] = {
            "collapse": collapse,
            "bins": bins,
            "preserve_feedback": preserve_feedback,
            "include_selfloops": include_selfloops,
            "min_frequency": min_frequency,
            "drop_isolates": drop_isolates,
            "node_style": node_style,
            "family_attr": family_attr,
            "edge_label": edge_label,
            "edge_style": edge_style,
            "program": program,
            "graph_attr": graph_attr,
            "node_attr": node_attr,
            "edge_attr": edge_attr,
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
        preserve_feedback=False,
        include_selfloops=False,
        min_frequency=0.5,
        drop_isolates=True,
        node_style="stability",
        family_attr={"shape": "oval"},
        edge_label="frequency",
        edge_style=None,
        program="neato",
        graph_attr={"rankdir": "LR"},
        node_attr={"fontsize": 20},
        edge_attr={"fontcolor": "gray"},
        width="640px",
    )

    assert calls["to_pydot"] == {
        "collapse": "family",
        "bins": (0.0, 0.5, 1.0),
        "preserve_feedback": False,
        "include_selfloops": False,
        "min_frequency": 0.5,
        "drop_isolates": True,
        "node_style": "stability",
        "family_attr": {"shape": "oval"},
        "edge_label": "frequency",
        "edge_style": None,
        "program": "neato",
        "graph_attr": {"rankdir": "LR"},
        "node_attr": {"fontsize": 20},
        "edge_attr": {"fontcolor": "gray"},
        "kwargs": {},
    }
    assert isinstance(calls["display"], FakeSVG)
    assert calls["display"].svg == '<svg width="640px"></svg>'
