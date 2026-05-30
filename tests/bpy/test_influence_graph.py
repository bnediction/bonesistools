#!/usr/bin/env python

import builtins
import sys
from types import ModuleType

import networkx as nx
import pytest

import bonesistools as bt


def _edge_signs(graph):
    return sorted(
        (source, target, sign) for source, target, sign in graph.edges(data="sign")
    )


def _cycle_sets(circuits):
    return {(frozenset(cycle), sign) for cycle, sign in circuits}


def test_influence_graph_add_edge_normalizes_signs_and_rejects_duplicates():
    ig = bt.bpy.ig.InfluenceGraph()

    ig.add_edge("A", "B", sign=1)
    ig.add_edge("A", "B", sign="-")

    assert str(ig) == "InfluenceGraph(nodes=2, edges=2)"
    assert repr(ig) == str(ig)
    assert _edge_signs(ig) == [("A", "B", -1), ("A", "B", 1)]

    with pytest.raises(ValueError, match="duplicated edge sign"):
        ig.add_edge("A", "B", sign=1)

    with pytest.raises(ValueError, match="unsupported edge sign"):
        ig.add_edge("B", "C", sign=0)

    assert _edge_signs(ig) == [("A", "B", -1), ("A", "B", 1)]

    key_style = bt.bpy.ig.InfluenceGraph()
    key_style.add_edge("A", "B", 1)

    assert _edge_signs(key_style) == [("A", "B", 1)]

    with pytest.raises(TypeError, match="missing required argument"):
        key_style.add_edge("B", "C")


def test_influence_graph_constructor_validates_existing_graph():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "B", sign="+")
    graph.add_edge("B", "C", sign="negative")

    ig = bt.bpy.ig.InfluenceGraph(graph)

    assert _edge_signs(ig) == [("A", "B", 1), ("B", "C", -1)]

    missing_sign = nx.MultiDiGraph()
    missing_sign.add_edge("A", "B")
    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        bt.bpy.ig.InfluenceGraph(missing_sign)

    duplicate_sign = nx.MultiDiGraph()
    duplicate_sign.add_edge("A", "B", sign=1)
    duplicate_sign.add_edge("A", "B", sign="+")
    with pytest.raises(ValueError, match="duplicated edge sign"):
        bt.bpy.ig.InfluenceGraph(duplicate_sign)


def test_influence_graph_add_edges_from_is_transactional():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": -1}),
        ]
    )

    assert _edge_signs(ig) == [("A", "B", 1), ("B", "C", -1)]

    ig.add_edges_from(
        [
            ("C", "D"),
            ("D", "E"),
        ],
        sign=1,
    )

    assert _edge_signs(ig) == [
        ("A", "B", 1),
        ("B", "C", -1),
        ("C", "D", 1),
        ("D", "E", 1),
    ]

    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        ig.add_edges_from([("E", "F")])

    assert "F" not in ig
    assert _edge_signs(ig) == [
        ("A", "B", 1),
        ("B", "C", -1),
        ("C", "D", 1),
        ("D", "E", 1),
    ]


def test_update_is_transactional_and_views_are_disabled():
    g1 = bt.bpy.ig.InfluenceGraph()
    g1.add_edge("A", "B", sign=1)

    g2 = bt.bpy.ig.InfluenceGraph()
    g2.add_edge("B", "C", sign=-1)
    g2.add_edge("C", "A", sign=1)

    g1.update(g2)

    assert _edge_signs(g1) == [("A", "B", 1), ("B", "C", -1), ("C", "A", 1)]

    adjacency = {
        "A": ["B", "C"],
        "B": ["C"],
    }
    edges = [
        (source, target, {"sign": 1})
        for source, targets in adjacency.items()
        for target in targets
    ]

    ig = bt.bpy.ig.InfluenceGraph()
    ig.update(edges=edges, nodes=adjacency)

    assert _edge_signs(ig) == [("A", "B", 1), ("A", "C", 1), ("B", "C", 1)]

    previous_edges = _edge_signs(ig)
    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        ig.update(edges=[("B", "D")])

    assert "D" not in ig
    assert _edge_signs(ig) == previous_edges

    copied = ig.copy()
    directed = ig.to_directed()

    assert isinstance(copied, bt.bpy.ig.InfluenceGraph)
    assert isinstance(directed, bt.bpy.ig.InfluenceGraph)
    assert copied is not ig
    assert directed is not ig

    with pytest.raises(NotImplementedError, match="view copies"):
        ig.copy(as_view=True)

    with pytest.raises(NotImplementedError, match="directed views"):
        ig.to_directed(as_view=True)

    with pytest.raises(NotImplementedError, match="conversion to undirected"):
        ig.to_undirected()

    with pytest.raises(NotImplementedError, match="add_weighted_edges_from"):
        ig.add_weighted_edges_from([("A", "B", 1.0)])


def test_scc_feedback_regulators_and_targets():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": 1}),
            ("C", "A", {"sign": 1}),
            ("C", "D", {"sign": 1}),
            ("D", "E", {"sign": 1}),
            ("E", "D", {"sign": 1}),
            ("F", "G", {"sign": 1}),
            ("H", "H", {"sign": 1}),
        ]
    )

    assert {frozenset(scc) for scc in ig.strongly_connected_components()} == {
        frozenset({"A", "B", "C"}),
        frozenset({"D", "E"}),
        frozenset({"H"}),
    }
    assert {
        frozenset(scc)
        for scc in ig.strongly_connected_components(include_singleton_selfloops=False)
    } == {frozenset({"A", "B", "C"}), frozenset({"D", "E"})}
    assert ig.feedback_nodes() == {"A", "B", "C", "D", "E", "H"}

    regulator_graph = bt.bpy.ig.InfluenceGraph()
    regulator_graph.add_edges_from(
        [
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": 1}),
            ("D", "B", {"sign": -1}),
        ]
    )
    regulator_graph.add_node("E")

    assert set(regulator_graph.regulators()) == {"A", "D"}
    assert set(regulator_graph.targets()) == {"C"}
    assert "E" not in regulator_graph.regulators()
    assert "E" not in regulator_graph.targets()


def test_structural_families_and_family_collapsed_graph_from_docstring():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("TF1", "g1", {"sign": 1}),
            ("TF1", "g2", {"sign": 1}),
            ("TF1", "g3", {"sign": 1}),
            ("TF1", "g4", {"sign": 1}),
            ("TF2", "g3", {"sign": -1}),
            ("TF2", "g4", {"sign": -1}),
        ]
    )

    families = ig.structural_families()

    assert set(map(frozenset, families.values())) == {
        frozenset({"g1", "g2"}),
        frozenset({"g3", "g4"}),
    }

    collapsed = ig.family_collapsed_graph()

    assert isinstance(collapsed, bt.bpy.ig.InfluenceGraph)
    assert sorted(collapsed.nodes()) == ["TF1", "TF2", "g1|g2", "g3|g4"]
    assert collapsed.nodes["g1|g2"]["members"] == {"g1", "g2"}
    assert collapsed.nodes["g3|g4"]["members"] == {"g3", "g4"}
    assert _edge_signs(collapsed) == [
        ("TF1", "g1|g2", 1),
        ("TF1", "g3|g4", 1),
        ("TF2", "g3|g4", -1),
    ]

    with pytest.warns(DeprecationWarning, match="family_collapsed_graph"):
        deprecated = ig.family_compressed_graph()

    assert sorted(deprecated.nodes()) == sorted(collapsed.nodes())


def test_structural_families_can_ignore_successors_and_feedback_nodes():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("TF", "g1", {"sign": 1}),
            ("TF", "g2", {"sign": 1}),
            ("g1", "out1", {"sign": 1}),
            ("g2", "out2", {"sign": -1}),
            ("A", "B", {"sign": 1}),
            ("B", "A", {"sign": 1}),
        ]
    )

    families = ig.structural_families(
        include_successors=False,
        exclude_feedback_nodes=True,
    )

    assert set(map(frozenset, families.values())) == {frozenset({"g1", "g2"})}


def test_edge_sign_autoregulations_and_path_sign_from_docstrings():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edge("A", "B", sign=1)
    ig.add_edge("A", "B", sign=-1)
    ig.add_edge("B", "C", sign=-1)
    ig.add_edge("C", "D", sign=1)
    ig.add_edge("D", "E", sign=1)
    ig.add_edge("D", "E", sign=-1)
    ig.add_edge("F", "F", sign=1)
    ig.add_edge("G", "G", sign=-1)

    assert ig.edge_sign("A", "B") == 0
    assert ig.edge_sign("C", "D") == 1
    assert ig.path_sign(["A", "B"]) == 0
    assert ig.path_sign(["B", "C", "D"]) == -1
    assert ig.path_sign(["D", "E"]) == 0
    assert sorted(ig.autoregulations()) == [("F", 1), ("G", -1)]
    assert ig.autoregulations(sign=1) == [("F", 1)]

    with pytest.raises(KeyError, match="no edge found"):
        ig.edge_sign("X", "Y")

    with pytest.raises(ValueError, match="at least two nodes"):
        ig.path_sign(["A"])


def test_signed_path_string_formats_signed_edges():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edge("A", "B", sign=1)
    ig.add_edge("B", "C", sign=-1)
    ig.add_edge("C", "D", sign=1)
    ig.add_edge("E", "F", sign=1)
    ig.add_edge("E", "F", sign=-1)

    assert ig.signed_path_string("A", "B", "C", "D") == "A -> B -| C -> D"
    assert ig.signed_path_string("B", "C", "D") == "B -| C -> D"
    assert ig.signed_path_string(["B", "C", "D"]) == "B -| C -> D"
    assert ig.signed_path_string(("C", "D")) == "C -> D"
    assert ig.signed_path_string("E", "F") == "E -- F"

    with pytest.raises(ValueError, match="at least two nodes"):
        ig.signed_path_string("A")

    with pytest.raises(KeyError, match="no edge found"):
        ig.signed_path_string("A", "C")


def test_marker_paths_from_docstring():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("marker1", "A", {"sign": 1}),
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": 1}),
            ("C", "A", {"sign": 1}),
            ("C", "D", {"sign": 1}),
            ("D", "E", {"sign": 1}),
            ("E", "D", {"sign": 1}),
            ("D", "marker2", {"sign": 1}),
            ("E", "marker2", {"sign": 1}),
        ]
    )

    observed = ig.marker_paths(["marker2"], direction="downstream")
    expected = [
        {
            "marker": "marker2",
            "scc": frozenset({"A", "B", "C"}),
            "direction": "downstream",
            "paths": [
                {
                    "path": ["C", "D", "marker2"],
                    "sign": 1,
                },
            ],
        },
        {
            "marker": "marker2",
            "scc": frozenset({"D", "E"}),
            "direction": "downstream",
            "paths": [
                {
                    "path": ["D", "marker2"],
                    "sign": 1,
                },
                {
                    "path": ["E", "marker2"],
                    "sign": 1,
                },
            ],
        },
    ]

    def normalize(records):
        return sorted(
            (
                {
                    **record,
                    "paths": sorted(record["paths"], key=lambda path: path["path"]),
                }
                for record in records
            ),
            key=lambda record: sorted(record["scc"]),
        )

    assert normalize(observed) == normalize(expected)

    assert ig.marker_paths(["missing"]) == []

    with pytest.raises(ValueError, match="unsupported direction"):
        ig.marker_paths(["marker2"], direction="sideways")


def test_marker_paths_supports_upstream_paths_and_sign_filter():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("marker", "A", {"sign": -1}),
            ("A", "B", {"sign": 1}),
            ("B", "A", {"sign": 1}),
            ("C", "D", {"sign": 1}),
            ("D", "C", {"sign": 1}),
        ]
    )

    observed = ig.marker_paths(
        ["marker"],
        direction="upstream",
        sccs=[["A", "B"], ["C", "D"]],
        sign=-1,
    )

    assert observed == [
        {
            "marker": "marker",
            "scc": frozenset({"A", "B"}),
            "direction": "upstream",
            "paths": [{"path": ["marker", "A"], "sign": -1}],
        }
    ]
    assert (
        ig.marker_paths(
            ["marker"],
            direction="downstream",
            sccs=[["A", "B"]],
        )
        == []
    )


def test_circuits_positive_and_negative_circuits_from_docstring():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": -1}),
            ("C", "A", {"sign": 1}),
            ("D", "E", {"sign": -1}),
            ("E", "D", {"sign": -1}),
        ]
    )

    assert _cycle_sets(ig.circuits()) == {
        (frozenset({"A", "B", "C"}), -1),
        (frozenset({"D", "E"}), 1),
    }
    assert {frozenset(cycle) for cycle in ig.positive_circuits()} == {
        frozenset({"D", "E"})
    }
    assert {frozenset(cycle) for cycle in ig.negative_circuits()} == {
        frozenset({"A", "B", "C"})
    }
    assert _cycle_sets(ig.circuits(sign=1)) == {(frozenset({"D", "E"}), 1)}
    assert _cycle_sets(ig.circuits(sign=-1)) == {(frozenset({"A", "B", "C"}), -1)}

    ig.add_edge("A", "B", sign=-1)
    assert _cycle_sets(ig.circuits()) == {
        (frozenset({"A", "B", "C"}), -1),
        (frozenset({"A", "B", "C"}), 1),
        (frozenset({"D", "E"}), 1),
    }


def test_feedback_induced_graph_and_collapsed_graph_from_docstring():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edges_from(
        [
            ("A", "B", {"sign": 1}),
            ("B", "C", {"sign": 1}),
            ("C", "A", {"sign": 1}),
            ("C", "D", {"sign": 1}),
            ("D", "E", {"sign": 1}),
            ("E", "D", {"sign": 1}),
            ("C", "X", {"sign": 1}),
            ("X", "Y", {"sign": 1}),
            ("F", "G", {"sign": 1}),
            ("H", "H", {"sign": 1}),
        ]
    )

    feedback = ig.feedback_induced_graph()

    assert isinstance(feedback, bt.bpy.ig.InfluenceGraph)
    assert sorted(feedback.nodes()) == ["A", "B", "C", "D", "E", "H"]
    assert sorted(feedback.edges()) == [
        ("A", "B"),
        ("B", "C"),
        ("C", "A"),
        ("C", "D"),
        ("D", "E"),
        ("E", "D"),
        ("H", "H"),
    ]

    collapsed = ig.collapsed_graph()

    assert isinstance(collapsed, bt.bpy.ig.InfluenceGraph)
    assert set(collapsed.nodes()) <= set(feedback.nodes())

    with pytest.warns(DeprecationWarning, match="collapsed_graph"):
        deprecated = ig.compressed_graph()

    assert sorted(deprecated.nodes()) == sorted(collapsed.nodes())


def test_to_graphviz_applies_signed_edge_styles_and_custom_options(fake_graphviz):
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edge("A", "B", sign=1)
    ig.add_edge("B", "C", sign=-1)

    graph = ig.to_graphviz(
        program="neato",
        edge_style=lambda data: {"label": data["sign"]},
        rankdir="LR",
    )

    edge_styles = {(source, target): attrs for source, target, attrs in graph.edges}

    assert isinstance(graph, fake_graphviz)
    assert graph.engine == "neato"
    assert graph.graph_attr["rankdir"] == "LR"
    assert edge_styles[("A", "B")] == {
        "sign": "1",
        "color": "green4",
        "arrowhead": "normal",
        "penwidth": "2",
        "label": "1",
    }
    assert edge_styles[("B", "C")] == {
        "sign": "-1",
        "color": "red2",
        "arrowhead": "tee",
        "penwidth": "2",
        "label": "-1",
    }


def test_to_pydot_applies_signed_edge_styles_and_custom_options():
    ig = bt.bpy.ig.InfluenceGraph()
    ig.add_edge("A", "B", sign=1)
    ig.add_edge("B", "C", sign=-1)

    dot = ig.to_pydot(
        program="dot",
        edge_style=lambda data: {"label": str(data["sign"])},
        rankdir="LR",
    )

    edge_styles = {
        (edge.get_source().strip('"'), edge.get_destination().strip('"')): {
            "color": edge.get_color(),
            "arrowhead": edge.get_arrowhead(),
            "label": edge.get_label(),
        }
        for edge in dot.get_edges()
    }

    assert dot.prog == "dot"
    assert dot.get_rankdir() == "LR"
    assert edge_styles[("A", "B")] == {
        "color": "green4",
        "arrowhead": "normal",
        "label": "1",
    }
    assert edge_styles[("B", "C")] == {
        "color": "red2",
        "arrowhead": "tee",
        "label": "-1",
    }


def test_show_uses_ipython_display(monkeypatch):
    ig = bt.bpy.ig.InfluenceGraph()
    displayed = []

    class FakeDot:
        def create_svg(self):
            return b"<svg></svg>"

    display_module = ModuleType("IPython.display")
    display_module.SVG = lambda svg: ("SVG", svg)
    display_module.display = displayed.append

    monkeypatch.setitem(sys.modules, "IPython", ModuleType("IPython"))
    monkeypatch.setitem(sys.modules, "IPython.display", display_module)
    monkeypatch.setattr(
        bt.bpy.ig.InfluenceGraph,
        "to_pydot",
        lambda self, **kwargs: FakeDot(),
    )

    ig.show()

    assert displayed == [("SVG", "<svg></svg>")]


def test_show_requires_ipython(monkeypatch):
    ig = bt.bpy.ig.InfluenceGraph()
    original_import = builtins.__import__

    def import_without_ipython(name, *args, **kwargs):
        if name == "IPython.display":
            raise ImportError

        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_ipython)

    with pytest.raises(RuntimeError, match="requires an IPython/Jupyter"):
        ig.show()


def test_private_sign_and_graph_validation_guards():
    ig = bt.bpy.ig.InfluenceGraph()

    with pytest.raises(ValueError, match="unsupported edge sign"):
        ig.add_edge("A", "B", sign=True)

    nx.MultiDiGraph.add_edge(ig, "A", "B")
    with pytest.raises(ValueError, match="missing edge attribute 'sign'"):
        ig._validate_graph()

    duplicate = bt.bpy.ig.InfluenceGraph()
    nx.MultiDiGraph.add_edge(duplicate, "A", "B", sign="+")
    nx.MultiDiGraph.add_edge(duplicate, "A", "B", sign=1)

    with pytest.raises(ValueError, match="duplicated edge sign"):
        duplicate._validate_graph()


def test_walks_from_returns_non_simple_walks_from_docstring():
    graph = nx.DiGraph()
    graph.add_edges_from(
        [
            ("A", "B"),
            ("B", "C"),
            ("C", "B"),
        ]
    )

    assert bt.bpy.ig.walks_from(graph, "A", max_depth=4) == [
        ["A", "B"],
        ["A", "B", "C"],
        ["A", "B", "C", "B"],
        ["A", "B", "C", "B", "C"],
    ]
    assert list(nx.all_simple_paths(graph, "A", "C")) == [["A", "B", "C"]]
    assert bt.bpy.ig.walks_from(graph, "A", max_depth=0) == []
