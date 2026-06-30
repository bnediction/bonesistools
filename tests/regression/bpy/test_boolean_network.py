#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, cast

import pytest
from boolean import BooleanAlgebra, Expression

import bonesistools as bt


def _pydot_get(obj, method):
    return getattr(cast(Any, obj), method)()


def _pydot_get_string(obj, method):
    return cast(str, _pydot_get(obj, method)).strip('"')


def _edge_data(graph, source, target, key=0):
    return cast(Any, graph[source][target])[key]


def test_boolean_network_coerces_string_rules():

    ba = BooleanAlgebra()
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & !C",
            "B": 0,
            "C": 1,
        },
        ba=ba,
    )

    assert bn.rule("A") == "B & ~C"
    assert bn.rule("B") == "0"
    assert bn.rule("C") == "1"


def test_boolean_network_coerces_boolean_constants():

    bn = bt.bpy.bn.BooleanNetwork({"A": True, "B": False})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_integer_constants():

    bn = bt.bpy.bn.BooleanNetwork({"A": 1, "B": 0})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_string_constants():

    bn = bt.bpy.bn.BooleanNetwork({"A": "1", "B": "0"})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_accepts_boolean_algebra_constants():

    ba = BooleanAlgebra()

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": cast(Expression, ba.TRUE),
            "B": cast(Expression, ba.FALSE),
        },
        ba=ba,
    )

    assert bn["A"] is ba.TRUE
    assert bn["B"] is ba.FALSE


def test_boolean_network_rejects_invalid_rule_type():

    with pytest.raises(TypeError):
        bt.bpy.bn.BooleanNetwork({"A": cast(Any, object())})


def test_boolean_network_rejects_undefined_symbols():
    with pytest.raises(ValueError):
        bt.bpy.bn.BooleanNetwork({"A": "B"})


def test_boolean_network_can_be_created_unchecked():

    bn = bt.bpy.bn.BooleanNetwork({"A": "B"}, check=False)

    assert bn.components == {"A"}
    assert bn.symbols == {"B"}
    assert bn.undefined_symbols == {"B"}
    assert not bn.is_closed


def test_boolean_network_validate_raises_on_undefined_symbols():

    bn = bt.bpy.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError):
        bn.validate()


def test_boolean_network_components():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.components == {"A", "B", "C"}


def test_boolean_network_symbols():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.symbols == {"B", "C"}


def test_boolean_network_undefined_symbols_empty_for_closed_network():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.undefined_symbols == set()
    assert bn.is_closed


def test_boolean_network_rules_property():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.rules == {
        "A": "B & ~C",
        "B": "0",
        "C": "1",
    }


def test_boolean_network_setitem_coerces_rules():

    bn = bt.bpy.bn.BooleanNetwork({"A": 1})

    bn["B"] = "A & ~C"
    bn["C"] = 0

    assert bn.rule("B") == "A & ~C"
    assert bn["C"] is bn.ba.FALSE

    with pytest.raises(TypeError, match="unsupported argument type for 'component'"):
        bn[cast(Any, 1)] = 0


def test_boolean_network_string_representation():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert str(bn) == "A <- B & ~C\nB <- 0\nC <- 1"


def test_boolean_network_repr_is_compact():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert repr(bn) == "BooleanNetwork(components=3)"
    assert "\n" not in repr(bn)


def test_boolean_network_copy_preserves_type_algebra_and_unchecked_rules():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B"}, check=False)

    copied = bn.copy()

    assert isinstance(copied, bt.bpy.bn.BooleanNetwork)
    assert copied == bn
    assert copied is not bn
    assert copied.ba is bn.ba
    assert copied.undefined_symbols == {"B"}


def test_is_boolean_network_like_validates_mapping_rules():
    assert bt.bpy.bn.typing.is_boolean_network_like({"A": "B", "B": 1})
    assert bt.bpy.bn.typing.is_boolean_network_like(bt.bpy.bn.BooleanNetwork({"A": 1}))

    assert not bt.bpy.bn.typing.is_boolean_network_like({"A": object()})
    assert not bt.bpy.bn.typing.is_boolean_network_like({1: "A"})
    assert not bt.bpy.bn.typing.is_boolean_network_like(object())


def test_boolean_network_to_bnet_returns_string():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.to_bnet() == "A, B&!C\nB, 0\nC, 1\n"


def test_boolean_network_rename_validates_inputs_and_collisions():
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "A & C",
            "C": 1,
        }
    )

    assert bn.rename("A", "A") is None
    assert bn.rules["A"] == "B"

    with pytest.raises(KeyError, match="component 'missing' not found"):
        bn.rename("missing", "Y")

    with pytest.raises(ValueError, match="already exists"):
        bn.rename("A", "B")

    with pytest.raises(TypeError, match="unsupported argument type for 'old'"):
        bn.rename(cast(Any, 1), "Y")

    with pytest.raises(TypeError, match="unsupported argument type for 'new'"):
        bn.rename("A", cast(Any, 1))


def test_boolean_network_rename_validates_candidate_before_mutating():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": "C"}, check=False)
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="undefined components"):
        bn.rename("A", "X")

    assert bn.rules == rules


def test_boolean_network_rename_accepts_named_old_new_arguments():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})

    bn.rename(old="A", new="X")

    assert bn.rules == {"X": "B", "B": "1"}


def test_boolean_network_relabel_renames_several_components_and_ignores_missing():
    bn = bt.bpy.bn.BooleanNetwork({"Trp53": "Sox2", "Sox2": 1})

    assert bn.relabel({"Trp53": "TP53", "Sox2": "SOX2", "missing": "X"}) is None

    assert bn.rules == {"TP53": "SOX2", "SOX2": "1"}


def test_boolean_network_relabel_supports_atomic_swaps():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": "A"})

    bn.relabel({"A": "B", "B": "A"})

    assert bn.rules == {"B": "A", "A": "B"}


def test_boolean_network_relabel_rejects_component_merges_without_mutating():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="merge Boolean network components"):
        bn.relabel({"A": "B"})

    assert bn.rules == rules


def test_boolean_network_relabel_validates_mapping():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})

    with pytest.raises(TypeError, match="unsupported argument type for 'mapping'"):
        bn.relabel(cast(Any, [("A", "X")]))

    with pytest.raises(TypeError, match="unsupported mapping key type"):
        bn.relabel(cast(Any, {1: "X"}))

    with pytest.raises(TypeError, match="unsupported mapping value type"):
        bn.relabel(cast(Any, {"A": 1}))

    assert bn.rules == {"A": "B", "B": "1"}


def test_boolean_network_from_bnet_and_to_bnet_file(tmp_path):
    infile = tmp_path / "network.bnet"
    outfile = tmp_path / "roundtrip.bnet"
    infile.write_text("# ignored\n\nA, B\nB, 1\n")

    bn = bt.bpy.bn.BooleanNetwork.from_bnet(infile)
    read_bn = bt.bpy.bn.read_bnet(infile)

    assert bn.rules == {"A": "B", "B": "1"}
    assert read_bn.rules == bn.rules
    assert bn.to_bnet(outfile) is None
    assert outfile.read_text() == "A, B\nB, 1\n"


def test_boolean_network_structural_equality():

    bn1 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "C | B",
            "B": 0,
            "C": 1,
        }
    )

    assert bn1 == bn2


def test_boolean_network_structural_inequality():

    bn1 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & (C | ~C)",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 0,
        }
    )

    assert bn1 != bn2
    assert bn1.__eq__(object()) is NotImplemented
    assert bn1.__ne__(object()) is NotImplemented


def test_boolean_network_equivalence_truth_table_only():

    bn1 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "(B & C) | (~B & D) | (C & D)",
            "B": 0,
            "C": 1,
            "D": 0,
        }
    )

    bn2 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "(B & C) | (~B & D)",
            "B": 0,
            "C": 1,
            "D": 0,
        }
    )

    assert not bn1.equivalent(bn2, method="simplify")
    assert bn1.equivalent(bn2, method="truth_table")


def test_boolean_network_not_equivalent_truth_table():

    bn1 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    assert not bn1.equivalent(bn2, method="truth_table")


def test_boolean_network_equivalence_requires_same_components():

    bn1 = bt.bpy.bn.BooleanNetwork({"A": 1})
    bn2 = bt.bpy.bn.BooleanNetwork({"A": 1, "B": 0})

    assert not bn1.equivalent(bn2)
    assert bn1.equivalent(object()) is NotImplemented


def test_boolean_network_equivalence_rejects_unknown_method():

    bn = bt.bpy.bn.BooleanNetwork({"A": 1})

    with pytest.raises(ValueError):
        bn.equivalent(bn, method=cast(Any, "unknown"))


def test_boolean_network_influences():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.influences() == {
        ("B", "A", 1),
        ("C", "A", -1),
    }


def test_boolean_network_fixed_points_and_predicate():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.fixed_points() == [
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
    ]
    assert bn.fixed_points(limit=1) == [{"A": 0, "B": 0}]
    assert bn.is_fixed_point({"A": True, "B": True}) is True
    assert bn.is_fixed_point({"A": 1, "B": 0}) is False

    no_fixed_point = bt.bpy.bn.BooleanNetwork({"A": "~A"})
    assert no_fixed_point.fixed_points() == []


def test_boolean_network_next_methods():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B & ~C", "B": "A", "C": 0})

    assert bn.next_state("A", {"A": 0, "B": 1, "C": 0}) == 1
    assert bn.next_state("A", {"A": 0, "B": 1, "C": 1}) == 0

    assert bn.next_configuration({"A": 1, "B": 1, "C": 0}) == {
        "A": 1,
        "B": 1,
        "C": 0,
    }
    assert bn.next_configuration({"A": 0, "B": 1, "C": 1}) == {
        "A": 0,
        "B": 0,
        "C": 0,
    }

    with pytest.raises(ValueError, match="expected components"):
        bn.next_state("A", {"A": 0, "B": 1, "C": 0, "D": 1})

    unchecked = bt.bpy.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="still depends on"):
        unchecked.next_state("A", {"A": 0})

    with pytest.raises(ValueError, match="still depends on"):
        unchecked.next_configuration({"A": 0})


def test_boolean_network_fixed_points_validate_inputs():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'limit'"):
        bn.fixed_points(limit=-1)

    with pytest.raises(ValueError, match="expected components"):
        bn.is_fixed_point({"A": 1})

    with pytest.raises(ValueError, match="expected components"):
        bn.is_fixed_point({"A": 1, "B": 1, "C": 0})

    assert bn.is_fixed_point({"A": 1, "B": "*"}) is False

    assert bn.is_fixed_point({"A": 1, "B": bt.bpy.ba.PartialBoolean("*")}) is False

    with pytest.raises(TypeError, match="unsupported argument type for 'state'"):
        bn.is_fixed_point(cast(Any, object()))

    with pytest.raises(ValueError, match="expected string component names"):
        bn.is_fixed_point(cast(Any, {"A": 1, 2: 0}))


def test_boolean_network_to_influence_graph():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_influence_graph()

    assert isinstance(graph, bt.bpy.ig.InfluenceGraph)
    assert set(graph.nodes) == {"A", "B", "C"}
    assert graph.has_edge("B", "A")
    assert graph.has_edge("C", "A")

    assert _edge_data(graph, "B", "A")["sign"] == 1
    assert _edge_data(graph, "C", "A")["sign"] == -1


def test_boolean_network_to_graphviz(fake_graphviz):
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_graphviz(rankdir="LR")

    edges = {(source, target): attrs for source, target, attrs in graph.edges}

    assert isinstance(graph, fake_graphviz)
    assert graph.graph_attr["rankdir"] == "LR"
    assert ("B", "A") in edges
    assert ("C", "A") in edges
    assert edges[("B", "A")]["color"] == "green4"
    assert edges[("B", "A")]["arrowhead"] == "normal"
    assert edges[("C", "A")]["color"] == "red2"
    assert edges[("C", "A")]["arrowhead"] == "tee"


def test_boolean_network_to_pydot():

    pytest.importorskip("pydot")

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    dot = bn.to_pydot(rankdir="LR")

    assert cast(Any, dot).get_rankdir() == "LR"

    edges = {
        (
            _pydot_get_string(edge, "get_source"),
            _pydot_get_string(edge, "get_destination"),
        ): edge
        for edge in dot.get_edges()
    }

    assert ("B", "A") in edges
    assert ("C", "A") in edges

    assert _pydot_get(edges[("B", "A")], "get_color") == "green4"
    assert _pydot_get(edges[("B", "A")], "get_arrowhead") == "normal"

    assert _pydot_get(edges[("C", "A")], "get_color") == "red2"
    assert _pydot_get(edges[("C", "A")], "get_arrowhead") == "tee"


def test_boolean_network_show(monkeypatch):
    calls = {}

    class FakeDot:
        def create_svg(self):
            return b'<svg width="100pt" height="200pt"></svg>'

    def fake_to_pydot(self, program="dot", edge_style=None, **kwargs):
        calls["to_pydot"] = {
            "program": program,
            "edge_style": edge_style,
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
    monkeypatch.setattr(bt.bpy.bn.BooleanNetwork, "to_pydot", fake_to_pydot)

    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})

    bn.show(
        program="neato",
        edge_style=lambda data: data,
        width=600,
        rankdir="LR",
    )

    assert calls["to_pydot"]["program"] == "neato"
    assert calls["to_pydot"]["kwargs"] == {"rankdir": "LR"}
    assert isinstance(calls["display"], FakeSVG)
    assert calls["display"].svg == '<svg width="600"></svg>'
