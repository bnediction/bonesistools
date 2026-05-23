#!/usr/bin/env python

import pytest

from boolean import BooleanAlgebra

import bonesistools as bt


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
            "A": ba.TRUE,
            "B": ba.FALSE,
        },
        ba=ba,
    )

    assert bn["A"] is ba.TRUE
    assert bn["B"] is ba.FALSE


def test_boolean_network_rejects_invalid_rule_type():

    with pytest.raises(TypeError):
        bt.bpy.bn.BooleanNetwork({"A": object()})


def test_boolean_network_validity():
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
        bn[1] = 0


def test_boolean_network_string_representation():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert str(bn) == "A <- B & ~C\nB <- 0\nC <- 1"


def test_boolean_network_repr_is_string_representation():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert repr(bn) == str(bn)


def test_boolean_network_copy_preserves_type_algebra_and_unchecked_rules():
    bn = bt.bpy.bn.BooleanNetwork({"A": "B"}, check=False)

    copied = bn.copy()

    assert isinstance(copied, bt.bpy.bn.BooleanNetwork)
    assert copied == bn
    assert copied is not bn
    assert copied.ba is bn.ba
    assert copied.undefined_symbols == {"B"}


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

    with pytest.raises(TypeError, match="unsupported argument type for 'old_name'"):
        bn.rename(1, "Y")

    with pytest.raises(TypeError, match="unsupported argument type for 'new_name'"):
        bn.rename("A", 1)


def test_boolean_network_from_bnet_and_to_bnet_file(tmp_path):
    infile = tmp_path / "network.bnet"
    outfile = tmp_path / "roundtrip.bnet"
    infile.write_text("# ignored\n\nA, B\nB, 1\n")

    bn = bt.bpy.bn.BooleanNetwork.from_bnet(infile)

    assert bn.rules == {"A": "B", "B": "1"}
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
        bn.equivalent(bn, method="unknown")


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


def test_boolean_network_to_networkx():

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_networkx()

    assert set(graph.nodes) == {"A", "B", "C"}
    assert graph.has_edge("B", "A")
    assert graph.has_edge("C", "A")

    assert graph["B"]["A"][0]["sign"] == 1
    assert graph["C"]["A"][0]["sign"] == -1


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

    assert dot.get_rankdir() == "LR"

    edges = {
        (edge.get_source().strip('"'), edge.get_destination().strip('"')): edge
        for edge in dot.get_edges()
    }

    assert ("B", "A") in edges
    assert ("C", "A") in edges

    assert edges[("B", "A")].get_color() == "green4"
    assert edges[("B", "A")].get_arrowhead() == "normal"

    assert edges[("C", "A")].get_color() == "red2"
    assert edges[("C", "A")].get_arrowhead() == "tee"
