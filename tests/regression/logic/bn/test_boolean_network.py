#!/usr/bin/env python

from typing import Any, cast

import pytest
from boolean import BooleanAlgebra, Expression

import bonesistools as bt
from bonesistools.logic.boolean_network import _typing


def _edge_data(graph, source, target, key=0):
    return cast(Any, graph[source][target])[key]


def test_boolean_network_coerces_string_rules():

    ba = BooleanAlgebra()
    bn = bt.logic.bn.BooleanNetwork(
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

    bn = bt.logic.bn.BooleanNetwork({"A": True, "B": False})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_integer_constants():

    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": 0})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_string_constants():

    bn = bt.logic.bn.BooleanNetwork({"A": "1", "B": "0"})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_accepts_boolean_algebra_constants():

    ba = BooleanAlgebra()

    bn = bt.logic.bn.BooleanNetwork(
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
        bt.logic.bn.BooleanNetwork({"A": cast(Any, object())})


def test_boolean_network_rejects_undefined_symbols():
    with pytest.raises(ValueError):
        bt.logic.bn.BooleanNetwork({"A": "B"})


def test_boolean_network_can_be_created_unchecked():

    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    assert bn.components == {"A"}
    assert bn.symbols == {"B"}
    assert bn.undefined_symbols == {"B"}
    assert not bn.is_closed


def test_boolean_network_validate_raises_on_undefined_symbols():

    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError):
        bn.validate()


def test_boolean_network_components():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.components == {"A", "B", "C"}


def test_boolean_network_symbols():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.symbols == {"B", "C"}


def test_boolean_network_undefined_symbols_empty_for_closed_network():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.undefined_symbols == set()
    assert bn.is_closed


def test_boolean_network_rules_property():

    bn = bt.logic.bn.BooleanNetwork(
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

    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    bn["B"] = "A & ~C"
    bn["C"] = 0

    assert bn.rule("B") == "A & ~C"
    assert bn["C"] is bn.ba.FALSE

    with pytest.raises(TypeError, match="unsupported argument type for 'component'"):
        bn[cast(Any, 1)] = 0


def test_boolean_network_string_representation():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert str(bn) == "A <- B & ~C\nB <- 0\nC <- 1"


def test_boolean_network_repr_is_compact():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert repr(bn) == "BooleanNetwork(components=3)"
    assert "\n" not in repr(bn)


def test_boolean_network_copy_preserves_type_algebra_and_unchecked_rules():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    copied = bn.copy()

    assert isinstance(copied, bt.logic.bn.BooleanNetwork)
    assert copied == bn
    assert copied is not bn
    assert copied.ba is bn.ba
    assert copied.undefined_symbols == {"B"}


def test_is_boolean_network_like_validates_mapping_rules():
    assert _typing.is_boolean_network_like({"A": "B", "B": 1})
    assert _typing.is_boolean_network_like(bt.logic.bn.BooleanNetwork({"A": 1}))

    assert not _typing.is_boolean_network_like({"A": object()})
    assert not _typing.is_boolean_network_like({1: "A"})
    assert not _typing.is_boolean_network_like(object())


def test_boolean_network_typing_namespace_is_not_exposed_publicly():
    assert "typing" not in dir(bt.logic.bn)

    with pytest.raises(AttributeError):
        getattr(bt.logic.bn, "typing")


def test_boolean_network_to_bnet_returns_string():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.to_bnet() == "A, B&!C\nB, 0\nC, 1\n"


def test_boolean_network_convert_to_mpbn():
    mpbn = pytest.importorskip("mpbn")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    converted = bn.convert("mpbn")

    assert isinstance(converted, mpbn.MPBooleanNetwork)
    assert set(converted) == {"A", "B", "C"}


def test_boolean_network_convert_to_biolqm_supports_model_analysis():
    biolqm = pytest.importorskip("biolqm")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    converted = bn.convert("biolqm")
    fixpoints = biolqm.fixpoints(converted)

    assert len(fixpoints) == 1
    assert dict(fixpoints[0]) == {"A": 0, "B": 0, "C": 1}


def test_boolean_network_convert_to_pyboolnet_uses_prime_implicants():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.convert("pyboolnet") == {
        "A": [
            [{"B": 0}, {"C": 1}],
            [{"B": 1, "C": 0}],
        ],
        "B": [[{}], []],
        "C": [[], [{}]],
    }


def test_boolean_network_convert_to_pyboolnet_matches_pyboolnet_primes():
    file_exchange = pytest.importorskip("pyboolnet.file_exchange")
    pyboolnet_primes = pytest.importorskip("pyboolnet.prime_implicants")
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "(B & ~C) | (~B & C)",
            "B": "A | B",
            "C": 0,
            "D": 1,
        }
    )

    converted = bn.convert("pyboolnet")
    expected = file_exchange.bnet2primes(cast(str, bn.to_bnet()))

    assert pyboolnet_primes.primes_are_equal(converted, expected)


def test_boolean_network_convert_to_minibn():
    minibn = pytest.importorskip("colomoto.minibn")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    converted = bn.convert("minibn")

    assert isinstance(converted, minibn.BooleanNetwork)
    assert set(converted) == {"A", "B", "C"}


def test_boolean_network_convert_validates_target():
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(TypeError, match="unsupported argument type for 'target'"):
        bn.convert(cast(Any, 1))

    with pytest.raises(ValueError, match="unsupported conversion target"):
        bn.convert(cast(Any, "unknown"))

    with pytest.raises(ValueError, match="unsupported conversion target"):
        bn.convert("minibn.BooleanNetwork")


def test_boolean_network_rename_validates_inputs_and_collisions():
    bn = bt.logic.bn.BooleanNetwork(
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
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "C"}, check=False)
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="undefined components"):
        bn.rename("A", "X")

    assert bn.rules == rules


def test_boolean_network_rename_accepts_named_old_new_arguments():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    bn.rename(old="A", new="X")

    assert bn.rules == {"X": "B", "B": "1"}


def test_boolean_network_relabel_renames_several_components_and_ignores_missing():
    bn = bt.logic.bn.BooleanNetwork({"Trp53": "Sox2", "Sox2": 1})

    assert bn.relabel({"Trp53": "TP53", "Sox2": "SOX2", "missing": "X"}) is None

    assert bn.rules == {"TP53": "SOX2", "SOX2": "1"}


def test_boolean_network_relabel_supports_atomic_swaps():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    bn.relabel({"A": "B", "B": "A"})

    assert bn.rules == {"B": "A", "A": "B"}


def test_boolean_network_relabel_rejects_component_merges_without_mutating():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="merge Boolean network components"):
        bn.relabel({"A": "B"})

    assert bn.rules == rules


def test_boolean_network_relabel_validates_mapping():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    with pytest.raises(TypeError, match="unsupported argument type for 'mapping'"):
        bn.relabel(cast(Any, [("A", "X")]))

    with pytest.raises(TypeError, match="unsupported mapping key type"):
        bn.relabel(cast(Any, {1: "X"}))

    with pytest.raises(TypeError, match="unsupported mapping value type"):
        bn.relabel(cast(Any, {"A": 1}))

    assert bn.rules == {"A": "B", "B": "1"}


def test_boolean_network_to_bnet_file(tmp_path):
    outfile = tmp_path / "roundtrip.bnet"
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    assert bn.rules == {"A": "B", "B": "1"}
    assert bn.to_bnet(outfile) is None
    assert outfile.read_text() == "A, B\nB, 1\n"


def test_deprecated_read_bnet_routes_to_io(tmp_path):
    infile = tmp_path / "network.bnet"
    infile.write_text("A, B\nB, 1\n")

    with pytest.warns(FutureWarning, match="bt.logic.bn.read_bnet"):
        bn = bt.logic.bn.read_bnet(infile)

    assert bn.rules == bt.logic.io.read_bnet(infile).rules


def test_boolean_network_structural_equality():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "C | B",
            "B": 0,
            "C": 1,
        }
    )

    assert bn1 == bn2


def test_boolean_network_structural_inequality():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & (C | ~C)",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 0,
        }
    )

    assert bn1 != bn2
    assert bn1.__eq__(object()) is NotImplemented
    assert bn1.__ne__(object()) is NotImplemented


def test_boolean_network_equivalence_truth_table_only():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "(B & C) | (~B & D) | (C & D)",
            "B": 0,
            "C": 1,
            "D": 0,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
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

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    assert not bn1.equivalent(bn2, method="truth_table")


def test_boolean_network_equivalence_requires_same_components():

    bn1 = bt.logic.bn.BooleanNetwork({"A": 1})
    bn2 = bt.logic.bn.BooleanNetwork({"A": 1, "B": 0})

    assert not bn1.equivalent(bn2)
    assert bn1.equivalent(object()) is NotImplemented


def test_boolean_network_equivalence_rejects_unknown_method():

    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(ValueError):
        bn.equivalent(bn, method=cast(Any, "unknown"))


def test_boolean_network_influences():

    bn = bt.logic.bn.BooleanNetwork(
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


def test_boolean_network_to_influence_graph():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_influence_graph()

    assert isinstance(graph, bt.logic.ig.InfluenceGraph)
    assert set(graph.nodes) == {"A", "B", "C"}
    assert graph.has_edge("B", "A")
    assert graph.has_edge("C", "A")

    assert _edge_data(graph, "B", "A")["sign"] == 1
    assert _edge_data(graph, "C", "A")["sign"] == -1
