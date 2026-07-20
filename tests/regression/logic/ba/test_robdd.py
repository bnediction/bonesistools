#!/usr/bin/env python

from itertools import product
from typing import Any, cast

import pytest
from boolean import BooleanAlgebra

import bonesistools as bt


def test_robdd_disjunction_structure_and_evaluation():
    robdd = bt.logic.ba.ROBDD("A | B")

    assert repr(robdd) == "ROBDD(variables=('A', 'B'), n_nodes=2)"
    assert robdd.variables == ("A", "B")
    assert robdd.n_nodes == 2
    assert robdd.evaluate({"A": 0, "B": 0}) is False
    assert robdd.evaluate({"A": 0, "B": 1}) is True
    assert robdd.evaluate({"A": 1, "B": 0, "unused": 2}) is True


def test_robdd_can_take_value_in_partial_assignment():
    robdd = bt.logic.ba.ROBDD("A | B")

    assert robdd.can_take_value(configuration={"A": 0}, value=0)
    assert robdd.can_take_value({"A": 0}, 1)
    assert not robdd.can_take_value({"A": 0, "B": 0}, 1)
    assert robdd.can_take_value({"A": "*", "B": 0}, 1)


def test_robdd_configurations_and_count():
    robdd = bt.logic.ba.ROBDD("A | B")

    assert robdd.configurations(1).enumerate() == (
        {"A": 0, "B": 1},
        {"A": 1, "B": 0},
        {"A": 1, "B": 1},
    )
    assert robdd.configurations(0).enumerate() == ({"A": 0, "B": 0},)
    assert robdd.count(1) == 3
    assert robdd.count(0) == 1


def test_robdd_configurations_builds_disjoint_encoding_without_add(monkeypatch):
    robdd = bt.logic.ba.ROBDD("(A & B) | (~A & C)")

    def fail(*args, **kwargs):
        raise AssertionError("ROBDD paths are already disjoint")

    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "add", fail)

    assert robdd.configurations().enumerate() == (
        {"A": 0, "B": 0, "C": 1},
        {"A": 0, "B": 1, "C": 1},
        {"A": 1, "B": 1, "C": 0},
        {"A": 1, "B": 1, "C": 1},
    )


def test_robdd_constants():
    false = bt.logic.ba.ROBDD(0)
    true = bt.logic.ba.ROBDD(1)

    assert false.variables == ()
    assert false.n_nodes == 0
    assert false.count(0) == 1
    assert false.count(1) == 0
    assert true.count(0) == 0
    assert true.count(1) == 1
    assert true.configurations().enumerate() == ({},)


def test_robdd_variable_order_changes_structure():
    rule = "(A & B) | (~A & C)"
    lexicographic = bt.logic.ba.ROBDD(rule)
    reordered = bt.logic.ba.ROBDD(rule, order=("B", "A", "C"))

    assert lexicographic.variables == ("A", "B", "C")
    assert reordered.variables == ("B", "A", "C")
    assert lexicographic.n_nodes == 3
    assert reordered.n_nodes == 4
    assert lexicographic.count() == reordered.count() == 4


def test_robdd_boolean_operators_combine_independent_diagrams():
    left = bt.logic.ba.ROBDD("A")
    right = bt.logic.ba.ROBDD("B")
    conjunction = left & right
    disjunction = left | right
    exclusive_disjunction = left ^ right
    negation = ~left

    assert conjunction.variables == ("A", "B")
    assert disjunction.variables == ("A", "B")
    assert exclusive_disjunction.variables == ("A", "B")
    assert negation.variables == ("A",)

    for a, b in product((0, 1), repeat=2):
        configuration = {"A": a, "B": b}
        assert conjunction.evaluate(configuration) is bool(a and b)
        assert disjunction.evaluate(configuration) is bool(a or b)
        assert exclusive_disjunction.evaluate(configuration) is bool(a ^ b)
        assert negation.evaluate(configuration) is (not bool(a))


def test_robdd_boolean_operators_are_immutable_and_reduced():
    left = bt.logic.ba.ROBDD("A | B")
    right = bt.logic.ba.ROBDD("A & B")
    left_representation = repr(left)
    right_representation = repr(right)

    contradiction = left ^ left
    tautology = left | ~left
    _ = (left & right) | (left ^ right)

    assert repr(left) == left_representation
    assert repr(right) == right_representation
    assert contradiction.variables == ("A", "B")
    assert contradiction.n_nodes == 0
    assert contradiction.count(0) == 4
    assert tautology.variables == ("A", "B")
    assert tautology.n_nodes == 0
    assert tautology.count(1) == 4


def test_robdd_boolean_operators_handle_constants():
    variable = bt.logic.ba.ROBDD("A")
    false = bt.logic.ba.ROBDD(False)
    true = bt.logic.ba.ROBDD(True)

    assert (variable & true).configurations().enumerate() == ({"A": 1},)
    assert (variable | false).configurations().enumerate() == ({"A": 1},)
    assert (variable ^ true).configurations().enumerate() == ({"A": 0},)
    assert (~false).configurations().enumerate() == ({},)


def test_robdd_boolean_operators_merge_compatible_custom_orders():
    left = bt.logic.ba.ROBDD("A | C", order=("C", "A"))
    right = bt.logic.ba.ROBDD("A & B", order=("B", "A"))
    result = left ^ right

    assert result.variables == ("B", "C", "A")
    for a, b, c in product((0, 1), repeat=3):
        configuration = {"A": a, "B": b, "C": c}
        assert result.evaluate(configuration) is bool((a or c) ^ (a and b))


def test_robdd_boolean_operators_reject_incompatible_custom_orders():
    left = bt.logic.ba.ROBDD("A & B", order=("A", "B"))
    right = bt.logic.ba.ROBDD("A | B", order=("B", "A"))

    with pytest.raises(ValueError, match="incompatible variable orders"):
        left & right


def test_robdd_accepts_boolean_expression():
    ba = BooleanAlgebra()
    robdd = bt.logic.ba.ROBDD(ba.parse("A & ~B"))

    assert robdd.evaluate({"A": 1, "B": 0})


def test_robdd_to_networkx_preserves_decisions_and_branches():
    graph = bt.logic.ba.ROBDD("A | B").to_networkx()
    root = graph.graph["root"]
    low = next(
        target for _, target, data in graph.edges(root, data=True) if not data["value"]
    )
    high = next(
        target for _, target, data in graph.edges(root, data=True) if data["value"]
    )

    assert graph.nodes[root]["variable"] == "A"
    assert graph.nodes[low]["variable"] == "B"
    assert graph.nodes[high] == {"terminal": True, "value": 1, "label": "1"}


def test_robdd_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="duplicated variables"):
        bt.logic.ba.ROBDD("A | B", order=("A", "A"))

    with pytest.raises(ValueError, match=r"missing=\['B'\], unknown=\['C'\]"):
        bt.logic.ba.ROBDD("A | B", order=("A", "C"))

    with pytest.raises(TypeError, match="only strings"):
        bt.logic.ba.ROBDD("A | B", order=cast(Any, ("A", 1)))

    robdd = bt.logic.ba.ROBDD("A | B")
    with pytest.raises(ValueError, match="fixed values for all ROBDD variables"):
        robdd.evaluate({"A": 0})

    with pytest.raises(ValueError, match="value must be 0 or 1"):
        robdd.count(2)
