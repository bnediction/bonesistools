#!/usr/bin/env python

from collections import Counter

import pytest
from boolean import BooleanAlgebra

import bonesistools as bt


def test_partial_boolean_differential_uses_partial_boolean_order():
    differential = bt.bpy.ba.PartialBooleanDifferential

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pbstar = bt.bpy.ba.PartialBoolean("*")
    pb1 = bt.bpy.ba.PartialBoolean(1)

    assert differential.differential(pb0, pb0) == 0
    assert differential.differential(pbstar, pbstar) == 0
    assert differential.differential(pb1, pb1) == 0

    assert differential.differential(pb0, pbstar) == 1
    assert differential.differential(pbstar, pb1) == 1
    assert differential.differential(pb0, pb1) == 1

    assert differential.differential(pbstar, pb0) == -1
    assert differential.differential(pb1, pbstar) == -1
    assert differential.differential(pb1, pb0) == -1

    assert differential.differential(False, 0) == 0
    assert differential.differential(True, 1) == 0
    assert differential.differential(float("nan"), pbstar) == 0


@pytest.mark.parametrize(
    "source_value,target_values,sign,expected",
    [
        (1, (0, 1), 1, True),
        (1, (1, 0), 1, False),
        (0, (0, 1), 1, False),
        (0, (1, 0), 1, True),
        (1, (1, 0), -1, True),
        (1, (0, 1), -1, False),
        (0, (1, 0), -1, False),
        (0, (0, 1), -1, True),
    ],
)
def test_boolean_predecessor_inference_pairwise_table(
    source_value,
    target_values,
    sign,
    expected,
):
    inference = bt.bpy.ba.BooleanPredecessorInference

    assert (
        inference.pairwise_predecessor_test(
            source_value,
            source_value,
            target_values[0],
            target_values[1],
            sign,
        )
        is expected
    )


@pytest.mark.parametrize(
    "source_values,target_values,sign",
    [
        ((0, 1), (0, 1), 1),
        ((1, 0), (0, 1), 1),
        ((float("nan"), float("nan")), (0, 1), 1),
        ((1, 1), (0, 0), 1),
        ((0, 0), (1, 1), -1),
    ],
)
def test_boolean_predecessor_inference_pairwise_no_conclusion(
    source_values,
    target_values,
    sign,
):
    inference = bt.bpy.ba.BooleanPredecessorInference

    assert (
        inference.pairwise_predecessor_test(
            source_values[0],
            source_values[1],
            target_values[0],
            target_values[1],
            sign,
        )
        is None
    )


def test_predecessor_inference_scores_single_influence():
    inference = bt.bpy.ba.BooleanPredecessorInference

    cell1 = {"A": 1, "B": 0}
    cell2 = {"A": 1, "B": 1}
    interactions = [("A", "B", {"sign": 1})]

    assert inference.predecessor_votes(cell1, cell2, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 1,
            inference.SECOND_PRECEDES_FIRST: 0,
            inference.INCONCLUSIVE: 0,
        }
    )
    assert inference.predecessor(cell1, cell2, interactions) == (
        inference.FIRST_PRECEDES_SECOND
    )
    assert inference.predecessor_score(cell1, cell2, interactions) == 1.0

    assert inference.predecessor_votes(cell2, cell1, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 0,
            inference.SECOND_PRECEDES_FIRST: 1,
            inference.INCONCLUSIVE: 0,
        }
    )
    assert inference.predecessor(cell2, cell1, interactions) == (
        inference.SECOND_PRECEDES_FIRST
    )
    assert inference.predecessor_score(cell2, cell1, interactions) == -1.0


def test_boolean_predecessor_inference_high_level_scbooldiff_example():
    inference = bt.bpy.ba.BooleanPredecessorInference

    cell1 = {"A": 1, "B": 0, "C": 1, "D": 1}
    cell2 = {"A": 1, "B": 1, "C": 1, "D": 0}
    interactions = [
        ("A", "B", {"sign": 1}),
        ("C", "D", {"sign": -1}),
    ]

    assert inference.predecessor_votes(cell1, cell2, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 2,
            inference.SECOND_PRECEDES_FIRST: 0,
            inference.INCONCLUSIVE: 0,
        }
    )
    assert inference.predecessor(cell1, cell2, interactions) == (
        inference.FIRST_PRECEDES_SECOND
    )
    assert inference.predecessor_score(cell1, cell2, interactions) == 1.0


def test_boolean_predecessor_inference_counts_inconclusive_interactions():
    inference = bt.bpy.ba.BooleanPredecessorInference

    cell1 = {"A": 1, "B": 0, "C": "*", "D": 0}
    cell2 = {"A": 1, "B": 1, "C": "*", "D": 1}
    interactions = [
        ("A", "B", {"sign": 1}),
        ("C", "D", {"sign": 1}),
        ("A", "missing", {"sign": 1}),
        ("missing", "B", {"sign": 1}),
    ]

    assert inference.predecessor_votes(cell1, cell2, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 1,
            inference.SECOND_PRECEDES_FIRST: 0,
            inference.INCONCLUSIVE: 3,
        }
    )
    assert inference.predecessor(cell1, cell2, interactions) == (
        inference.FIRST_PRECEDES_SECOND
    )
    assert inference.predecessor_score(cell1, cell2, interactions) == 1.0


def test_predecessor_inference_returns_none_for_tied_votes():
    inference = bt.bpy.ba.BooleanPredecessorInference

    cell1 = {"A": 1, "B": 0, "C": 1, "D": 1}
    cell2 = {"A": 1, "B": 1, "C": 1, "D": 0}
    interactions = [
        ("A", "B", {"sign": 1}),
        ("C", "D", {"sign": 1}),
    ]

    assert inference.predecessor_votes(cell1, cell2, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 1,
            inference.SECOND_PRECEDES_FIRST: 1,
            inference.INCONCLUSIVE: 0,
        }
    )
    assert inference.predecessor(cell1, cell2, interactions) is None
    assert inference.predecessor_score(cell1, cell2, interactions) == 0.0


def test_predecessor_score_is_zero_without_conclusive_votes():
    inference = bt.bpy.ba.BooleanPredecessorInference

    cell1 = {"A": "*", "B": 0}
    cell2 = {"A": "*", "B": 1}
    interactions = [("A", "B", {"sign": 1})]

    assert inference.predecessor_votes(cell1, cell2, interactions) == Counter(
        {
            inference.FIRST_PRECEDES_SECOND: 0,
            inference.SECOND_PRECEDES_FIRST: 0,
            inference.INCONCLUSIVE: 1,
        }
    )
    assert inference.predecessor(cell1, cell2, interactions) is None
    assert inference.predecessor_score(cell1, cell2, interactions) == 0.0


def test_boolean_predecessor_inference_rejects_invalid_inputs():
    differential = bt.bpy.ba.PartialBooleanDifferential
    inference = bt.bpy.ba.BooleanPredecessorInference

    with pytest.raises(ValueError):
        differential.differential(2, 1)

    with pytest.raises(ValueError):
        differential.differential("0", 1)

    with pytest.raises(TypeError):
        differential.differential(object(), 1)

    with pytest.raises(ValueError):
        inference.pairwise_predecessor_test(1, 1, 0, 1, sign=0)

    with pytest.raises(ValueError):
        inference.predecessor_votes(
            {"A": 1, "B": 0},
            {"A": 1, "B": 1},
            [("A", "B", {"sign": 0})],
        )


def test_configuration_like_is_concrete_and_distinct_from_hypercube_like():
    assert bt.bpy.ba.is_configuration_like({"A": 0, "B": True})
    assert bt.bpy.ba.is_configuration_like({"A": bt.bpy.ba.PartialBoolean(1)})

    assert not bt.bpy.ba.is_configuration_like({"A": "*"})
    assert not bt.bpy.ba.is_configuration_like({"A": bt.bpy.ba.PartialBoolean("*")})
    assert not bt.bpy.ba.is_configuration_like({"A": 2})
    assert not bt.bpy.ba.is_configuration_like({1: 0})
    assert not bt.bpy.ba.is_configuration_like([("A", 0)])

    assert bt.bpy.ba.is_hypercube_like({"A": "*"})


def test_dnf_to_structure_conjunction():
    ba = BooleanAlgebra()
    expr = ba.parse("x & ~y")

    assert bt.bpy.ba.dnf_to_structure(ba, expr) == frozenset(
        {
            frozenset(
                {
                    ("x", True),
                    ("y", False),
                }
            )
        }
    )


def test_dnf_to_structure_disjunction():
    ba = BooleanAlgebra()
    expr = ba.parse("x | ~y")

    assert bt.bpy.ba.dnf_to_structure(ba, expr) == frozenset(
        {
            frozenset({("x", True)}),
            frozenset({("y", False)}),
        }
    )


def test_dnf_to_structure_constants():
    ba = BooleanAlgebra()

    assert bt.bpy.ba.dnf_to_structure(ba, ba.TRUE) is True
    assert bt.bpy.ba.dnf_to_structure(ba, ba.FALSE) is False


def test_dnf_to_structure_list_container():
    ba = BooleanAlgebra()
    expr = ba.parse("x | ~y")

    assert bt.bpy.ba.dnf_to_structure(ba, expr, container=list, sort=True) == [
        [("x", True)],
        [("y", False)],
    ]


def test_dnf_to_structure_rejects_non_dnf():
    ba = BooleanAlgebra()
    expr = ba.parse("x & (y | z)")

    with pytest.raises(ValueError):
        bt.bpy.ba.dnf_to_structure(ba, expr)
