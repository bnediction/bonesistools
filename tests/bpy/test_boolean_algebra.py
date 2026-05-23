#!/usr/bin/env python

import pytest

from boolean import BooleanAlgebra

import bonesistools as bt


def test_boolean_differential_calculus_differential_uses_partial_boolean_order():
    calculus = bt.bpy.ba.BooleanDifferentialCalculus()

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pbstar = bt.bpy.ba.PartialBoolean("*")
    pb1 = bt.bpy.ba.PartialBoolean(1)

    assert calculus.differential(pb0, pb0) == 0
    assert calculus.differential(pbstar, pbstar) == 0
    assert calculus.differential(pb1, pb1) == 0

    assert calculus.differential(pb0, pbstar) == 1
    assert calculus.differential(pbstar, pb1) == 1
    assert calculus.differential(pb0, pb1) == 1

    assert calculus.differential(pbstar, pb0) == -1
    assert calculus.differential(pb1, pbstar) == -1
    assert calculus.differential(pb1, pb0) == -1

    assert calculus.differential(False, 0) == 0
    assert calculus.differential(True, 1) == 0
    assert calculus.differential(float("nan"), pbstar) == 0


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
def test_boolean_differential_calculus_pairwise_predecessor_table(
    source_value,
    target_values,
    sign,
    expected,
):
    calculus = bt.bpy.ba.BooleanDifferentialCalculus()

    assert (
        calculus.pairwise_predecessor_test(
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
def test_boolean_differential_calculus_pairwise_predecessor_no_conclusion(
    source_values,
    target_values,
    sign,
):
    calculus = bt.bpy.ba.BooleanDifferentialCalculus()

    assert (
        calculus.pairwise_predecessor_test(
            source_values[0],
            source_values[1],
            target_values[0],
            target_values[1],
            sign,
        )
        is None
    )


def test_boolean_differential_calculus_scores_cell_order_from_single_influence():
    calculus = bt.bpy.ba.BooleanDifferentialCalculus()
    interactions = [("A", "B", {"sign": 1})]

    def predecessor_score(first_cell, second_cell):
        score = 0

        for source, target, data in interactions:
            conclusion = calculus.pairwise_predecessor_test(
                first_cell[source],
                second_cell[source],
                first_cell[target],
                second_cell[target],
                data["sign"],
            )

            if conclusion is True:
                score += 1
            elif conclusion is False:
                score -= 1

        return score

    cell1 = {"A": 1, "B": 0}
    cell2 = {"A": 1, "B": 1}

    assert predecessor_score(cell1, cell2) == 1
    assert predecessor_score(cell2, cell1) == -1


def test_boolean_differential_calculus_rejects_invalid_inputs():
    calculus = bt.bpy.ba.BooleanDifferentialCalculus()

    with pytest.raises(ValueError):
        calculus.differential(2, 1)

    with pytest.raises(TypeError):
        calculus.differential("0", 1)

    with pytest.raises(ValueError):
        calculus.pairwise_predecessor_test(1, 1, 0, 1, sign=0)


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
