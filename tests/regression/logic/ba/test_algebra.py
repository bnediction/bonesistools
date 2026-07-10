#!/usr/bin/env python

from collections import Counter

import pytest

import bonesistools as bt


def test_static_boolean_algebra_utilities_cannot_be_instantiated():
    with pytest.raises(TypeError, match="BooleanPredecessorInference"):
        bt.logic.ba.BooleanPredecessorInference()


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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

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
    inference = bt.logic.ba.BooleanPredecessorInference

    with pytest.raises(ValueError):
        inference.pairwise_predecessor_test(1, 1, 0, 1, sign=0)

    with pytest.raises(ValueError):
        inference.predecessor_votes(
            {"A": 1, "B": 0},
            {"A": 1, "B": 1},
            [("A", "B", {"sign": 0})],
        )
