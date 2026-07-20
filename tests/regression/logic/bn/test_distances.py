#!/usr/bin/env python

from typing import Any, cast

import pytest
from boolean import BooleanAlgebra, Expression
from boolean.boolean import _FALSE

import bonesistools as bt
from bonesistools.logic.boolean_network import _distances


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_identical_component_sets_are_independent_of_order(metric):

    bn1 = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    bn2 = bt.logic.bn.BooleanNetwork({"B": "A", "A": "B"})

    assert bt.logic.bn.distance(bn1, bn2, metric=metric) == 0.0
    assert bt.logic.bn.similarity(bn1, bn2, metric=metric) == 1.0


def test_network_like_rule_mappings_are_supported():

    bn1 = {"A": "A & B", "B": "B"}
    bn2 = {"B": "B", "A": "A"}

    assert bt.logic.bn.distance(bn1, bn2, metric="hamming") == 1 / 8
    assert bt.logic.bn.similarity(bn1, bn2, metric="hamming") == 7 / 8


@pytest.mark.parametrize(
    "comparison",
    [bt.logic.bn.distance, bt.logic.bn.similarity],
)
@pytest.mark.parametrize(
    "bn1,bn2",
    [
        (
            bt.logic.bn.BooleanNetwork({"A": "A"}),
            bt.logic.bn.BooleanNetwork(),
        ),
        (
            bt.logic.bn.BooleanNetwork({"A": "A"}),
            bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"}),
        ),
    ],
)
def test_different_component_sets_are_rejected(comparison, bn1, bn2):

    with pytest.raises(ValueError, match="must have the same component set"):
        comparison(bn1, bn2)


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_empty_networks_are_vacuously_identical(metric):

    bn1 = bt.logic.bn.BooleanNetwork()
    bn2 = bt.logic.bn.BooleanNetwork()

    assert bt.logic.bn.distance(bn1, bn2, metric=metric) == 0.0
    assert bt.logic.bn.similarity(bn1, bn2, metric=metric) == 1.0


@pytest.mark.parametrize(
    "rule1,rule2",
    [
        ("B & C", "C & B"),
        ("B", "B & (B | C)"),
        ("B | False", "B"),
        ("(B & C) | (~B & C)", "C"),
        ("B & (C | ~C)", "B"),
    ],
)
@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_equivalent_rule_syntax_has_zero_distance(rule1, rule2, metric):

    bn1 = bt.logic.bn.BooleanNetwork({"A": rule1, "B": "B", "C": "C"})
    bn2 = bt.logic.bn.BooleanNetwork({"A": rule2, "B": "B", "C": "C"})

    assert bt.logic.bn.distance(bn1, bn2, metric=metric) == 0.0
    assert bt.logic.bn.distance(bn2, bn1, metric=metric) == 0.0
    assert bt.logic.bn.similarity(bn1, bn2, metric=metric) == 1.0
    assert bt.logic.bn.similarity(bn2, bn1, metric=metric) == 1.0


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_equivalent_constant_rules_have_zero_distance(metric):

    bn1 = bt.logic.bn.BooleanNetwork({"A": False, "B": "B"})
    bn2 = bt.logic.bn.BooleanNetwork({"A": "B & ~B", "B": "B"})

    assert bt.logic.bn.distance(bn1, bn2, metric=metric) == 0.0
    assert bt.logic.bn.similarity(bn1, bn2, metric=metric) == 1.0


def test_hamming_distance_rebases_equivalent_noncanonical_constants():

    ba1 = BooleanAlgebra()
    ba2 = BooleanAlgebra()
    symbol1 = ba1.Symbol("A")
    symbol2 = ba2.Symbol("A")
    bn1 = bt.logic.bn.BooleanNetwork._from_expressions(
        {"A": symbol1},
        ba=ba1,
        check=True,
    )
    bn2 = bt.logic.bn.BooleanNetwork._from_expressions(
        {"A": ba2.OR(symbol2, _FALSE())},
        ba=ba2,
        check=True,
    )

    assert isinstance(bn1["A"], Expression)
    assert bt.logic.bn.distance(bn1, bn2) == 0.0
    assert bt.logic.bn.similarity(bn1, bn2) == 1.0


@pytest.mark.parametrize(
    "changed,expected_distance",
    [
        (0, 0.0),
        (1, 0.25),
        (2, 0.5),
        (4, 1.0),
    ],
)
def test_equivalence_distance_counts_changed_update_functions(
    changed,
    expected_distance,
):

    components = ("A", "B", "C", "D")
    bn1 = bt.logic.bn.BooleanNetwork({component: component for component in components})
    bn2 = bt.logic.bn.BooleanNetwork(
        {
            component: f"~{component}" if index < changed else component
            for index, component in enumerate(components)
        }
    )

    assert bt.logic.bn.distance(bn1, bn2, metric="equivalence") == expected_distance
    assert (
        bt.logic.bn.similarity(bn1, bn2, metric="equivalence")
        == 1.0 - expected_distance
    )


def test_hamming_distance_for_conjunction_and_literal_is_one_quarter():

    components = ("A", "B", "C")
    bn1 = bt.logic.bn.BooleanNetwork({component: "B & C" for component in components})
    bn2 = bt.logic.bn.BooleanNetwork({component: "B" for component in components})

    assert bt.logic.bn.distance(bn1, bn2) == 1 / 4
    assert bt.logic.bn.similarity(bn1, bn2) == 3 / 4


def test_hamming_distance_between_literal_and_negation_is_one():

    bn1 = bt.logic.bn.BooleanNetwork({"B": "B"})
    bn2 = bt.logic.bn.BooleanNetwork({"B": "~B"})

    assert bt.logic.bn.distance(bn1, bn2) == 1.0
    assert bt.logic.bn.similarity(bn1, bn2) == 0.0


def test_hamming_distance_between_false_and_true_is_one():

    bn1 = bt.logic.bn.BooleanNetwork({"A": False})
    bn2 = bt.logic.bn.BooleanNetwork({"A": True})

    assert bt.logic.bn.distance(bn1, bn2) == 1.0
    assert bt.logic.bn.similarity(bn1, bn2) == 0.0


def test_hamming_distance_averages_local_distances_uniformly():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & C",
            "B": "B",
            "C": "C",
        }
    )
    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "~B",
            "C": "C",
        }
    )

    assert bt.logic.bn.distance(bn1, bn2) == 5 / 12
    assert bt.logic.bn.similarity(bn1, bn2) == 7 / 12


def test_hamming_distance_handles_reduced_constant_subexpressions():

    bn1 = bt.logic.bn.BooleanNetwork({"A": "B | C", "B": "B", "C": "C", "D": "D"})
    bn2 = bt.logic.bn.BooleanNetwork({"A": "B | C | D", "B": "B", "C": "C", "D": "D"})

    assert bt.logic.bn.distance(bn1, bn2) == 1 / 32
    assert bt.logic.bn.similarity(bn1, bn2) == 31 / 32


def test_hamming_distance_is_invariant_to_ignored_variables():

    small_components = ("B", "C")
    large_components = ("B", "C", "D", "E")
    small1 = bt.logic.bn.BooleanNetwork(
        {component: "B & C" for component in small_components}
    )
    small2 = bt.logic.bn.BooleanNetwork(
        {component: "B" for component in small_components}
    )
    large1 = bt.logic.bn.BooleanNetwork(
        {component: "B & C" for component in large_components}
    )
    large2 = bt.logic.bn.BooleanNetwork(
        {component: "B" for component in large_components}
    )

    assert bt.logic.bn.distance(small1, small2) == 1 / 4
    assert bt.logic.bn.distance(large1, large2) == 1 / 4
    assert bt.logic.bn.similarity(small1, small2) == 3 / 4
    assert bt.logic.bn.similarity(large1, large2) == 3 / 4


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_distance_and_similarity_are_direct_exact_ratios(metric):

    bn1 = bt.logic.bn.BooleanNetwork({"A": 0, "B": 0, "C": 0})
    bn2 = bt.logic.bn.BooleanNetwork({"A": 1, "B": 0, "C": 0})

    distance = bt.logic.bn.distance(bn1, bn2, metric=metric)
    similarity = bt.logic.bn.similarity(bn1, bn2, metric=metric)

    assert distance == 1 / 3
    assert similarity == 2 / 3
    assert distance + similarity == pytest.approx(1.0)


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_comparisons_are_symmetric_normalized_and_reflexive(metric):

    bn1 = bt.logic.bn.BooleanNetwork({"A": "B & C", "B": "B", "C": 0})
    bn2 = bt.logic.bn.BooleanNetwork({"A": "B", "B": "~B", "C": 1})

    distance12 = bt.logic.bn.distance(bn1, bn2, metric=metric)
    distance21 = bt.logic.bn.distance(bn2, bn1, metric=metric)
    similarity12 = bt.logic.bn.similarity(bn1, bn2, metric=metric)
    similarity21 = bt.logic.bn.similarity(bn2, bn1, metric=metric)

    assert 0.0 <= distance12 <= 1.0
    assert 0.0 <= similarity12 <= 1.0
    assert distance12 == distance21
    assert similarity12 == similarity21
    assert distance12 + similarity12 == pytest.approx(1.0)
    assert bt.logic.bn.distance(bn1, bn1, metric=metric) == 0.0
    assert bt.logic.bn.similarity(bn1, bn1, metric=metric) == 1.0


def test_equivalence_and_hamming_measure_different_rule_changes():

    components = ("A", "B", "C")
    bn1 = bt.logic.bn.BooleanNetwork({component: "B & C" for component in components})
    bn2 = bt.logic.bn.BooleanNetwork({component: "B" for component in components})

    assert bt.logic.bn.distance(bn1, bn2, metric="equivalence") == 1.0
    assert bt.logic.bn.distance(bn1, bn2, metric="hamming") == 1 / 4


def test_semantic_distance_does_not_compare_influence_graphs():

    bn1 = bt.logic.bn.BooleanNetwork({"A": "A & B", "B": "B"})
    bn2 = bt.logic.bn.BooleanNetwork({"A": "A | B", "B": "B"})

    assert (
        bt.logic.ig.distance(
            bn1.to_influence_graph(),
            bn2.to_influence_graph(),
        )
        == 0.0
    )
    assert bt.logic.bn.distance(bn1, bn2) == 1 / 4


@pytest.mark.parametrize(
    "comparison",
    [bt.logic.bn.distance, bt.logic.bn.similarity],
)
def test_unsupported_metric_is_rejected(comparison):

    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    with pytest.raises(
        ValueError,
        match="unsupported Boolean-network comparison metric",
    ):
        comparison(bn, bn, metric=cast(Any, "cosine"))


def test_distance_does_not_delegate_to_similarity(monkeypatch):

    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    def fail(*args, **kwargs):
        raise AssertionError("distance delegated to similarity")

    monkeypatch.setattr(_distances, "similarity", fail)

    assert _distances.distance(bn, bn) == 0.0


def test_similarity_does_not_delegate_to_distance(monkeypatch):

    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    def fail(*args, **kwargs):
        raise AssertionError("similarity delegated to distance")

    monkeypatch.setattr(_distances, "distance", fail)

    assert _distances.similarity(bn, bn) == 1.0
