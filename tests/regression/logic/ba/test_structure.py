#!/usr/bin/env python

from typing import Any, cast

import pytest
from boolean import BooleanAlgebra

import bonesistools as bt
from bonesistools.logic.boolean_algebra import _structure


def test_expressions_equivalent_asp_detects_equivalence():
    ba = BooleanAlgebra()

    assert bt.logic.ba.expressions_equivalent(
        ba.parse("(B & C) | (~B & D) | (C & D)"),
        ba.parse("(B & C) | (~B & D)"),
        method="asp",
        ba=ba,
    )


def test_expressions_equivalent_asp_detects_nonequivalence():
    ba = BooleanAlgebra()

    assert not bt.logic.ba.expressions_equivalent(
        ba.parse("B & C"),
        ba.parse("B | C"),
        method="asp",
        ba=ba,
    )


def test_dnf_implicants_conjunction():
    assert bt.logic.ba.dnf_implicants("x & ~y") == (
        bt.logic.ba.Hypercube({"x": 1, "y": 0}),
    )


def test_dnf_implicants_disjunction():
    assert bt.logic.ba.dnf_implicants("x | ~y") == (
        bt.logic.ba.Hypercube({"x": 1}),
        bt.logic.ba.Hypercube({"y": 0}),
    )


def test_dnf_implicants_constants():
    assert bt.logic.ba.dnf_implicants(1) == (bt.logic.ba.Hypercube({}),)
    assert bt.logic.ba.dnf_implicants(0) == ()
    assert bt.logic.ba.dnf_implicants(0, value=0) == (bt.logic.ba.Hypercube({}),)


def test_dnf_implicants_distributes_conjunction_over_disjunction():
    assert bt.logic.ba.dnf_implicants("A & (B | C)") == (
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
        bt.logic.ba.Hypercube({"A": 1, "C": 1}),
    )


def test_dnf_implicants_for_false_value():
    assert bt.logic.ba.dnf_implicants("A & !B", value=0) == (
        bt.logic.ba.Hypercube({"A": 0}),
        bt.logic.ba.Hypercube({"B": 1}),
    )


def test_dnf_implicants_do_not_absorb_redundant_clauses():
    rule = "(A & B) | (A & !B)"

    assert bt.logic.ba.dnf_implicants(rule) == (
        bt.logic.ba.Hypercube({"A": 1, "B": 0}),
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
    )
    assert bt.logic.ba.prime_implicants(rule) == (bt.logic.ba.Hypercube({"A": 1}),)


def test_dnf_implicants_rejects_invalid_value():
    with pytest.raises(ValueError, match="expected 0 or 1"):
        bt.logic.ba.dnf_implicants("A", value=cast(Any, 2))


def test_prime_implicants_disjunction():
    assert bt.logic.ba.prime_implicants("B | C") == (
        bt.logic.ba.Hypercube({"B": 1}),
        bt.logic.ba.Hypercube({"C": 1}),
    )


def test_prime_implicants_conjunction_with_negation():
    assert bt.logic.ba.prime_implicants("B & !C") == (
        bt.logic.ba.Hypercube({"B": 1, "C": 0}),
    )


def test_prime_implicants_for_false_value():
    assert bt.logic.ba.prime_implicants("B & !C", value=0) == (
        bt.logic.ba.Hypercube({"B": 0}),
        bt.logic.ba.Hypercube({"C": 1}),
    )


def test_prime_implicants_minimize_redundant_conditions():
    assert bt.logic.ba.prime_implicants("(A & B) | (A & !B)") == (
        bt.logic.ba.Hypercube({"A": 1}),
    )


def test_prime_implicants_constants():
    assert bt.logic.ba.prime_implicants(1) == (bt.logic.ba.Hypercube({}),)
    assert bt.logic.ba.prime_implicants(0) == ()
    assert bt.logic.ba.prime_implicants(0, value=0) == (bt.logic.ba.Hypercube({}),)


def test_prime_implicants_evaluates_unsimplified_constant_negation():
    ba = BooleanAlgebra()

    assert bt.logic.ba.prime_implicants(ba.NOT(ba.FALSE), ba=ba) == (
        bt.logic.ba.Hypercube({}),
    )


def test_prime_implicants_rejects_invalid_value():
    with pytest.raises(ValueError, match="expected 0 or 1"):
        bt.logic.ba.prime_implicants("A", value=cast(Any, 2))


def test_prime_implicants_uses_truth_table_for_small_rules(monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("small rules must not use ASP")

    monkeypatch.setattr(
        _structure,
        "_compute_prime_cubes_with_clingo",
        fail_if_called,
    )

    assert bt.logic.ba.prime_implicants("A | B") == (
        bt.logic.ba.Hypercube({"A": 1}),
        bt.logic.ba.Hypercube({"B": 1}),
    )


def test_prime_implicants_uses_asp_for_large_rules(monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("large rules must not enumerate minterms")

    monkeypatch.setattr(_structure, "_matching_minterms", fail_if_called)

    assert bt.logic.ba.prime_implicants(
        "(A & B) | (!A & C) | (D & E)",
    ) == (
        bt.logic.ba.Hypercube({"A": 0, "C": 1}),
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
        bt.logic.ba.Hypercube({"B": 1, "C": 1}),
        bt.logic.ba.Hypercube({"D": 1, "E": 1}),
    )
