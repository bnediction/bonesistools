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


def test_expressions_equivalent_asp_reports_missing_clingo(monkeypatch):
    ba = BooleanAlgebra()

    def missing_clingo(name: str) -> Any:
        if name == "clingo":
            raise ImportError("missing clingo")

        return __import__(name)

    monkeypatch.setattr(_structure, "import_module", missing_clingo)

    with pytest.raises(ImportError, match="requires `clingo` to be installed"):
        bt.logic.ba.expressions_equivalent(
            ba.parse("A"),
            ba.parse("A"),
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


def test_prime_implicants_rejects_invalid_backend():
    with pytest.raises(ValueError, match="expected one of"):
        bt.logic.ba.prime_implicants("A", backend=cast(Any, "unknown"))


def test_prime_implicants_asp_backend_matches_truth_table():
    rule = "(A & B) | (A & !C) | (!A & B & C)"

    assert bt.logic.ba.prime_implicants(rule, backend="asp") == (
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
        bt.logic.ba.Hypercube({"A": 1, "C": 0}),
        bt.logic.ba.Hypercube({"B": 1, "C": 1}),
    )
    assert bt.logic.ba.prime_implicants(
        rule,
        value=0,
        backend="asp",
    ) == bt.logic.ba.prime_implicants(rule, value=0)


def test_prime_implicants_asp_backend_reports_missing_clingo(monkeypatch):
    def missing_clingo(name: str) -> Any:
        if name == "clingo":
            raise ImportError("missing clingo")

        return __import__(name)

    monkeypatch.setattr(_structure, "import_module", missing_clingo)

    with pytest.raises(ImportError, match="requires `clingo` to be installed"):
        bt.logic.ba.prime_implicants("A", backend="asp")
