#!/usr/bin/env python

from typing import Any, cast

import pytest
from boolean import BooleanAlgebra

import bonesistools as bt
from bonesistools.boolpy.boolean_algebra import _structure


def test_expressions_equivalent_asp_detects_equivalence():
    ba = BooleanAlgebra()

    assert bt.bpy.ba.expressions_equivalent(
        ba.parse("(B & C) | (~B & D) | (C & D)"),
        ba.parse("(B & C) | (~B & D)"),
        method="asp",
        ba=ba,
    )


def test_expressions_equivalent_asp_detects_nonequivalence():
    ba = BooleanAlgebra()

    assert not bt.bpy.ba.expressions_equivalent(
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
        bt.bpy.ba.expressions_equivalent(
            ba.parse("A"),
            ba.parse("A"),
            method="asp",
            ba=ba,
        )


def test_dnf_implicants_conjunction():
    assert bt.bpy.ba.dnf_implicants("x & ~y") == (
        bt.bpy.ba.Hypercube({"x": 1, "y": 0}),
    )


def test_dnf_implicants_disjunction():
    assert bt.bpy.ba.dnf_implicants("x | ~y") == (
        bt.bpy.ba.Hypercube({"x": 1}),
        bt.bpy.ba.Hypercube({"y": 0}),
    )


def test_dnf_implicants_constants():
    assert bt.bpy.ba.dnf_implicants(1) == (bt.bpy.ba.Hypercube({}),)
    assert bt.bpy.ba.dnf_implicants(0) == ()
    assert bt.bpy.ba.dnf_implicants(0, value=0) == (bt.bpy.ba.Hypercube({}),)


def test_dnf_implicants_distributes_conjunction_over_disjunction():
    assert bt.bpy.ba.dnf_implicants("A & (B | C)") == (
        bt.bpy.ba.Hypercube({"A": 1, "B": 1}),
        bt.bpy.ba.Hypercube({"A": 1, "C": 1}),
    )


def test_dnf_implicants_for_false_value():
    assert bt.bpy.ba.dnf_implicants("A & !B", value=0) == (
        bt.bpy.ba.Hypercube({"A": 0}),
        bt.bpy.ba.Hypercube({"B": 1}),
    )


def test_dnf_implicants_do_not_absorb_redundant_clauses():
    rule = "(A & B) | (A & !B)"

    assert bt.bpy.ba.dnf_implicants(rule) == (
        bt.bpy.ba.Hypercube({"A": 1, "B": 0}),
        bt.bpy.ba.Hypercube({"A": 1, "B": 1}),
    )
    assert bt.bpy.ba.prime_implicants(rule) == (bt.bpy.ba.Hypercube({"A": 1}),)


def test_dnf_implicants_rejects_invalid_value():
    with pytest.raises(ValueError, match="expected 0 or 1"):
        bt.bpy.ba.dnf_implicants("A", value=cast(Any, 2))


def test_prime_implicants_disjunction():
    assert bt.bpy.ba.prime_implicants("B | C") == (
        bt.bpy.ba.Hypercube({"B": 1}),
        bt.bpy.ba.Hypercube({"C": 1}),
    )


def test_prime_implicants_conjunction_with_negation():
    assert bt.bpy.ba.prime_implicants("B & !C") == (
        bt.bpy.ba.Hypercube({"B": 1, "C": 0}),
    )


def test_prime_implicants_for_false_value():
    assert bt.bpy.ba.prime_implicants("B & !C", value=0) == (
        bt.bpy.ba.Hypercube({"B": 0}),
        bt.bpy.ba.Hypercube({"C": 1}),
    )


def test_prime_implicants_minimize_redundant_conditions():
    assert bt.bpy.ba.prime_implicants("(A & B) | (A & !B)") == (
        bt.bpy.ba.Hypercube({"A": 1}),
    )


def test_prime_implicants_constants():
    assert bt.bpy.ba.prime_implicants(1) == (bt.bpy.ba.Hypercube({}),)
    assert bt.bpy.ba.prime_implicants(0) == ()
    assert bt.bpy.ba.prime_implicants(0, value=0) == (bt.bpy.ba.Hypercube({}),)


def test_prime_implicants_evaluates_unsimplified_constant_negation():
    ba = BooleanAlgebra()

    assert bt.bpy.ba.prime_implicants(ba.NOT(ba.FALSE), ba=ba) == (
        bt.bpy.ba.Hypercube({}),
    )


def test_prime_implicants_rejects_invalid_value():
    with pytest.raises(ValueError, match="expected 0 or 1"):
        bt.bpy.ba.prime_implicants("A", value=cast(Any, 2))


def test_prime_implicants_rejects_invalid_backend():
    with pytest.raises(ValueError, match="expected one of"):
        bt.bpy.ba.prime_implicants("A", backend=cast(Any, "unknown"))


def test_prime_implicants_asp_backend_matches_truth_table():
    rule = "(A & B) | (A & !C) | (!A & B & C)"

    assert bt.bpy.ba.prime_implicants(rule, backend="asp") == (
        bt.bpy.ba.Hypercube({"A": 1, "B": 1}),
        bt.bpy.ba.Hypercube({"A": 1, "C": 0}),
        bt.bpy.ba.Hypercube({"B": 1, "C": 1}),
    )
    assert bt.bpy.ba.prime_implicants(
        rule,
        value=0,
        backend="asp",
    ) == bt.bpy.ba.prime_implicants(rule, value=0)


def test_prime_implicants_asp_backend_reports_missing_clingo(monkeypatch):
    def missing_clingo(name: str) -> Any:
        if name == "clingo":
            raise ImportError("missing clingo")

        return __import__(name)

    monkeypatch.setattr(_structure, "import_module", missing_clingo)

    with pytest.raises(ImportError, match="requires `clingo` to be installed"):
        bt.bpy.ba.prime_implicants("A", backend="asp")
