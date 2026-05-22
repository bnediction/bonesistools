#!/usr/bin/env python

import pytest

from boolean import BooleanAlgebra

import bonesistools as bt


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
