#!/usr/bin/env python

import pytest

import bonesistools as bt


def test_partial_boolean_equality_and_contains():

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pb1 = bt.bpy.ba.PartialBoolean(1)
    pbstar = bt.bpy.ba.PartialBoolean("*")

    assert pb0 == 0
    assert pb1 == 1
    assert pbstar == "*"

    assert pbstar.contains(0)
    assert pbstar.contains(1)
    assert pbstar.contains("*")
    assert pb0 in pbstar
    assert pb1 in pbstar
    assert pbstar in pbstar

    assert pb0.contains(0)
    assert not pb0.contains(1)
    assert not pb0.contains("*")
    assert pbstar not in pb0
    assert pb1 not in pb0
    assert pb0.__eq__(object()) is NotImplemented
    assert pb0.__ne__(object()) is NotImplemented


def test_partial_boolean_properties_and_conversion():

    assert bt.bpy.ba.PartialBoolean(0).is_fixed is True
    assert bt.bpy.ba.PartialBoolean(1).is_fixed is True
    assert bt.bpy.ba.PartialBoolean("*").is_fixed is False

    assert bt.bpy.ba.PartialBoolean("*").is_free is True
    assert bt.bpy.ba.PartialBoolean(0).is_free is False

    assert bool(bt.bpy.ba.PartialBoolean(0)) is False
    assert bool(bt.bpy.ba.PartialBoolean(1)) is True

    with pytest.raises(ValueError):
        bool(bt.bpy.ba.PartialBoolean("*"))

    assert bt.bpy.ba.PartialBoolean(0).as_set() == frozenset({0})
    assert bt.bpy.ba.PartialBoolean(1).as_set() == frozenset({1})
    assert bt.bpy.ba.PartialBoolean("*").as_set() == frozenset({0, 1})
    assert bt.bpy.ba.PartialBoolean("*").possibilities() == frozenset({0, 1})

    kleene = bt.bpy.ba.PartialBoolean("*").to_kleene()

    assert isinstance(kleene, bt.bpy.ba.KleeneValue)
    assert kleene == bt.bpy.ba.KleeneValue("*")
    assert kleene.to_partial_boolean() == bt.bpy.ba.PartialBoolean("*")


def test_partial_boolean_representation_hashing_and_immutability():

    pb = bt.bpy.ba.PartialBoolean("*")

    assert pb.value == "*"
    assert str(pb) == "*"
    assert repr(pb) == "PartialBoolean('*')"

    assert pb in {bt.bpy.ba.PartialBoolean("*")}

    x = bt.bpy.ba.PartialBoolean(0)
    y = x
    x = bt.bpy.ba.PartialBoolean(1)

    assert x == 1
    assert y == 0


def test_partial_boolean_rejects_invalid_values():

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean(3)

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean("A")

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean(None)

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean("*").contains(object())


def test_partial_boolean_order():

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pbstar = bt.bpy.ba.PartialBoolean("*")
    pb1 = bt.bpy.ba.PartialBoolean(1)

    assert pb0 < pbstar
    assert pb1 < pbstar
    assert pbstar > pb0
    assert pbstar > pb1

    assert pb0 <= pbstar
    assert pb1 <= pbstar
    assert pbstar >= pb0
    assert pbstar >= pb1

    assert not (pb0 < pb1)
    assert not (pb1 < pb0)
    assert not (pb0 <= pb1)
    assert not (pb1 <= pb0)
    assert not (pbstar < pb0)
    assert not (pbstar < pb1)


def test_partial_boolean_order_is_not_kleene_truth_order():

    partial_one = bt.bpy.ba.PartialBoolean(1)
    partial_star = bt.bpy.ba.PartialBoolean("*")
    kleene_one = bt.bpy.ba.KleeneValue(1)
    kleene_star = bt.bpy.ba.KleeneValue("*")

    assert partial_one < partial_star
    assert not (partial_star < partial_one)

    assert kleene_star < kleene_one
    assert not (kleene_one < kleene_star)
