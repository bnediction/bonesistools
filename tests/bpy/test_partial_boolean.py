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

    assert pb0.contains(0)
    assert not pb0.contains(1)
    assert not pb0.contains("*")
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

    assert pb0 < pbstar < pb1
    assert pb0 <= pbstar <= pb1
    assert pb1 > pbstar > pb0
    assert pb1 >= pbstar >= pb0
    assert sorted([pb1, pb0, pbstar]) == [pb0, pbstar, pb1]
