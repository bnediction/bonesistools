#!/usr/bin/env python

import pytest

import bonesistools as bt


def test_partial_boolean_ordering():

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pbstar = bt.bpy.ba.PartialBoolean("*")
    pb1 = bt.bpy.ba.PartialBoolean(1)

    assert pb0 == 0
    assert pbstar == "*"
    assert pb1 == 1

    assert pb0 < pbstar
    assert pbstar < pb1

    assert pb0 <= pbstar
    assert pbstar <= pb1

    assert pb1 > pbstar
    assert pbstar > pb0

    assert pb1 >= pbstar
    assert pbstar >= pb0

    assert pb0 != pb1

    assert sorted([pb1, pb0, pbstar]) == [
        pb0,
        pbstar,
        pb1,
    ]


def test_partial_boolean_properties():

    assert bt.bpy.ba.PartialBoolean(0).is_fixed is True
    assert bt.bpy.ba.PartialBoolean(1).is_fixed is True
    assert bt.bpy.ba.PartialBoolean("*").is_fixed is False

    assert bt.bpy.ba.PartialBoolean("*").is_free is True
    assert bt.bpy.ba.PartialBoolean(0).is_free is False
    assert bt.bpy.ba.PartialBoolean(1).is_free is False


def test_partial_boolean_representation_and_value():
    pb = bt.bpy.ba.PartialBoolean("*")

    assert pb.value == "*"
    assert str(pb) == "*"
    assert repr(pb) == "PartialBoolean('*')"


def test_partial_boolean_boolean_inputs():
    assert bt.bpy.ba.PartialBoolean(False) == 0
    assert bt.bpy.ba.PartialBoolean(True) == 1


def test_partial_boolean_bool_conversion():

    assert bool(bt.bpy.ba.PartialBoolean(0)) is False
    assert bool(bt.bpy.ba.PartialBoolean(1)) is True

    with pytest.raises(ValueError):
        bool(bt.bpy.ba.PartialBoolean("*"))


def test_partial_boolean_hashable():
    values = {
        bt.bpy.ba.PartialBoolean(0),
        bt.bpy.ba.PartialBoolean("*"),
        bt.bpy.ba.PartialBoolean(1),
    }

    assert bt.bpy.ba.PartialBoolean(0) in values
    assert bt.bpy.ba.PartialBoolean("*") in values
    assert bt.bpy.ba.PartialBoolean(1) in values


def test_partial_boolean_immutability():

    x = bt.bpy.ba.PartialBoolean(0)
    y = x

    x = bt.bpy.ba.PartialBoolean(1)

    assert x == 1
    assert y == 0


def test_partial_boolean_invalid_value():

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean(3)

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean("A")

    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean(None)


def test_partial_boolean_unsupported_equality():
    assert (bt.bpy.ba.PartialBoolean(0) == object()) is False
    assert (bt.bpy.ba.PartialBoolean(0) != object()) is True


def test_partial_boolean_unsupported_ordering():
    with pytest.raises(ValueError):
        bt.bpy.ba.PartialBoolean(0) < object()
