#!/usr/bin/env python

from typing import Any, cast

import pytest

import bonesistools as bt


def test_kleene_differential_uses_truth_order():

    pb0 = bt.bpy.ba.PartialBoolean(0)
    pbstar = bt.bpy.ba.PartialBoolean("*")
    pb1 = bt.bpy.ba.PartialBoolean(1)

    assert bt.bpy.ba.diff(pb0, pb0) == 0
    assert bt.bpy.ba.diff(pbstar, pbstar) == 0
    assert bt.bpy.ba.diff(pb1, pb1) == 0

    assert bt.bpy.ba.diff(pb0, pbstar) == 1
    assert bt.bpy.ba.diff(pbstar, pb1) == 1
    assert bt.bpy.ba.diff(pb0, pb1) == 1

    assert bt.bpy.ba.diff(pbstar, pb0) == -1
    assert bt.bpy.ba.diff(pb1, pbstar) == -1
    assert bt.bpy.ba.diff(pb1, pb0) == -1

    assert bt.bpy.ba.diff(False, 0) == 0
    assert bt.bpy.ba.diff(True, 1) == 0
    assert bt.bpy.ba.diff(float("nan"), pbstar) == 0


def test_kleene_value_order_and_logical_operators():

    false = bt.bpy.ba.KleeneValue(0)
    unknown = bt.bpy.ba.KleeneValue("*")
    true = bt.bpy.ba.KleeneValue(1)

    assert false < unknown < true
    assert false.rank == 0
    assert unknown.rank == 1
    assert true.rank == 2

    assert ~false == true
    assert ~true == false
    assert ~unknown == unknown

    assert true.meet(unknown) == unknown
    assert true.join(unknown) == true
    assert false.meet(true) == false
    assert false.join(true) == true

    assert false & unknown == false
    assert true & unknown == unknown
    assert unknown & unknown == unknown
    assert (true & unknown) == true.meet(unknown)

    assert true | unknown == true
    assert false | unknown == unknown
    assert unknown | unknown == unknown
    assert (true | unknown) == true.join(unknown)

    assert unknown.to_partial_boolean() == bt.bpy.ba.PartialBoolean("*")


def test_kleene_value_representation_hashing_and_unknown_property():
    unknown = bt.bpy.ba.KleeneValue("*")

    assert repr(unknown) == "KleeneValue('*')"
    assert str(unknown) == "*"
    assert hash(unknown) == hash("*")
    assert unknown.is_unknown is True
    assert bt.bpy.ba.KleeneValue(0).is_unknown is False
    assert unknown in {bt.bpy.ba.KleeneValue("*")}


def test_kleene_value_invalid_comparisons_return_not_implemented():
    value = bt.bpy.ba.KleeneValue(0)
    other = object()

    assert value.__eq__(other) is NotImplemented
    assert value.__ne__(other) is NotImplemented
    assert value.__lt__(other) is NotImplemented
    assert value.__le__(other) is NotImplemented
    assert value.__gt__(other) is NotImplemented
    assert value.__ge__(other) is NotImplemented


def test_kleene_value_valid_inequality_methods():
    false = bt.bpy.ba.KleeneValue(0)
    true = bt.bpy.ba.KleeneValue(1)

    assert false.__ne__(true) is True
    assert true.__gt__(false) is True


def test_kleene_meet_and_join_module_wrappers():
    assert bt.bpy.ba.meet(1, "*") == bt.bpy.ba.KleeneValue("*")
    assert bt.bpy.ba.meet(0, 1) == bt.bpy.ba.KleeneValue(0)
    assert bt.bpy.ba.join(0, "*") == bt.bpy.ba.KleeneValue("*")
    assert bt.bpy.ba.join(0, 1) == bt.bpy.ba.KleeneValue(1)


@pytest.mark.parametrize("value", [0, 1, "*"])
def test_kleene_value_cannot_be_converted_to_bool(value):

    with pytest.raises(TypeError):
        bool(bt.bpy.ba.KleeneValue(value))


def test_kleene_diff_rejects_invalid_inputs():

    with pytest.raises(ValueError):
        bt.bpy.ba.diff(2, 1)

    with pytest.raises(ValueError):
        bt.bpy.ba.diff("0", 1)

    with pytest.raises(TypeError):
        bt.bpy.ba.diff(cast(Any, object()), 1)
