#!/usr/bin/env python

import bonesistools as bt


def test_configuration_like_is_concrete_and_distinct_from_hypercube_like():
    assert bt.bpy.ba.is_configuration_like({"A": 0, "B": True})
    assert bt.bpy.ba.is_configuration_like({"A": bt.bpy.ba.PartialBoolean(1)})

    assert not bt.bpy.ba.is_configuration_like({"A": "*"})
    assert not bt.bpy.ba.is_configuration_like({"A": bt.bpy.ba.PartialBoolean("*")})
    assert not bt.bpy.ba.is_configuration_like({"A": 2})
    assert not bt.bpy.ba.is_configuration_like({1: 0})
    assert not bt.bpy.ba.is_configuration_like([("A", 0)])

    assert bt.bpy.ba.is_hypercube_like({"A": "*"})
