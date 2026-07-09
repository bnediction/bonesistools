#!/usr/bin/env python

import bonesistools as bt
from bonesistools.boolpy.boolean_algebra import _typing


def test_typing_helpers_are_not_exposed_publicly():
    for name in [
        "BooleanRule",
        "ConfigurationLike",
        "HypercubeLike",
        "Implicants",
        "KleeneValueLike",
        "PartialBooleanLike",
        "is_boolean_expression_available",
        "is_boolean_expression_like",
        "is_boolean_rule_like",
        "is_configuration_like",
        "is_hypercube_like",
        "is_kleene_value_like",
        "is_partial_boolean_like",
    ]:
        assert name not in dir(bt.bpy.ba)
        assert not hasattr(bt.bpy.ba, name)


def test_bpy_namespace_exposes_short_aliases_only():
    assert {"ba", "bn", "ig", "io"} <= set(dir(bt.bpy))
    for name in [
        "boolean_algebra",
        "boolean_network",
        "influence_graph",
        "input_output",
        "plotting",
    ]:
        assert name not in dir(bt.bpy)
        assert not hasattr(bt.bpy, name)

    assert "pl" not in dir(bt.bpy)
    assert not hasattr(bt.bpy, "pl")


def test_configuration_like_is_concrete_and_distinct_from_hypercube_like():
    assert _typing.is_configuration_like({"A": 0, "B": True})
    assert _typing.is_configuration_like({"A": bt.bpy.ba.PartialBoolean(1)})

    assert not _typing.is_configuration_like({"A": "*"})
    assert not _typing.is_configuration_like({"A": bt.bpy.ba.PartialBoolean("*")})
    assert not _typing.is_configuration_like({"A": 2})
    assert not _typing.is_configuration_like({1: 0})
    assert not _typing.is_configuration_like([("A", 0)])

    assert _typing.is_hypercube_like({"A": "*"})
