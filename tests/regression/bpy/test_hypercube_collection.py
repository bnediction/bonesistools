#!/usr/bin/env python

from types import MappingProxyType
from typing import Any, cast

import pytest

import bonesistools as bt


def test_hypercube_collection_initialization_and_set_behavior():

    hcs = bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0},
            {"A": 0, "B": "*"},
            {"A": 1},
        ]
    )

    assert len(hcs) == 2
    assert {"A": 0} in hcs
    assert {"A": 0, "B": "*"} in hcs
    assert {"A": 1} in hcs
    assert {"A": 2} not in hcs
    assert object() not in hcs
    assert hcs.components == frozenset({"A"})

    equivalent = bt.bpy.ba.HypercubeCollection([{"A": 0}, {"A": 1}])
    assert hcs == equivalent
    assert hcs.__eq__(object()) is NotImplemented
    assert hcs.__eq__([{"A": 2}]) is NotImplemented


def test_hypercube_collection_add_discard_and_copy():

    hcs = bt.bpy.ba.HypercubeCollection()

    hcs.add({"A": 0})
    hcs.add(bt.bpy.ba.Hypercube({"A": 1}))
    hcs.add(MappingProxyType({"B": "*"}))

    assert len(hcs) == 3

    copied = hcs.copy()

    assert copied == hcs
    assert copied is not hcs

    hcs.discard({"A": 0})
    hcs.discard({"A": 3})
    hcs.discard(cast(Any, object()))

    assert {"A": 0} not in hcs
    assert {"A": 1} in hcs
    assert {"B": "*"} in hcs
    assert {"A": 1} in copied
    assert {"B": "*"} in copied


def test_hypercube_collection_filtering():

    hcs = bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 0},
            {"A": 1},
            {"A": "*"},
        ]
    )

    assert hcs.fully_specified() == bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
        ]
    )

    assert hcs.smaller_than({"A": "*"}) == bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 0},
            {"A": 1},
            {"A": "*"},
        ]
    )

    larger = hcs.larger_than({"A": 0, "B": 1})

    assert {"A": 0, "B": 1} in larger
    assert {"A": 0} in larger
    assert {"A": "*"} in larger
    assert {"A": 1} not in larger
    assert len(larger) == 3


def test_hypercube_collection_minimal_and_maximal():

    hcs = bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 0},
            {"A": 1},
            {"A": "*"},
        ]
    )

    assert hcs.minimal() == bt.bpy.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 1},
        ]
    )

    assert hcs.maximal() == bt.bpy.ba.HypercubeCollection(
        [
            {"A": "*"},
        ]
    )


def test_hypercube_collection_rejects_invalid_inputs():

    with pytest.raises(ValueError):
        bt.bpy.ba.HypercubeCollection([{"A": 2}])

    hcs = bt.bpy.ba.HypercubeCollection()

    with pytest.raises(TypeError):
        hcs.add(cast(Any, object()))

    with pytest.raises(ValueError):
        hcs.add({"A": 2})

    with pytest.raises(TypeError):
        hcs.smaller_than(cast(Any, object()))

    with pytest.raises(ValueError):
        hcs.larger_than({"A": 2})
