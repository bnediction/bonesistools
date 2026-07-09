#!/usr/bin/env python

import pytest

import bonesistools as bt


def test_configuration_set_adds_configurations_and_iterates_complete_states():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"])
    configurations.add({"A": 0, "B": 0})
    configurations.add({"A": 0, "B": 1})

    assert len(configurations) == 2
    assert configurations.count() == 2
    assert list(configurations) == [
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
    ]
    assert configurations.enumerate() == (
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
    )


def test_configuration_set_adds_subspaces_without_public_hypercube_semantics():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"])
    configurations.add({"A": "*", "B": 0})

    assert len(configurations) == 2
    assert {"A": 0, "B": 0} in configurations
    assert {"A": 1, "B": 0} in configurations
    assert {"A": 0, "B": 1} not in configurations
    assert {"B": 0} in configurations


def test_configuration_set_replaces_redundant_specific_subspaces():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"])
    configurations.add({"A": 0, "B": 0})
    configurations.add({"A": "*", "B": 0})

    assert len(configurations) == 2
    assert len(configurations._hypercubes) == 1
    assert configurations.enumerate() == (
        {"A": 0, "B": 0},
        {"A": 1, "B": 0},
    )


def test_configuration_set_keeps_exact_union_for_partially_overlapping_subspaces():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"])
    configurations.add({"A": 0})
    configurations.add({"B": 0})

    assert len(configurations) == 3
    assert set(tuple(state.items()) for state in configurations) == {
        (("A", 0), ("B", 0)),
        (("A", 0), ("B", 1)),
        (("A", 1), ("B", 0)),
    }
    assert {"A": 1, "B": 1} not in configurations


def test_configuration_set_compresses_adjacent_subspaces():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"])
    configurations.add({"A": 0, "B": 0})
    configurations.add({"A": 1, "B": 0})

    assert len(configurations._hypercubes) == 2

    configurations.compress()

    assert len(configurations) == 2
    assert len(configurations._hypercubes) == 1
    assert configurations.enumerate() == (
        {"A": 0, "B": 0},
        {"A": 1, "B": 0},
    )


def test_configuration_set_equality_is_semantic_not_representational():

    compact = bt.bpy.ba.ConfigurationSet(["A", "B"], [{"B": 0}])

    expanded = bt.bpy.ba.ConfigurationSet(["A", "B"])
    expanded.add({"A": 0, "B": 0})
    expanded.add({"A": 1, "B": 0})

    reordered = bt.bpy.ba.ConfigurationSet(["B", "A"], [{"B": 0}])
    different = bt.bpy.ba.ConfigurationSet(["A", "B"], [{"A": 0}])
    different_components = bt.bpy.ba.ConfigurationSet(["A", "C"], [{"A": 0}])

    assert compact == expanded
    assert expanded == compact
    assert compact == reordered
    assert compact != different
    assert compact != different_components
    assert compact.__eq__(object()) is NotImplemented


def test_configuration_set_samples_reproducibly():

    configurations = bt.bpy.ba.ConfigurationSet(["A", "B"], [{"A": "*"}])

    assert configurations.sample(seed=0) == {"A": 0, "B": 0}
    assert configurations.sample(3, seed=0) == (
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
        {"A": 0, "B": 1},
    )


def test_configuration_set_rejects_invalid_inputs():

    with pytest.raises(ValueError, match="duplicated"):
        bt.bpy.ba.ConfigurationSet(["A", "A"])

    configurations = bt.bpy.ba.ConfigurationSet(["A"])

    with pytest.raises(ValueError, match="unknown components"):
        configurations.add({"B": 0})

    with pytest.raises(ValueError):
        configurations.add({"A": 2})

    with pytest.raises(ValueError, match="empty"):
        configurations.sample()
