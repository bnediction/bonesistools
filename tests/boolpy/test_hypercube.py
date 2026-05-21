#!/usr/bin/env python

import json

import bonesistools as bt

import pytest


def test_hypercube_mapping_interface():

    hc = bt.bpy.ba.Hypercube({"A": 0, "B": "*"})
    hc["C"] = True

    assert dict(hc) == {
        "A": bt.bpy.ba.PartialBoolean(0),
        "B": bt.bpy.ba.PartialBoolean("*"),
        "C": bt.bpy.ba.PartialBoolean(1),
    }

    assert hc["A"] == 0
    assert hc["B"] == "*"
    assert hc["C"] == 1

    assert hc.components == frozenset({"A", "B", "C"})
    assert hc.is_fully_specified is False

    hc["B"] = 1

    assert hc.is_fully_specified is True


def test_hypercube_copy_drop_and_update():

    hc = bt.bpy.ba.Hypercube({"A": 0, "B": "*"})
    copied = hc.copy()

    assert copied == hc
    assert copied is not hc

    dropped = hc.drop(["B"])

    assert dropped == {"A": 0}
    assert hc == {"A": 0, "B": "*"}

    assert hc.drop(["B"], inplace=True) is None
    assert hc == {"A": 0}

    hc.update({"C": True, "D": "*"})

    assert hc == {"A": 0, "C": 1, "D": "*"}


def test_hypercube_comparisons():

    fixed = bt.bpy.ba.Hypercube({"A": 0, "B": 1})
    partial = bt.bpy.ba.Hypercube({"A": 0})
    explicit_partial = bt.bpy.ba.Hypercube({"A": 0, "B": "*"})

    assert partial == explicit_partial
    assert fixed < partial
    assert fixed <= partial
    assert partial > fixed
    assert partial >= fixed

    assert not partial < fixed
    assert not fixed > partial

    assert bt.bpy.ba.Hypercube({"A": 0}) != bt.bpy.ba.Hypercube({"A": 1})


def test_hypercube_identical_and_different():

    hc1 = bt.bpy.ba.Hypercube({"A": 0, "B": 1, "C": "*"})
    hc2 = bt.bpy.ba.Hypercube({"A": 0, "B": "*", "D": 1})

    assert hc1.identical(hc2) == {"A", "C"}
    assert hc1.different(hc2) == {"B", "D"}


def test_hypercube_invalid_values():

    with pytest.raises(ValueError):
        bt.bpy.ba.Hypercube({"A": 2})

    with pytest.raises(ValueError):
        bt.bpy.ba.Hypercube({"A": {"nested": 1}})

    hc = bt.bpy.ba.Hypercube({"A": 0})

    with pytest.raises(ValueError):
        hc["B"] = 3


def test_read_hypercube(tmp_path):

    file = tmp_path / "hypercube.json"

    with open(file, "w") as fp:
        json.dump(
            {
                "A": 0,
                "B": 1,
                "C": "*",
            },
            fp,
        )

    hc = bt.bpy.ba.read_hypercube(file)

    assert isinstance(hc, bt.bpy.ba.Hypercube)
    assert hc == {
        "A": 0,
        "B": 1,
        "C": "*",
    }


def test_read_hypercube_rejects_nested_json(tmp_path):

    file = tmp_path / "hypercube.json"

    with open(file, "w") as fp:
        json.dump(
            {
                "A": {
                    "nested": 1,
                }
            },
            fp,
        )

    with pytest.raises(ValueError):
        bt.bpy.ba.read_hypercube(file)
