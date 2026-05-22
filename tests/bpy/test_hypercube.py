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


def test_read_hypercubes_from_json_converts_null_to_free_value(tmp_path):
    file = tmp_path / "hypercubes.json"

    with open(file, "w") as fp:
        json.dump(
            {
                "hc1": {"A": 0, "B": None},
                "hc2": {"A": 1, "B": "*"},
            },
            fp,
        )

    hypercubes = bt.bpy.ba.read_hypercubes(file)

    assert set(hypercubes) == {"hc1", "hc2"}
    assert hypercubes["hc1"] == {"A": 0, "B": "*"}
    assert hypercubes["hc2"] == {"A": 1, "B": "*"}


def test_read_hypercubes_from_csv_columns_and_rows(tmp_path):
    file = tmp_path / "hypercubes.csv"
    file.write_text(
        ",hc1,hc2\n" "A,0,1\n" "B,,0\n",
    )

    by_columns = bt.bpy.ba.read_hypercubes(file, axis="columns")
    by_rows = bt.bpy.ba.read_hypercubes(file, axis="rows")

    assert by_columns["hc1"] == {"A": 0, "B": "*"}
    assert by_columns["hc2"] == {"A": 1, "B": 0}
    assert by_rows["A"] == {"hc1": 0, "hc2": 1}
    assert by_rows["B"] == {"hc1": "*", "hc2": 0}


def test_read_hypercubes_from_tsv_and_rejects_invalid_inputs(tmp_path):
    tsv_file = tmp_path / "hypercubes.tsv"
    tsv_file.write_text(
        "\thc1\n" "A\t1\n",
    )

    hypercubes = bt.bpy.ba.read_hypercubes(tsv_file)
    assert hypercubes["hc1"] == {"A": 1}

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.bpy.ba.read_hypercubes(tsv_file, axis="diagonal")

    unsupported_file = tmp_path / "hypercubes.txt"
    unsupported_file.write_text("A=1")

    with pytest.raises(ValueError, match="unsupported input format"):
        bt.bpy.ba.read_hypercubes(unsupported_file)
