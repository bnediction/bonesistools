#!/usr/bin/env python

import json
from types import MappingProxyType
from typing import Any, cast

import pytest

import bonesistools as bt


def test_hypercube_mapping_interface():

    hc = bt.logic.ba.Hypercube({"A": 0, "B": "*"})
    hc["C"] = True

    assert dict(hc) == {
        "A": bt.logic.ba.PartialBoolean(0),
        "B": bt.logic.ba.PartialBoolean("*"),
        "C": bt.logic.ba.PartialBoolean(1),
    }

    assert hc["A"] == 0
    assert hc["B"] == "*"
    assert hc["C"] == 1
    assert len(hc) == 3

    assert hc.components == frozenset({"A", "B", "C"})
    assert hc.is_fully_specified is False

    hc["B"] = 1

    assert hc.is_fully_specified is True


def test_hypercube_copy_drop_and_update():

    hc = bt.logic.ba.Hypercube({"A": 0, "B": "*"})
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


def test_hypercube_rename_and_relabel():
    hc = bt.logic.ba.Hypercube({"Trp53": 1, "Myc": 0, "free": "*"})

    hc.rename("Trp53", "TP53")
    assert hc == {"TP53": 1, "Myc": 0, "free": "*"}

    hc.relabel({"Myc": "MYC", "missing": "ignored"})
    assert hc == {"TP53": 1, "MYC": 0, "free": "*"}

    hc.relabel({"TP53": "MYC_old", "MYC": "TP53"})
    assert hc == {"MYC_old": 1, "TP53": 0, "free": "*"}

    hc.rename("free", "free")
    assert hc == {"MYC_old": 1, "TP53": 0, "free": "*"}


def test_hypercube_relabel_rejects_invalid_inputs_and_component_merges():
    hc = bt.logic.ba.Hypercube({"A": 0, "B": 1, "C": "*"})

    with pytest.raises(KeyError, match="component 'missing' not found"):
        hc.rename("missing", "D")

    with pytest.raises(ValueError, match="merge hypercube components"):
        hc.rename("A", "B")

    with pytest.raises(ValueError, match="merge hypercube components"):
        hc.relabel({"A": "D", "B": "D"})

    with pytest.raises(TypeError, match="unsupported argument type for 'old'"):
        hc.rename(cast(Any, 1), "D")

    with pytest.raises(TypeError, match="unsupported argument type for 'new'"):
        hc.rename("A", cast(Any, 1))

    with pytest.raises(TypeError, match="unsupported argument type for 'mapping'"):
        hc.relabel(cast(Any, object()))

    with pytest.raises(TypeError, match="unsupported mapping key type"):
        hc.relabel(cast(Any, {1: "D"}))

    with pytest.raises(TypeError, match="unsupported mapping value type"):
        hc.relabel(cast(Any, {"A": 1}))


def test_hypercube_comparisons():

    fixed = bt.logic.ba.Hypercube({"A": 0, "B": 1})
    partial = bt.logic.ba.Hypercube({"A": 0})
    explicit_partial = bt.logic.ba.Hypercube({"A": 0, "B": "*"})
    readonly_partial = MappingProxyType({"A": 0, "B": "*"})

    assert partial == explicit_partial
    assert partial == readonly_partial
    assert partial.contains(fixed)
    assert fixed < partial
    assert fixed <= partial
    assert partial > fixed
    assert partial >= fixed
    assert fixed.is_smaller_than(partial)
    assert partial.is_larger_than(fixed)
    assert fixed.is_strictly_smaller_than(partial)
    assert partial.is_strictly_larger_than(fixed)

    assert not partial < fixed
    assert not fixed > partial

    assert bt.logic.ba.Hypercube({"A": 0}) != bt.logic.ba.Hypercube({"A": 1})
    assert repr(fixed) == "Hypercube(A=0, B=1)"
    assert repr(explicit_partial) == "Hypercube(A=0, B=*)"
    assert repr(bt.logic.ba.Hypercube({"A-B": 1})) == "Hypercube(A-B=1)"
    assert repr(bt.logic.ba.Hypercube({"A,B": 1})) == "Hypercube('A,B'=1)"


def test_hypercube_identical_and_different():

    hc1 = bt.logic.ba.Hypercube({"A": 0, "B": 1, "C": "*"})
    hc2 = bt.logic.ba.Hypercube({"A": 0, "B": "*", "D": 1})

    assert hc1.identical(hc2) == {"A", "C"}
    assert hc1.different(hc2) == {"B", "D"}


def test_hypercube_invalid_values():

    with pytest.raises(ValueError):
        bt.logic.ba.Hypercube({"A": 2})

    with pytest.raises(ValueError):
        bt.logic.ba.Hypercube({"A": cast(Any, {"nested": 1})})

    hc = bt.logic.ba.Hypercube({"A": 0})

    with pytest.raises(ValueError):
        hc["B"] = 3

    with pytest.raises(TypeError):
        hc.contains(cast(Any, object()))


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

    hypercubes = bt.logic.io.read_hypercubes(file)

    assert set(hypercubes) == {"hc1", "hc2"}
    assert hypercubes["hc1"] == {"A": 0, "B": "*"}
    assert hypercubes["hc2"] == {"A": 1, "B": "*"}


def test_read_hypercubes_from_csv_columns_and_rows(tmp_path):
    file = tmp_path / "hypercubes.csv"
    file.write_text(
        ",hc1,hc2\n" "A,0,1\n" "B,,0\n",
    )

    by_columns = bt.logic.io.read_hypercubes(file, orientation="columns")
    by_rows = bt.logic.io.read_hypercubes(file, orientation="rows")

    assert by_columns["hc1"] == {"A": 0, "B": "*"}
    assert by_columns["hc2"] == {"A": 1, "B": 0}
    assert by_rows["A"] == {"hc1": 0, "hc2": 1}
    assert by_rows["B"] == {"hc1": "*", "hc2": 0}


def test_read_hypercubes_from_tsv_and_rejects_invalid_inputs(tmp_path):
    tsv_file = tmp_path / "hypercubes.tsv"
    tsv_file.write_text(
        "\thc1\n" "A\t1\n",
    )

    hypercubes = bt.logic.io.read_hypercubes(tsv_file)
    assert hypercubes["hc1"] == {"A": 1}

    with pytest.warns(FutureWarning):
        with pytest.raises(ValueError):
            bt.logic.io.read_hypercubes(tsv_file, axis=cast(Any, "diagonal"))

    unsupported_file = tmp_path / "hypercubes.txt"
    unsupported_file.write_text("A=1")

    with pytest.raises(ValueError, match="unsupported input format"):
        bt.logic.io.read_hypercubes(unsupported_file)


def test_deprecated_read_hypercubes_routes_to_io(tmp_path):
    file = tmp_path / "hypercubes.csv"
    file.write_text(",hc1\nA,1\n")

    with pytest.warns(FutureWarning, match="bt.logic.ba.read_hypercubes"):
        hypercubes = getattr(bt.logic.ba, "read_hypercubes")(file)

    assert hypercubes == bt.logic.io.read_hypercubes(file)


def test_hypercube_collection_initialization_and_set_behavior():

    hcs = bt.logic.ba.HypercubeCollection(
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

    equivalent = bt.logic.ba.HypercubeCollection([{"A": 0}, {"A": 1}])
    assert hcs == equivalent
    assert hcs.__eq__(object()) is NotImplemented
    assert hcs.__eq__([{"A": 2}]) is NotImplemented


def test_hypercube_collection_add_discard_and_copy():

    hcs = bt.logic.ba.HypercubeCollection()

    hcs.add({"A": 0})
    hcs.add(bt.logic.ba.Hypercube({"A": 1}))
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

    hcs = bt.logic.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 0},
            {"A": 1},
            {"A": "*"},
        ]
    )

    assert hcs.fully_specified() == bt.logic.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
        ]
    )

    assert hcs.smaller_than({"A": "*"}) == bt.logic.ba.HypercubeCollection(
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

    hcs = bt.logic.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 0},
            {"A": 1},
            {"A": "*"},
        ]
    )

    assert hcs.minimal() == bt.logic.ba.HypercubeCollection(
        [
            {"A": 0, "B": 1},
            {"A": 1},
        ]
    )

    assert hcs.maximal() == bt.logic.ba.HypercubeCollection(
        [
            {"A": "*"},
        ]
    )


def test_hypercube_collection_rejects_invalid_inputs():

    with pytest.raises(ValueError):
        bt.logic.ba.HypercubeCollection([{"A": 2}])

    hcs = bt.logic.ba.HypercubeCollection()

    with pytest.raises(TypeError):
        hcs.add(cast(Any, object()))

    with pytest.raises(ValueError):
        hcs.add({"A": 2})

    with pytest.raises(TypeError):
        hcs.smaller_than(cast(Any, object()))

    with pytest.raises(ValueError):
        hcs.larger_than({"A": 2})
