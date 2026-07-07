#!/usr/bin/env python

from typing import Any, cast

import pytest

import bonesistools as bt


@pytest.fixture
def bnet_directory(tmp_path):
    models = {
        "model_1.bnet": "A, B\nB, 1\nC, 0\n",
        "model_2.bnet": "A, B&C\nB, 1\nC, 0\n",
        "model_3.bnet": "A, B&!C\nB, 1\nC, 0\n",
    }

    for filename, content in models.items():
        (tmp_path / filename).write_text(content)

    return tmp_path


@pytest.fixture
def bnet_ensemble(bnet_directory):
    return bt.bpy.bn.read_bnet_directory(bnet_directory)


def test_read_bnet_directory_and_ensemble(bnet_ensemble):
    assert len(bnet_ensemble) == 3
    assert bnet_ensemble.components == frozenset({"A", "B", "C"})


def test_read_bnet_directory_validation_and_recursive_loading(tmp_path):
    with pytest.raises(FileNotFoundError, match="directory does not exist"):
        bt.bpy.bn.read_bnet_directory(tmp_path / "missing")

    not_directory = tmp_path / "network.bnet"
    not_directory.write_text("A, 1\n")
    with pytest.raises(NotADirectoryError, match="path is not a directory"):
        bt.bpy.bn.read_bnet_directory(not_directory)

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no '.bnet' file found"):
        bt.bpy.bn.read_bnet_directory(empty)

    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "sub").mkdir()
    (nested / "sub" / "model.bnet").write_text("A, 1\n")

    ensemble = bt.bpy.bn.read_bnet_directory(nested, recursive=True)
    assert len(ensemble) == 1
    assert ensemble.components == frozenset({"A"})


def test_boolean_network_ensemble_initialization_and_slice_mutation_errors():
    with pytest.raises(TypeError, match="either 'components' or 'bns'"):
        bt.bpy.bn.BooleanNetworkEnsemble()

    with pytest.raises(TypeError, match="mutually exclusive"):
        bt.bpy.bn.BooleanNetworkEnsemble(components=["A"], bns=[{"A": 1}])

    with pytest.raises(ValueError, match="empty Boolean network collection"):
        bt.bpy.bn.BooleanNetworkEnsemble(bns=[])

    with pytest.raises(TypeError, match="Boolean network-like"):
        bt.bpy.bn.BooleanNetworkEnsemble(bns=[cast(Any, object())])

    ensemble = bt.bpy.bn.BooleanNetworkEnsemble(components=["A", "B"])
    assert len(ensemble) == 0
    assert ensemble.components == frozenset({"A", "B"})
    assert ensemble.ba is getattr(ensemble, "_BooleanNetworkEnsemble__ba")
    with pytest.raises(AttributeError):
        setattr(ensemble, "ba", object())

    first = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})
    second = bt.bpy.bn.BooleanNetwork({"A": 0, "B": "A"})
    ensemble.insert(0, first)
    ensemble[0:1] = [second]

    assert len(ensemble) == 1
    assert ensemble[0].rules == {"A": "0", "B": "A"}

    del ensemble[0]
    assert len(ensemble) == 0


def test_boolean_network_ensemble_regulator_counts(bnet_ensemble):
    counts = bnet_ensemble.regulator_counts()

    assert counts["A"]["B"][True] == 3
    assert counts["A"]["C"][True] == 1
    assert counts["A"]["C"][False] == 1


def test_boolean_network_ensemble_influence_counts(bnet_ensemble):
    influences = bnet_ensemble.influence_counts()

    assert influences["B"]["A"][True] == 3
    assert influences["C"]["A"][True] == 1
    assert influences["C"]["A"][False] == 1


def test_boolean_network_ensemble_rule_structures(bnet_ensemble):
    structures = bnet_ensemble.rule_structures()

    assert set(structures) == {"A", "B", "C"}

    assert len(structures["A"]) == 3
    assert structures["B"] == [True, True, True]
    assert structures["C"] == [False, False, False]

    assert len(set(structures["A"])) == 3


def test_boolean_network_ensemble_to_networkx(bnet_ensemble):
    graph = bnet_ensemble.to_networkx()

    assert set(graph.nodes) == {"A", "B", "C"}

    assert graph.nodes["A"]["function_count"] == 3
    assert graph.nodes["B"]["function_count"] == 1
    assert graph.nodes["C"]["function_count"] == 1

    assert graph.nodes["A"]["function_stability"] == pytest.approx(1 / 3)
    assert graph.nodes["B"]["function_stability"] == 1
    assert graph.nodes["C"]["function_stability"] == 1

    assert graph.has_edge("B", "A")
    assert graph.has_edge("C", "A")

    edge_data = list(graph.get_edge_data("B", "A").values())
    assert any(
        data["sign"] is True
        and data["count"] == 3
        and data["frequency"] == 1
        and "ratio" not in data
        for data in edge_data
    )

    edge_data = list(graph.get_edge_data("C", "A").values())
    assert any(
        data["sign"] is True
        and data["count"] == 1
        and data["frequency"] == pytest.approx(1 / 3)
        for data in edge_data
    )
    assert any(
        data["sign"] is False
        and data["count"] == 1
        and data["frequency"] == pytest.approx(1 / 3)
        for data in edge_data
    )


def test_boolean_network_ensemble_to_networkx_drop_isolates(bnet_ensemble):
    graph = bnet_ensemble.to_networkx(drop_isolates=True)

    assert set(graph.nodes) == {"A", "B", "C"}


def test_boolean_network_ensemble_to_graphviz_redirects_to_aggregated_graph(
    bnet_ensemble,
):
    with pytest.raises(NotImplementedError, match="from_boolean_networks"):
        bnet_ensemble.to_graphviz(rankdir="LR")


def test_boolean_network_ensemble_to_pydot_redirects_to_aggregated_graph(
    bnet_ensemble,
):
    with pytest.raises(NotImplementedError, match="from_boolean_networks"):
        bnet_ensemble.to_pydot(rankdir="LR")


def test_boolean_network_ensemble_show(monkeypatch, bnet_ensemble):
    calls = {}

    def fake_show(self, **kwargs):
        calls["graph"] = self
        calls["show"] = kwargs

    monkeypatch.setattr(bt.bpy.ig.AggregatedInfluenceGraph, "show", fake_show)

    bnet_ensemble.show(
        collapse="family",
        bins=None,
        preserve_feedback=False,
        include_selfloops=False,
        drop_isolates=True,
        min_frequency=0.5,
        edge_label="frequency",
        edge_style=None,
        program="neato",
        graph_attr={"rankdir": "LR"},
        node_attr={"shape": "box"},
        edge_attr={"fontsize": "10"},
        width="640px",
        height="480px",
    )

    graph = calls["graph"]
    assert isinstance(graph, bt.bpy.ig.AggregatedInfluenceGraph)
    assert graph.total == len(bnet_ensemble)

    assert calls["show"] == {
        "collapse": "family",
        "bins": None,
        "preserve_feedback": False,
        "include_selfloops": False,
        "min_frequency": 0.5,
        "drop_isolates": True,
        "node_style": "stability",
        "edge_label": "frequency",
        "edge_style": None,
        "program": "neato",
        "graph_attr": {"rankdir": "LR"},
        "node_attr": {"shape": "box"},
        "edge_attr": {"fontsize": "10"},
        "width": "640px",
        "height": "480px",
    }


def test_ensemble_allows_external_regulators_unchecked():
    bn = bt.bpy.bn.BooleanNetwork({"A": "X"}, check=False)
    ensemble = bt.bpy.bn.BooleanNetworkEnsemble(bns=[bn])

    assert ensemble.influence_counts()["X"]["A"][True] == 1


def test_boolean_network_ensemble_sequence_mutation(bnet_ensemble):
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "C",
            "B": 1,
            "C": 0,
        }
    )

    bnet_ensemble[0] = bn

    assert len(bnet_ensemble) == 3
    assert bnet_ensemble[0].rule("A") == "C"


def test_boolean_network_ensemble_rejects_invalid_components(bnet_ensemble):
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 1,
        }
    )

    with pytest.raises(ValueError):
        bnet_ensemble.append(bn)


def test_boolean_network_ensemble_rejects_invalid_network_type(bnet_ensemble):
    with pytest.raises(TypeError):
        bnet_ensemble.append(object())


def test_boolean_network_ensemble_copies_inserted_networks(bnet_ensemble):
    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 1,
            "C": 0,
        }
    )

    bnet_ensemble[0] = bn
    del bn["C"]

    assert "C" in bnet_ensemble[0]
