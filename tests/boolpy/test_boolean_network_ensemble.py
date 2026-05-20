#!/usr/bin/env python

#!/usr/bin/env python

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
    assert any(data["sign"] is True and data["count"] == 3 for data in edge_data)

    edge_data = list(graph.get_edge_data("C", "A").values())
    assert any(data["sign"] is True and data["count"] == 1 for data in edge_data)
    assert any(data["sign"] is False and data["count"] == 1 for data in edge_data)


def test_boolean_network_ensemble_to_networkx_remove_isolated_nodes(bnet_ensemble):
    graph = bnet_ensemble.to_networkx(remove_isolated_nodes=True)

    assert set(graph.nodes) == {"A", "B", "C"}


def test_boolean_network_ensemble_to_pydot(bnet_ensemble):
    pytest.importorskip("pydot")

    dot = bnet_ensemble.to_pydot(
        remove_isolated_nodes=True,
        show_edge_labels=False,
        node_style="stability",
    )

    edges = {
        (edge.get_source().strip('"'), edge.get_destination().strip('"'))
        for edge in dot.get_edges()
    }

    assert ("B", "A") in edges
    assert ("C", "A") in edges

    for edge in dot.get_edges():
        assert edge.get_label() is None


def test_boolean_network_ensemble_to_pydot_with_edge_labels(bnet_ensemble):
    pytest.importorskip("pydot")

    dot = bnet_ensemble.to_pydot(show_edge_labels=True)

    labels = {edge.get_label().strip('"') for edge in dot.get_edges()}

    assert "3" in labels
    assert "1" in labels


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
