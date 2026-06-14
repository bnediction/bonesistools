#!/usr/bin/env python

#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, cast

import pytest

import bonesistools as bt


def _pydot_get(obj, method):
    return getattr(cast(Any, obj), method)()


def _pydot_get_string(obj, method):
    return cast(str, _pydot_get(obj, method)).strip('"')


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
    assert any(data["sign"] is True and data["count"] == 3 for data in edge_data)

    edge_data = list(graph.get_edge_data("C", "A").values())
    assert any(data["sign"] is True and data["count"] == 1 for data in edge_data)
    assert any(data["sign"] is False and data["count"] == 1 for data in edge_data)


def test_boolean_network_ensemble_to_networkx_remove_isolated_nodes(bnet_ensemble):
    graph = bnet_ensemble.to_networkx(remove_isolated_nodes=True)

    assert set(graph.nodes) == {"A", "B", "C"}


def test_boolean_network_ensemble_to_graphviz(bnet_ensemble, fake_graphviz):
    graph = bnet_ensemble.to_graphviz(
        remove_isolated_nodes=True,
        show_edge_labels=False,
        node_style="stability",
        rankdir="LR",
    )

    edges = {(source, target): attrs for source, target, attrs in graph.edges}
    nodes = {node: attrs for node, attrs in graph.nodes}

    assert isinstance(graph, fake_graphviz)
    assert graph.graph_attr["rankdir"] == "LR"
    assert ("B", "A") in edges
    assert ("C", "A") in edges
    assert all("label" not in attrs for attrs in edges.values())
    assert all("fillcolor" in attrs for attrs in nodes.values())


def test_ensemble_to_graphviz_styles_thresholds(
    bnet_ensemble,
    fake_graphviz,
):
    with pytest.raises(ValueError, match="invalid argument value for 'min_ratio'"):
        bnet_ensemble.to_graphviz(min_ratio=2)

    graph = bnet_ensemble.to_graphviz(
        node_style="count",
        edge_style=True,
        min_ratio=0.4,
        rankdir="LR",
    )

    assert isinstance(graph, fake_graphviz)
    assert graph.graph_attr["rankdir"] == "LR"
    assert all("style" in attrs for _, _, attrs in graph.edges)
    assert all(attrs["shape"] == "oval" for _, attrs in graph.nodes)


def test_boolean_network_ensemble_to_pydot(bnet_ensemble):
    pytest.importorskip("pydot")

    dot = bnet_ensemble.to_pydot(
        remove_isolated_nodes=True,
        show_edge_labels=False,
        node_style="stability",
    )

    edges = {
        (
            _pydot_get_string(edge, "get_source"),
            _pydot_get_string(edge, "get_destination"),
        )
        for edge in dot.get_edges()
    }

    assert ("B", "A") in edges
    assert ("C", "A") in edges

    for edge in dot.get_edges():
        assert _pydot_get(edge, "get_label") is None


def test_boolean_network_ensemble_to_pydot_styles_thresholds_and_options(bnet_ensemble):
    pytest.importorskip("pydot")

    with pytest.raises(ValueError, match="invalid argument value for 'min_ratio'"):
        bnet_ensemble.to_pydot(min_ratio=2)

    dot = bnet_ensemble.to_pydot(
        node_style="count",
        edge_style=True,
        min_ratio=0.4,
        rankdir="LR",
    )

    assert cast(Any, dot).get_rankdir() == "LR"
    assert all(_pydot_get(edge, "get_style") is not None for edge in dot.get_edges())
    assert all(_pydot_get(node, "get_shape") == "oval" for node in dot.get_nodes())


def test_boolean_network_ensemble_to_pydot_with_edge_labels(bnet_ensemble):
    pytest.importorskip("pydot")

    dot = bnet_ensemble.to_pydot(show_edge_labels=True)

    labels = {_pydot_get_string(edge, "get_label") for edge in dot.get_edges()}

    assert "3" in labels
    assert "1" in labels


def test_boolean_network_ensemble_show(monkeypatch, bnet_ensemble):
    calls = {}

    class FakeDot:
        def create_svg(self):
            return b'<svg width="100pt" height="200pt"></svg>'

    def fake_to_pydot(
        self,
        remove_isolated_nodes=False,
        node_style=None,
        min_ratio=0.0,
        show_edge_labels=True,
        edge_style=None,
        program="dot",
        **kwargs,
    ):
        calls["to_pydot"] = {
            "remove_isolated_nodes": remove_isolated_nodes,
            "node_style": node_style,
            "min_ratio": min_ratio,
            "show_edge_labels": show_edge_labels,
            "edge_style": edge_style,
            "program": program,
            "kwargs": kwargs,
        }
        return FakeDot()

    class FakeSVG:
        def __init__(self, svg):
            self.svg = svg

    def fake_display(svg):
        calls["display"] = svg

    display_module = ModuleType("IPython.display")
    setattr(display_module, "SVG", FakeSVG)
    setattr(display_module, "display", fake_display)
    monkeypatch.setitem(sys.modules, "IPython.display", display_module)
    monkeypatch.setattr(
        bt.bpy.bn.BooleanNetworkEnsemble,
        "to_pydot",
        fake_to_pydot,
    )

    bnet_ensemble.show(
        remove_isolated_nodes=True,
        node_style="count",
        min_ratio=0.5,
        show_edge_labels=False,
        edge_style=False,
        program="neato",
        height="480px",
        rankdir="LR",
    )

    assert calls["to_pydot"] == {
        "remove_isolated_nodes": True,
        "node_style": "count",
        "min_ratio": 0.5,
        "show_edge_labels": False,
        "edge_style": False,
        "program": "neato",
        "kwargs": {"rankdir": "LR"},
    }
    assert isinstance(calls["display"], FakeSVG)
    assert calls["display"].svg == '<svg height="480px"></svg>'


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
