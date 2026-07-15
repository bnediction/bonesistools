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
    return bt.logic.io.read_bnet_directory(bnet_directory)


def test_read_bnet_directory_and_ensemble(bnet_ensemble):
    assert len(bnet_ensemble) == 3
    assert bnet_ensemble.components == frozenset({"A", "B", "C"})


def test_boolean_network_ensemble_convert_to_mpbn(bnet_ensemble):
    mpbn = pytest.importorskip("mpbn")

    converted = bnet_ensemble.convert("mpbn")

    assert isinstance(converted, list)
    assert len(converted) == len(bnet_ensemble)
    assert all(isinstance(bn, mpbn.MPBooleanNetwork) for bn in converted)
    assert all(set(bn) == {"A", "B", "C"} for bn in converted)


def test_boolean_network_ensemble_convert_to_biolqm(bnet_ensemble):
    biolqm = pytest.importorskip("biolqm")

    converted = bnet_ensemble.convert("biolqm")

    assert len(converted) == len(bnet_ensemble)
    assert all(biolqm.is_biolqm_object(model) for model in converted)
    assert all(
        {str(component.getNodeID()) for component in model.getComponents()}
        == {"A", "B", "C"}
        for model in converted
    )


def test_boolean_network_ensemble_convert_to_pyboolnet(bnet_ensemble):
    converted = bnet_ensemble.convert("pyboolnet")

    assert len(converted) == len(bnet_ensemble)
    assert all(set(primes) == {"A", "B", "C"} for primes in converted)


def test_boolean_network_ensemble_convert_to_minibn(bnet_ensemble):
    minibn = pytest.importorskip("colomoto.minibn")

    converted = bnet_ensemble.convert("minibn")

    assert isinstance(converted, list)
    assert len(converted) == len(bnet_ensemble)
    assert all(isinstance(bn, minibn.BooleanNetwork) for bn in converted)
    assert all(set(bn) == {"A", "B", "C"} for bn in converted)


def test_boolean_network_ensemble_convert_validates_target(bnet_ensemble):
    with pytest.raises(TypeError, match="unsupported argument type for 'target'"):
        bnet_ensemble.convert(cast(Any, 1))

    with pytest.raises(ValueError, match="unsupported conversion target"):
        bnet_ensemble.convert(cast(Any, "unknown"))


def test_read_bnet_directory_validation_and_recursive_loading(tmp_path):
    with pytest.raises(FileNotFoundError, match="directory does not exist"):
        bt.logic.io.read_bnet_directory(tmp_path / "missing")

    not_directory = tmp_path / "network.bnet"
    not_directory.write_text("A, 1\n")
    with pytest.raises(NotADirectoryError, match="path is not a directory"):
        bt.logic.io.read_bnet_directory(not_directory)

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no '.bnet' file found"):
        bt.logic.io.read_bnet_directory(empty)

    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "sub").mkdir()
    (nested / "sub" / "model.bnet").write_text("A, 1\n")

    ensemble = bt.logic.io.read_bnet_directory(nested, recursive=True)
    assert len(ensemble) == 1
    assert ensemble.components == frozenset({"A"})


def test_deprecated_read_bnet_directory_routes_to_io(bnet_directory):
    with pytest.warns(FutureWarning, match="bt.logic.bn.read_bnet_directory"):
        ensemble = bt.logic.bn.read_bnet_directory(bnet_directory)

    expected = bt.logic.io.read_bnet_directory(bnet_directory)
    assert [bn.rules for bn in ensemble] == [bn.rules for bn in expected]


def test_boolean_network_ensemble_initialization_and_slice_mutation_errors():
    with pytest.raises(TypeError, match="either 'components' or 'bns'"):
        bt.logic.bn.BooleanNetworkEnsemble()

    with pytest.raises(TypeError, match="mutually exclusive"):
        bt.logic.bn.BooleanNetworkEnsemble(components=["A"], bns=[{"A": 1}])

    with pytest.raises(ValueError, match="empty Boolean network collection"):
        bt.logic.bn.BooleanNetworkEnsemble(bns=[])

    with pytest.raises(TypeError, match="Boolean network-like"):
        bt.logic.bn.BooleanNetworkEnsemble(bns=[cast(Any, object())])

    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=["A", "B"])
    assert len(ensemble) == 0
    assert ensemble.components == frozenset({"A", "B"})
    assert ensemble.ba is getattr(ensemble, "_BooleanNetworkEnsemble__ba")
    with pytest.raises(AttributeError):
        setattr(ensemble, "ba", object())

    first = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})
    second = bt.logic.bn.BooleanNetwork({"A": 0, "B": "A"})
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


def test_boolean_network_ensemble_dnf_implicants(bnet_ensemble):
    implicants = bnet_ensemble.dnf_implicants()

    assert set(implicants) == {"A", "B", "C"}

    assert len(implicants["A"]) == 3
    assert implicants["B"] == [
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
    ]
    assert implicants["C"] == [(), (), ()]

    assert len(set(implicants["A"])) == 3

    negative_implicants = bnet_ensemble.dnf_implicants(value=0)
    assert negative_implicants["B"] == [(), (), ()]
    assert negative_implicants["C"] == [
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
    ]


def test_boolean_network_ensemble_prime_implicants(bnet_ensemble):
    implicants = bnet_ensemble.prime_implicants()

    assert implicants["A"] == [
        (bt.logic.ba.Hypercube({"B": 1}),),
        (bt.logic.ba.Hypercube({"B": 1, "C": 1}),),
        (bt.logic.ba.Hypercube({"B": 1, "C": 0}),),
    ]
    assert implicants["B"] == [
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
        (bt.logic.ba.Hypercube(),),
    ]
    assert implicants["C"] == [(), (), ()]


@pytest.mark.parametrize("method", ["dnf_implicants", "prime_implicants"])
def test_boolean_network_ensemble_implicants_validate_value(method):
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=("A",))

    with pytest.raises(ValueError, match="expected 0 or 1"):
        getattr(ensemble, method)(value=2)


def test_boolean_network_ensemble_to_influence_graph(bnet_ensemble):
    graph = bnet_ensemble.to_influence_graph()

    assert isinstance(graph, bt.logic.ig.AggregatedInfluenceGraph)
    assert graph.total == 3
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
        data["sign"] == 1
        and data["count"] == 3
        and "frequency" not in data
        and "ratio" not in data
        for data in edge_data
    )
    assert graph.edge_frequency("B", "A", sign=1) == 1

    edge_data = list(graph.get_edge_data("C", "A").values())
    assert any(
        data["sign"] == 1
        and data["count"] == 1
        and "frequency" not in data
        for data in edge_data
    )
    assert any(
        data["sign"] == -1
        and data["count"] == 1
        and "frequency" not in data
        for data in edge_data
    )
    assert graph.edge_frequency("C", "A", sign=1) == pytest.approx(1 / 3)
    assert graph.edge_frequency("C", "A", sign=-1) == pytest.approx(1 / 3)


def test_boolean_network_ensemble_to_influence_graph_drop_isolates(bnet_ensemble):
    graph = bnet_ensemble.to_influence_graph(drop_isolates=True)

    assert set(graph.nodes) == {"A", "B", "C"}


def test_boolean_network_ensemble_to_influence_graph_rejects_empty_ensemble():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=("A",))

    with pytest.raises(ValueError, match="empty Boolean network ensemble"):
        ensemble.to_influence_graph()


def test_ensemble_allows_external_regulators_unchecked():
    bn = bt.logic.bn.BooleanNetwork({"A": "X"}, check=False)
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(bns=[bn])

    assert ensemble.influence_counts()["X"]["A"][True] == 1


def test_boolean_network_ensemble_sequence_mutation(bnet_ensemble):
    bn = bt.logic.bn.BooleanNetwork(
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
    bn = bt.logic.bn.BooleanNetwork(
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
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 1,
            "C": 0,
        }
    )

    bnet_ensemble[0] = bn
    del bn["C"]

    assert "C" in bnet_ensemble[0]
