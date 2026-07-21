#!/usr/bin/env python

import os
import shutil
import subprocess
import sys
from textwrap import dedent
from typing import Any, cast

import numpy as np
import pytest

import bonesistools as bt
from bonesistools.logic.boolean_network import _distances


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


def test_boolean_network_ensemble_to_mpbn(bnet_ensemble):
    mpbn = pytest.importorskip("mpbn")

    converted = bnet_ensemble.to_mpbn()

    assert isinstance(converted, list)
    assert len(converted) == len(bnet_ensemble)
    assert all(isinstance(bn, mpbn.MPBooleanNetwork) for bn in converted)
    assert all(set(bn) == {"A", "B", "C"} for bn in converted)


def test_boolean_network_ensemble_to_biolqm(bnet_ensemble):
    biolqm = pytest.importorskip("biolqm")

    converted = bnet_ensemble.to_biolqm()

    assert len(converted) == len(bnet_ensemble)
    assert all(biolqm.is_biolqm_object(model) for model in converted)
    assert all(
        {str(component.getNodeID()) for component in model.getComponents()}
        == {"A", "B", "C"}
        for model in converted
    )


def test_boolean_network_ensemble_to_pyboolnet(bnet_ensemble):
    converted = bnet_ensemble.to_pyboolnet()

    assert len(converted) == len(bnet_ensemble)
    assert all(set(primes) == {"A", "B", "C"} for primes in converted)


def test_boolean_network_ensemble_to_minibn(bnet_ensemble):
    minibn = pytest.importorskip("colomoto.minibn")

    converted = bnet_ensemble.to_minibn()

    assert isinstance(converted, list)
    assert len(converted) == len(bnet_ensemble)
    assert all(isinstance(bn, minibn.BooleanNetwork) for bn in converted)
    assert all(set(bn) == {"A", "B", "C"} for bn in converted)


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
    with pytest.raises(TypeError, match="either 'components' or at least one network"):
        bt.logic.bn.BooleanNetworkEnsemble()

    with pytest.raises(TypeError, match="mutually exclusive"):
        bt.logic.bn.BooleanNetworkEnsemble({"A": 1}, components=["A"])

    with pytest.raises(TypeError, match="Boolean network-like"):
        bt.logic.bn.BooleanNetworkEnsemble(cast(Any, object()))

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


def test_boolean_network_ensemble_reuses_identical_coerced_rules():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(
        {"A": "B & ~C", "B": 1, "C": 0},
        {"A": "B & ~C", "B": 1, "C": 0},
    )

    assert ensemble[0]["A"] is ensemble[1]["A"]
    assert ensemble[0]["B"] is ensemble[1]["B"]
    assert ensemble[0]["C"] is ensemble[1]["C"]


def test_boolean_network_ensemble_simplify_reduces_all_rules_in_place():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(
        {"A": "B | (B & C)", "B": 0, "C": 1},
        {"A": "B & 1", "B": 1, "C": 0},
    )

    assert ensemble.simplify() is None
    assert [bn.rule("A") for bn in ensemble] == ["B", "B"]


def test_boolean_network_ensemble_distance_returns_expected_matrices():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(
        {"A": "A", "B": "B"},
        {"A": "A & B", "B": "B"},
        {"A": "~A", "B": "B"},
    )

    np.testing.assert_array_equal(
        ensemble.distance(metric="hamming"),
        np.array(
            [
                [0.0, 1 / 8, 1 / 2],
                [1 / 8, 0.0, 3 / 8],
                [1 / 2, 3 / 8, 0.0],
            ]
        ),
    )
    np.testing.assert_array_equal(
        ensemble.distance(metric="equivalence"),
        np.array(
            [
                [0.0, 1 / 2, 1 / 2],
                [1 / 2, 0.0, 1 / 2],
                [1 / 2, 1 / 2, 0.0],
            ]
        ),
    )


def test_boolean_network_ensemble_distance_is_semantic():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(
        {"A": "A", "B": "B"},
        {"A": "A | (A & B)", "B": "B"},
    )

    np.testing.assert_array_equal(
        ensemble.distance(metric="hamming"),
        np.zeros((2, 2)),
    )
    np.testing.assert_array_equal(
        ensemble.distance(metric="equivalence"),
        np.zeros((2, 2)),
    )


def test_boolean_network_ensemble_distance_preserves_invariant_denominator():
    variables = tuple(f"X{index}" for index in range(63))
    invariant_rules = {variable: variable for variable in variables}
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(
        {
            **invariant_rules,
            "T": " & ".join(variables),
        },
        {
            **invariant_rules,
            "T": 0,
        },
    )

    np.testing.assert_array_equal(
        ensemble.distance(metric="equivalence"),
        np.array([[0.0, 1 / 64], [1 / 64, 0.0]]),
    )
    np.testing.assert_array_equal(
        ensemble.distance(metric="hamming"),
        np.array([[0.0, 1 / (64 * 2**63)], [1 / (64 * 2**63), 0.0]]),
    )


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_boolean_network_ensemble_sparse_distance_matches_dense(
    bnet_ensemble,
    metric,
    monkeypatch,
):
    monkeypatch.setattr(
        _distances,
        "_use_sparse_equivalence_matrix",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        _distances,
        "_use_sparse_hamming_matrix",
        lambda *args, **kwargs: False,
    )
    dense = bnet_ensemble.distance(metric=metric)

    monkeypatch.setattr(
        _distances,
        "_use_sparse_equivalence_matrix",
        lambda *args, **kwargs: True,
    )
    monkeypatch.setattr(
        _distances,
        "_use_sparse_hamming_matrix",
        lambda *args, **kwargs: True,
    )
    sparse = bnet_ensemble.distance(metric=metric)

    np.testing.assert_array_equal(sparse, dense)


@pytest.mark.parametrize("metric", ["equivalence", "hamming"])
def test_boolean_network_ensemble_distance_matches_pairwise_function(
    bnet_ensemble,
    metric,
):
    matrix = bnet_ensemble.distance(metric=metric)

    assert matrix.dtype == np.float64
    assert matrix.shape == (len(bnet_ensemble), len(bnet_ensemble))
    for left_index, left in enumerate(bnet_ensemble):
        for right_index, right in enumerate(bnet_ensemble):
            assert matrix[left_index, right_index] == bt.logic.bn.distance(
                left,
                right,
                metric=metric,
            )


def test_boolean_network_ensemble_distance_handles_empty_and_singleton_ensembles():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=("A",))
    empty = ensemble.distance()

    assert isinstance(empty, np.ndarray)
    assert empty.dtype == np.float64
    assert empty.shape == (0, 0)
    assert empty.size == 0

    ensemble.insert(0, {"A": "A"})
    np.testing.assert_array_equal(
        ensemble.distance(),
        np.array([[0.0]]),
    )


def test_boolean_network_ensemble_distance_rejects_unsupported_metric():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=("A",))

    with pytest.raises(ValueError, match="unsupported Boolean-network"):
        ensemble.distance(metric=cast(Any, "unsupported"))


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
    assert list(graph.nodes) == ["A", "B", "C"]
    assert [
        (source, target, data["sign"])
        for source, target, data in graph.edges(data=True)
    ] == [
        ("B", "A", 1),
        ("C", "A", -1),
        ("C", "A", 1),
    ]

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
        data["sign"] == 1 and data["count"] == 1 and "frequency" not in data
        for data in edge_data
    )
    assert any(
        data["sign"] == -1 and data["count"] == 1 and "frequency" not in data
        for data in edge_data
    )
    assert graph.edge_frequency("C", "A", sign=1) == pytest.approx(1 / 3)
    assert graph.edge_frequency("C", "A", sign=-1) == pytest.approx(1 / 3)


def test_boolean_network_ensemble_graphviz_layout_is_hash_seed_independent():
    pytest.importorskip("pydot")

    if shutil.which("dot") is None:
        pytest.skip("Graphviz 'dot' executable is unavailable")

    script = dedent("""
        import bonesistools as bt

        networks = [
            {
                "APEX": "BCL2 | MCL1",
                "BCL2": "MCL1 & ~CASP3",
                "CASP3": "BAX & ~BCL2",
                "BAX": "TP53 | CASP3",
                "MCL1": "AKT & ~TP53",
                "AKT": "PI3K",
                "G1": "APEX",
                "G2": "APEX",
                "OUT": "G1 | G2",
                "PI3K": "APEX",
                "TP53": "CASP3 | MCL1",
            },
            {
                "APEX": "BCL2 & MCL1",
                "BCL2": "MCL1 | ~CASP3",
                "CASP3": "BAX & ~BCL2",
                "BAX": "TP53 & CASP3",
                "MCL1": "AKT | ~TP53",
                "AKT": "PI3K",
                "G1": "APEX",
                "G2": "APEX",
                "OUT": "G1 | G2",
                "PI3K": "APEX",
                "TP53": "CASP3 & MCL1",
            },
        ]

        graph = bt.logic.bn.BooleanNetworkEnsemble(*networks).to_influence_graph()

        for collapse in (None, "family", "feedback", "both"):
            dot = graph.to_pydot(
                collapse=collapse,
                edge_label=None,
                node_style=None,
                edge_style=None,
            )
            print(dot.create(format="plain").decode())
        """)

    layouts = []
    for seed in (1, 2):
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = str(seed)
        result = subprocess.run(
            [sys.executable, "-c", script],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
        layouts.append(result.stdout)

    assert layouts[0] == layouts[1]


def test_boolean_network_ensemble_to_influence_graph_drop_isolates(bnet_ensemble):
    graph = bnet_ensemble.to_influence_graph(drop_isolates=True)

    assert set(graph.nodes) == {"A", "B", "C"}


def test_boolean_network_ensemble_to_influence_graph_rejects_empty_ensemble():
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(components=("A",))

    with pytest.raises(ValueError, match="empty Boolean network ensemble"):
        ensemble.to_influence_graph()


def test_ensemble_allows_external_regulators_unchecked():
    bn = bt.logic.bn.BooleanNetwork({"A": "X"}, check=False)
    ensemble = bt.logic.bn.BooleanNetworkEnsemble(bn)

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
