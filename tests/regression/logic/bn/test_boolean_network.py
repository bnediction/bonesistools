#!/usr/bin/env python

import sys
from types import ModuleType
from typing import Any, cast

import pytest
from boolean import BooleanAlgebra, Expression

import bonesistools as bt
from bonesistools.logic.boolean_network import (
    _dynamics,
    _most_permissive,
    _network,
    _typing,
)


def _pydot_get(obj, method):
    return getattr(cast(Any, obj), method)()


def _pydot_get_string(obj, method):
    return cast(str, _pydot_get(obj, method)).strip('"')


def _edge_data(graph, source, target, key=0):
    return cast(Any, graph[source][target])[key]


def _configuration_from_bits(bits):
    return {f"x{index}": int(bit) for index, bit in enumerate(bits, start=1)}


def _configuration_bits(configuration):
    return "".join(str(configuration[f"x{index}"]) for index in range(1, 4))


def _assert_most_permissive_stg(bn, expected_targets):
    states = tuple(expected_targets)

    for source in states:
        enumerated = [
            _configuration_bits(configuration)
            for configuration in bn.reachable_configurations(
                _configuration_from_bits(source)
            )
        ]
        enumerated_targets = {target for target in enumerated if target != source}

        assert len(enumerated) == len(set(enumerated))
        assert set(enumerated) == expected_targets[source] | {source}
        assert enumerated_targets == expected_targets[source]

        for target in states:
            if target == source:
                continue

            for backend in ("auto", "asp", "hypercube"):
                assert bn.reachability(
                    _configuration_from_bits(source),
                    _configuration_from_bits(target),
                    backend=cast(Any, backend),
                ) is (target in expected_targets[source])


def test_boolean_network_coerces_string_rules():

    ba = BooleanAlgebra()
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & !C",
            "B": 0,
            "C": 1,
        },
        ba=ba,
    )

    assert bn.rule("A") == "B & ~C"
    assert bn.rule("B") == "0"
    assert bn.rule("C") == "1"


def test_boolean_network_coerces_boolean_constants():

    bn = bt.logic.bn.BooleanNetwork({"A": True, "B": False})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_integer_constants():

    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": 0})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_coerces_string_constants():

    bn = bt.logic.bn.BooleanNetwork({"A": "1", "B": "0"})

    assert bn["A"] is bn.ba.TRUE
    assert bn["B"] is bn.ba.FALSE


def test_boolean_network_accepts_boolean_algebra_constants():

    ba = BooleanAlgebra()

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": cast(Expression, ba.TRUE),
            "B": cast(Expression, ba.FALSE),
        },
        ba=ba,
    )

    assert bn["A"] is ba.TRUE
    assert bn["B"] is ba.FALSE


def test_boolean_network_rejects_invalid_rule_type():

    with pytest.raises(TypeError):
        bt.logic.bn.BooleanNetwork({"A": cast(Any, object())})


def test_boolean_network_rejects_undefined_symbols():
    with pytest.raises(ValueError):
        bt.logic.bn.BooleanNetwork({"A": "B"})


def test_boolean_network_can_be_created_unchecked():

    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    assert bn.components == {"A"}
    assert bn.symbols == {"B"}
    assert bn.undefined_symbols == {"B"}
    assert not bn.is_closed


def test_boolean_network_validate_raises_on_undefined_symbols():

    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError):
        bn.validate()


def test_boolean_network_components():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.components == {"A", "B", "C"}


def test_boolean_network_symbols():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.symbols == {"B", "C"}


def test_boolean_network_undefined_symbols_empty_for_closed_network():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.undefined_symbols == set()
    assert bn.is_closed


def test_boolean_network_rules_property():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.rules == {
        "A": "B & ~C",
        "B": "0",
        "C": "1",
    }


def test_boolean_network_setitem_coerces_rules():

    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    bn["B"] = "A & ~C"
    bn["C"] = 0

    assert bn.rule("B") == "A & ~C"
    assert bn["C"] is bn.ba.FALSE

    with pytest.raises(TypeError, match="unsupported argument type for 'component'"):
        bn[cast(Any, 1)] = 0


def test_boolean_network_string_representation():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert str(bn) == "A <- B & ~C\nB <- 0\nC <- 1"


def test_boolean_network_repr_is_compact():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert repr(bn) == "BooleanNetwork(components=3)"
    assert "\n" not in repr(bn)


def test_boolean_network_copy_preserves_type_algebra_and_unchecked_rules():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    copied = bn.copy()

    assert isinstance(copied, bt.logic.bn.BooleanNetwork)
    assert copied == bn
    assert copied is not bn
    assert copied.ba is bn.ba
    assert copied.undefined_symbols == {"B"}


def test_is_boolean_network_like_validates_mapping_rules():
    assert _typing.is_boolean_network_like({"A": "B", "B": 1})
    assert _typing.is_boolean_network_like(bt.logic.bn.BooleanNetwork({"A": 1}))

    assert not _typing.is_boolean_network_like({"A": object()})
    assert not _typing.is_boolean_network_like({1: "A"})
    assert not _typing.is_boolean_network_like(object())


def test_boolean_network_typing_namespace_is_not_exposed_publicly():
    assert "typing" not in dir(bt.logic.bn)

    with pytest.raises(AttributeError):
        getattr(bt.logic.bn, "typing")


def test_boolean_network_to_bnet_returns_string():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.to_bnet() == "A, B&!C\nB, 0\nC, 1\n"


def test_boolean_network_convert_to_mpbn():
    mpbn = pytest.importorskip("mpbn")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    converted = bn.convert("mpbn")

    assert isinstance(converted, mpbn.MPBooleanNetwork)
    assert set(converted) == {"A", "B", "C"}


def test_boolean_network_convert_to_minibn():
    minibn = pytest.importorskip("colomoto.minibn")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    converted = bn.convert("minibn.BooleanNetwork")
    converted_alias = bn.convert("minibn")

    assert isinstance(converted, minibn.BooleanNetwork)
    assert isinstance(converted_alias, minibn.BooleanNetwork)
    assert set(converted) == {"A", "B", "C"}
    assert set(converted_alias) == {"A", "B", "C"}


def test_boolean_network_convert_validates_target():
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(TypeError, match="unsupported argument type for 'target'"):
        bn.convert(cast(Any, 1))

    with pytest.raises(ValueError, match="unsupported conversion target"):
        bn.convert(cast(Any, "unknown"))


def test_boolean_network_rename_validates_inputs_and_collisions():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "A & C",
            "C": 1,
        }
    )

    assert bn.rename("A", "A") is None
    assert bn.rules["A"] == "B"

    with pytest.raises(KeyError, match="component 'missing' not found"):
        bn.rename("missing", "Y")

    with pytest.raises(ValueError, match="already exists"):
        bn.rename("A", "B")

    with pytest.raises(TypeError, match="unsupported argument type for 'old'"):
        bn.rename(cast(Any, 1), "Y")

    with pytest.raises(TypeError, match="unsupported argument type for 'new'"):
        bn.rename("A", cast(Any, 1))


def test_boolean_network_rename_validates_candidate_before_mutating():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "C"}, check=False)
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="undefined components"):
        bn.rename("A", "X")

    assert bn.rules == rules


def test_boolean_network_rename_accepts_named_old_new_arguments():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    bn.rename(old="A", new="X")

    assert bn.rules == {"X": "B", "B": "1"}


def test_boolean_network_relabel_renames_several_components_and_ignores_missing():
    bn = bt.logic.bn.BooleanNetwork({"Trp53": "Sox2", "Sox2": 1})

    assert bn.relabel({"Trp53": "TP53", "Sox2": "SOX2", "missing": "X"}) is None

    assert bn.rules == {"TP53": "SOX2", "SOX2": "1"}


def test_boolean_network_relabel_supports_atomic_swaps():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    bn.relabel({"A": "B", "B": "A"})

    assert bn.rules == {"B": "A", "A": "B"}


def test_boolean_network_relabel_rejects_component_merges_without_mutating():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})
    rules = bn.rules.copy()

    with pytest.raises(ValueError, match="merge Boolean network components"):
        bn.relabel({"A": "B"})

    assert bn.rules == rules


def test_boolean_network_relabel_validates_mapping():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    with pytest.raises(TypeError, match="unsupported argument type for 'mapping'"):
        bn.relabel(cast(Any, [("A", "X")]))

    with pytest.raises(TypeError, match="unsupported mapping key type"):
        bn.relabel(cast(Any, {1: "X"}))

    with pytest.raises(TypeError, match="unsupported mapping value type"):
        bn.relabel(cast(Any, {"A": 1}))

    assert bn.rules == {"A": "B", "B": "1"}


def test_boolean_network_from_bnet_and_to_bnet_file(tmp_path):
    infile = tmp_path / "network.bnet"
    outfile = tmp_path / "roundtrip.bnet"
    infile.write_text("# ignored\n\nA, B\nB, 1\n")

    bn = bt.logic.bn.BooleanNetwork.from_bnet(infile)
    read_bn = bt.logic.io.read_bnet(infile)

    assert bn.rules == {"A": "B", "B": "1"}
    assert read_bn.rules == bn.rules
    assert bn.to_bnet(outfile) is None
    assert outfile.read_text() == "A, B\nB, 1\n"


def test_deprecated_read_bnet_routes_to_io(tmp_path):
    infile = tmp_path / "network.bnet"
    infile.write_text("A, B\nB, 1\n")

    with pytest.warns(FutureWarning, match="bt.logic.bn.read_bnet"):
        bn = bt.logic.bn.read_bnet(infile)

    assert bn.rules == bt.logic.io.read_bnet(infile).rules


def test_boolean_network_structural_equality():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "C | B",
            "B": 0,
            "C": 1,
        }
    )

    assert bn1 == bn2


def test_boolean_network_structural_inequality():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & (C | ~C)",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": 0,
        }
    )

    assert bn1 != bn2
    assert bn1.__eq__(object()) is NotImplemented
    assert bn1.__ne__(object()) is NotImplemented


def test_boolean_network_equivalence_truth_table_only():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "(B & C) | (~B & D) | (C & D)",
            "B": 0,
            "C": 1,
            "D": 0,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "(B & C) | (~B & D)",
            "B": 0,
            "C": 1,
            "D": 0,
        }
    )

    assert not bn1.equivalent(bn2, method="simplify")
    assert bn1.equivalent(bn2, method="truth_table")


def test_boolean_network_not_equivalent_truth_table():

    bn1 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & C",
            "B": 0,
            "C": 1,
        }
    )

    bn2 = bt.logic.bn.BooleanNetwork(
        {
            "A": "B | C",
            "B": 0,
            "C": 1,
        }
    )

    assert not bn1.equivalent(bn2, method="truth_table")


def test_boolean_network_equivalence_requires_same_components():

    bn1 = bt.logic.bn.BooleanNetwork({"A": 1})
    bn2 = bt.logic.bn.BooleanNetwork({"A": 1, "B": 0})

    assert not bn1.equivalent(bn2)
    assert bn1.equivalent(object()) is NotImplemented


def test_boolean_network_equivalence_rejects_unknown_method():

    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(ValueError):
        bn.equivalent(bn, method=cast(Any, "unknown"))


def test_boolean_network_influences():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    assert bn.influences() == {
        ("B", "A", 1),
        ("C", "A", -1),
    }


def test_boolean_network_fixed_points_and_predicate():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.fixed_points() == [
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
    ]
    assert bn.fixed_points(limit=1) == [{"A": 0, "B": 0}]
    assert bn.is_fixed_point({"A": True, "B": True}) is True
    assert bn.is_fixed_point({"A": 1, "B": 0}) is False

    no_fixed_point = bt.logic.bn.BooleanNetwork({"A": "~A"})
    assert no_fixed_point.fixed_points() == []


def test_boolean_network_trapspaces_minimal_fixed_points():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.trapspaces() == (
        bt.logic.ba.Hypercube({"A": 0, "B": 0}),
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
    )


def test_boolean_network_trapspaces_universal_when_no_fixed_subspace():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    assert bn.trapspaces() == (bt.logic.ba.Hypercube({}),)


def test_boolean_network_trapspaces_rejects_invalid_options():
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(ValueError, match="invalid argument value for 'kind'"):
        bn.trapspaces(kind=cast(Any, "maximal"))

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.trapspaces(backend=cast(Any, "truth_table"))


def test_boolean_network_trapspaces_rejects_open_network():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="undefined components"):
        bn.trapspaces()


def test_boolean_network_trapspaces_asp_reports_missing_clingo(monkeypatch):
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    def missing_clingo(name: str) -> Any:
        if name == "clingo":
            raise ImportError("missing clingo")
        return __import__(name)

    monkeypatch.setattr(_network, "import_module", missing_clingo)

    with pytest.raises(ImportError, match="requires `clingo` to be installed"):
        bn.trapspaces()


def test_boolean_network_reachable_attractors_synchronous_cycle():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    attractors = bn.reachable_attractors(
        {"A": 0},
        update="synchronous",
        backend="explicit",
    )

    assert len(attractors) == 1
    assert isinstance(attractors[0], bt.logic.ba.ConfigurationSet)
    assert attractors[0].enumerate() == (
        {"A": 0},
        {"A": 1},
    )


def test_boolean_network_reachable_attractors_synchronous_partial_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "B"})

    attractors = bn.reachable_attractors(
        {"A": "*"},
        update="synchronous",
        backend="explicit",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        (
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
        ),
        (
            {"A": 0, "B": 1},
            {"A": 1, "B": 1},
        ),
    )


def test_boolean_network_reachable_attractors_synchronous_bdd_matches_explicit():
    pytest.importorskip("dd.autoref")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "~A",
            "B": "B",
            "C": "A | C",
        }
    )

    explicit = bn.reachable_attractors(
        {"A": "*", "B": 1, "C": 0},
        update="synchronous",
        backend="explicit",
    )
    bdd = bn.reachable_attractors(
        {"A": "*", "B": 1, "C": 0},
        update="synchronous",
        backend="bdd",
    )

    assert tuple(attractor.enumerate() for attractor in bdd) == tuple(
        attractor.enumerate() for attractor in explicit
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_attractors_bdd_matches_explicit(update):
    pytest.importorskip("dd.autoref")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "A",
            "C": "~C",
        }
    )

    explicit = bn.reachable_attractors(
        {"A": 0, "C": "*"},
        update=cast(Any, update),
        backend="explicit",
    )
    bdd = bn.reachable_attractors(
        {"A": 0, "C": "*"},
        update=cast(Any, update),
        backend="bdd",
    )

    assert tuple(attractor.enumerate() for attractor in bdd) == tuple(
        attractor.enumerate() for attractor in explicit
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_attractors_bdd_reports_missing_dd(
    monkeypatch,
    update,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    def missing_dd(name: str) -> Any:
        if name == "dd.autoref":
            raise ImportError("missing dd")
        return __import__(name)

    monkeypatch.setattr(_dynamics, "import_module", missing_dd)

    with pytest.raises(ImportError, match="backend='bdd'"):
        bn.reachable_attractors(
            {"A": 0},
            update=cast(Any, update),
            backend="bdd",
        )


def test_boolean_network_reachable_attractors_asynchronous_branching():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 1},
        update="asynchronous",
        backend="explicit",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 0, "B": 0},),
        ({"A": 1, "B": 1},),
    )


def test_boolean_network_reachable_attractors_asynchronous_partial_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0},
        update="asynchronous",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 0, "B": 0},),
        ({"A": 1, "B": 1},),
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_attractors_defaults_to_fully_free_state(update):
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    omitted = bn.reachable_attractors(
        update=cast(Any, update),
        backend="explicit",
    )
    explicit = bn.reachable_attractors(
        {},
        update=cast(Any, update),
        backend="explicit",
    )

    assert tuple(attractor.enumerate() for attractor in omitted) == tuple(
        attractor.enumerate() for attractor in explicit
    )


def test_boolean_network_reachable_attractors_explores_shared_paths_once(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    calls = 0
    successor_state_bits = _dynamics._successor_state_bits

    def counted_successor_state_bits(*args, **kwargs):
        nonlocal calls
        calls += 1
        return successor_state_bits(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_successor_state_bits",
        counted_successor_state_bits,
    )

    attractors = bn.reachable_attractors(
        {},
        update="asynchronous",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 0, "B": 0},),
        ({"A": 1, "B": 1},),
    )
    assert calls == 4


def test_boolean_network_reachable_attractors_general_uses_subsets_of_updates():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 1},
        update="general",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 0, "B": 0},),
        ({"A": 1, "B": 1},),
    )


def test_boolean_network_smallest_closed_hypercube_uses_component_subset():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        {"A": 0, "B": 1},
        relaxed_components=("A",),
    ) == bt.logic.ba.Hypercube({"B": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        {"A": 0, "B": 1},
        relaxed_components=("A", "B"),
    ) == bt.logic.ba.Hypercube({})


def test_boolean_network_smallest_closed_hypercube_matches_cmsb22_figure_2():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": 1,
            "x2": "x1",
            "x3": "(~x1 & x2) | x3",
        }
    )
    initial_state = {"x1": 0, "x2": 0, "x3": 1}

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_state,
        relaxed_components=("x1",),
    ) == bt.logic.ba.Hypercube({"x2": 0, "x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_state,
        relaxed_components=("x1", "x2"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_state,
        relaxed_components=("x1", "x2", "x3"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    initial_state = {"x1": 0, "x2": 1, "x3": 1}

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_state,
        relaxed_components=("x1", "x2", "x3"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_state,
        relaxed_components=("x2", "x3"),
    ) == bt.logic.ba.Hypercube({"x1": 0, "x3": 1})


def test_boolean_network_reachability_matches_cmsb22_figure_2():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": 1,
            "x2": "x1",
            "x3": "(~x1 & x2) | x3",
        }
    )
    initial_state = {"x1": 0, "x2": 0, "x3": 1}

    assert bn.reachability(
        initial_state,
        {"x1": 1, "x2": 0, "x3": 1},
    )
    assert bn.reachability(
        initial_state,
        {"x1": 1, "x2": 1, "x3": 1},
    )
    assert not bn.reachability(
        initial_state,
        {"x1": 0, "x2": 0, "x3": 0},
    )
    assert not bn.reachability(
        initial_state,
        {"x1": 0, "x2": 1, "x3": 1},
    )


def test_boolean_network_mp_transition_space_matches_cmsb22_figure_3_left():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": 1,
            "x2": "x1",
            "x3": "(~x1 & x2) | x3",
        }
    )
    expected_targets = {
        "000": {"100", "101", "110", "111"},
        "101": {"111"},
        "100": {"110"},
        "110": set(),
        "111": set(),
    }

    _assert_most_permissive_stg(bn, expected_targets)


def test_boolean_network_mp_transition_space_matches_cmsb22_figure_3_right():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": "x1 & ~x3",
            "x2": "x1",
            "x3": "~x1",
        }
    )
    expected_targets = {
        "000": {"001"},
        "001": set(),
        "010": {"000", "001", "011"},
        "011": {"001"},
        "100": {"110"},
        "101": {"000", "001", "010", "011", "100", "110", "111"},
        "110": set(),
        "111": {"000", "001", "010", "011", "100", "101", "110"},
    }

    _assert_most_permissive_stg(bn, expected_targets)


def test_boolean_network_smallest_closed_hypercube_accepts_depth_limit():
    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": "A"})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        {"A": 0, "B": 0},
        relaxed_components=("A", "B"),
        depth=1,
    ) == bt.logic.ba.Hypercube({"B": 0})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        {"A": 0, "B": 0},
        relaxed_components=("A", "B"),
        depth=None,
    ) == bt.logic.ba.Hypercube({})


def test_boolean_network_reachable_attractors_most_permissive():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": 1,
            "B": "A",
            "C": "(~A & B) | C",
        }
    )

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 0, "C": 0},
        update="most-permissive",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 1, "B": 1, "C": 0},),
        ({"A": 1, "B": 1, "C": 1},),
    )


def test_boolean_network_reachable_attractors_most_permissive_partial_state():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "AA": "AA",
            "BB": "AA | BB",
            "CC": "~AA | CC",
        }
    )

    attractors = bn.reachable_attractors(
        {"AA": 0},
        update="most-permissive",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"AA": 0, "BB": 0, "CC": 1},),
        ({"AA": 0, "BB": 1, "CC": 1},),
    )


def test_boolean_network_reachable_attractors_most_permissive_defaults_to_free_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"})

    omitted = bn.reachable_attractors(update="most-permissive")
    explicit = bn.reachable_attractors({}, update="most-permissive")

    assert tuple(attractor.enumerate() for attractor in omitted) == tuple(
        attractor.enumerate() for attractor in explicit
    )


def test_boolean_network_reachability_most_permissive_tracks_irreversibles():
    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": "A"})

    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 0},
    )
    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 0},
    )
    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
    )
    assert not bn.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
    )


def test_boolean_network_reachability_most_permissive_allows_reversible_space():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 0, "B": 0},
    )
    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 1, "B": 1},
    )
    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 1, "B": 0},
    )


def test_boolean_network_reachability_most_permissive_rejects_self_support():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    for backend in ("auto", "asp", "hypercube"):
        assert not bn.reachability(
            {"A": 0, "B": 0},
            {"A": 1, "B": 1},
            backend=cast(Any, backend),
        )


def test_boolean_network_reachability_validate_options_and_states():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachable_configurations(
            {"A": 0, "B": 0},
            update=cast(Any, "asynchronous"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.reachable_configurations(
            {"A": 0, "B": 0},
            backend=cast(Any, "asp"),
        )

    with pytest.raises(ValueError, match="initial_state must define fixed values"):
        list(bn.reachable_configurations({"A": 0}))

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachability(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            update=cast(Any, "asynchronous"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.reachability(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            backend=cast(Any, "explicit"),
        )

    with pytest.raises(ValueError, match="initial_state must define fixed values"):
        bn.reachability({"A": 0}, {"A": 1, "B": 0})

    with pytest.raises(ValueError, match="target_state must define fixed values"):
        bn.reachability({"A": 0, "B": 0}, {"A": 1})

    with pytest.raises(ValueError, match="unknown components in initial_state"):
        bn.reachability({"A": 0, "B": 0, "C": 0}, {"A": 1, "B": 0})


def test_boolean_network_reachable_attractors_compresses_large_terminal_scc():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "~B"})

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 0},
        update="asynchronous",
    )

    assert len(attractors) == 1
    assert attractors[0].enumerate() == (
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
        {"A": 1, "B": 0},
        {"A": 1, "B": 1},
    )
    assert len(attractors[0]._hypercubes) == 1


def test_boolean_network_reachable_attractors_validate_options():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachable_attractors({"A": 0}, update=cast(Any, "parallel"))

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.reachable_attractors(
            {"A": 0},
            update="synchronous",
            backend=cast(Any, "asp"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.reachable_attractors(
            {"A": 0},
            update="most-permissive",
            backend=cast(Any, "explicit"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'backend'"):
        bn.reachable_attractors(
            {"A": 0},
            update="asynchronous",
            backend=cast(Any, "asp"),
        )

    with pytest.raises(ValueError, match="unknown components"):
        bn.reachable_attractors({"B": 0})


def test_boolean_network_reachable_attractors_rejects_open_network():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="undefined components"):
        bn.reachable_attractors({"A": 0})


def test_boolean_network_next_methods():
    bn = bt.logic.bn.BooleanNetwork({"A": "B & ~C", "B": "A", "C": 0})

    assert bn.next_state("A", {"A": 0, "B": 1, "C": 0}) == 1
    assert bn.next_state("A", {"A": 0, "B": 1, "C": 1}) == 0

    assert bn.next_configuration({"A": 1, "B": 1, "C": 0}) == {
        "A": 1,
        "B": 1,
        "C": 0,
    }
    assert bn.next_configuration({"A": 0, "B": 1, "C": 1}) == {
        "A": 0,
        "B": 0,
        "C": 0,
    }

    with pytest.raises(ValueError, match="expected components"):
        bn.next_state("A", {"A": 0, "B": 1, "C": 0, "D": 1})

    unchecked = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="still depends on"):
        unchecked.next_state("A", {"A": 0})

    with pytest.raises(ValueError, match="still depends on"):
        unchecked.next_configuration({"A": 0})


def test_boolean_network_next_methods_evaluate_composed_constants():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "~((B & ~C) | (B & C))",
            "B": "B",
            "C": "C",
        }
    )

    assert bn.next_state("A", {"A": 0, "B": 0, "C": 0}) == 1
    assert bn.next_configuration({"A": 0, "B": 0, "C": 0}) == {
        "A": 1,
        "B": 0,
        "C": 0,
    }


def test_boolean_network_fixed_points_validate_inputs():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'limit'"):
        bn.fixed_points(limit=-1)

    with pytest.raises(ValueError, match="expected components"):
        bn.is_fixed_point({"A": 1})

    with pytest.raises(ValueError, match="expected components"):
        bn.is_fixed_point({"A": 1, "B": 1, "C": 0})

    assert bn.is_fixed_point({"A": 1, "B": "*"}) is False

    assert bn.is_fixed_point({"A": 1, "B": bt.logic.ba.PartialBoolean("*")}) is False

    with pytest.raises(TypeError, match="unsupported argument type for 'state'"):
        bn.is_fixed_point(cast(Any, object()))

    with pytest.raises(ValueError, match="expected string component names"):
        bn.is_fixed_point(cast(Any, {"A": 1, 2: 0}))


def test_boolean_network_to_influence_graph():

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_influence_graph()

    assert isinstance(graph, bt.logic.ig.InfluenceGraph)
    assert set(graph.nodes) == {"A", "B", "C"}
    assert graph.has_edge("B", "A")
    assert graph.has_edge("C", "A")

    assert _edge_data(graph, "B", "A")["sign"] == 1
    assert _edge_data(graph, "C", "A")["sign"] == -1


def test_boolean_network_to_graphviz(fake_graphviz):
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    graph = bn.to_graphviz(rankdir="LR")

    edges = {(source, target): attrs for source, target, attrs in graph.edges}

    assert isinstance(graph, fake_graphviz)
    assert graph.graph_attr["rankdir"] == "LR"
    assert ("B", "A") in edges
    assert ("C", "A") in edges
    assert edges[("B", "A")]["color"] == "green4"
    assert edges[("B", "A")]["arrowhead"] == "normal"
    assert edges[("C", "A")]["color"] == "red2"
    assert edges[("C", "A")]["arrowhead"] == "tee"


def test_boolean_network_to_pydot():

    pytest.importorskip("pydot")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B & ~C",
            "B": 0,
            "C": 1,
        }
    )

    dot = bn.to_pydot(rankdir="LR")

    assert cast(Any, dot).get_rankdir() == "LR"

    edges = {
        (
            _pydot_get_string(edge, "get_source"),
            _pydot_get_string(edge, "get_destination"),
        ): edge
        for edge in dot.get_edges()
    }

    assert ("B", "A") in edges
    assert ("C", "A") in edges

    assert _pydot_get(edges[("B", "A")], "get_color") == "green4"
    assert _pydot_get(edges[("B", "A")], "get_arrowhead") == "normal"

    assert _pydot_get(edges[("C", "A")], "get_color") == "red2"
    assert _pydot_get(edges[("C", "A")], "get_arrowhead") == "tee"


def test_deprecated_bn_to_pydot_is_not_promoted_from_bn_namespace():
    assert "bn_to_pydot" not in dir(bt.logic.bn)
    assert hasattr(bt.logic.bn, "bn_to_pydot")


def test_boolean_network_show(monkeypatch):
    calls = {}

    class FakeDot:
        def create_svg(self):
            return b'<svg width="100pt" height="200pt"></svg>'

    def fake_to_pydot(self, program="dot", edge_style=None, **kwargs):
        calls["to_pydot"] = {
            "program": program,
            "edge_style": edge_style,
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
    monkeypatch.setattr(bt.logic.bn.BooleanNetwork, "to_pydot", fake_to_pydot)

    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})

    bn.show(
        program="neato",
        edge_style=lambda sign: {"label": sign},
        width=600,
        rankdir="LR",
    )

    assert calls["to_pydot"]["program"] == "neato"
    assert calls["to_pydot"]["kwargs"] == {"rankdir": "LR"}
    assert isinstance(calls["display"], FakeSVG)
    assert calls["display"].svg == '<svg width="600"></svg>'
