#!/usr/bin/env python

from inspect import signature
from itertools import combinations, product
from typing import Any, cast

import pytest

import bonesistools as bt
from bonesistools.logic.boolean_network import _dynamics, _network


def _canonical_attractors(attractors):
    return frozenset(
        frozenset(tuple(sorted(configuration.items())) for configuration in attractor)
        for attractor in attractors
    )


def _canonical_configurations(configurations):
    return frozenset(
        tuple(sorted(configuration.items())) for configuration in configurations
    )


def _reachable_attractors_with_backend(bn, initial_state, *, update, backend):
    initial_states = bt.logic.ba.ConfigurationSet(tuple(bn.keys()), [initial_state])

    if backend == "bdd":
        return _dynamics._bdd_reachable_attractors(
            bn,
            initial_states,
            update=update,
        )

    if update == "synchronous":
        return _dynamics._synchronous_reachable_attractors(bn, initial_states)

    return _dynamics._explicit_reachable_attractors(
        bn,
        initial_states,
        update=update,
    )


def _explicit_successors(bn, state, update):
    next_state = bn.next_configuration(state)
    unstable = tuple(
        component for component in bn if state[component] != next_state[component]
    )

    if update == "synchronous":
        return (next_state,)

    successors = []
    updated_subsets = (
        ((component,) for component in unstable)
        if update == "asynchronous"
        else (
            subset
            for size in range(1, len(unstable) + 1)
            for subset in combinations(unstable, size)
        )
    )
    for updated in updated_subsets:
        successor = dict(state)
        for component in updated:
            successor[component] = next_state[component]
        successors.append(successor)

    return tuple(successors)


def _matches_hypercube(configuration, hypercube):
    return all(
        configuration[component] == value for component, value in hypercube.items()
    )


def _reference_reachable_configurations(bn, initial_state, update):
    initial_key = tuple(initial_state[component] for component in bn)
    pending = [dict(initial_state)]
    scheduled = {initial_key}
    configurations = []

    while pending:
        configuration = pending.pop()
        configurations.append(configuration)
        for successor in _explicit_successors(bn, configuration, update):
            key = tuple(successor[component] for component in bn)
            if key not in scheduled:
                scheduled.add(key)
                pending.append(successor)

    return configurations


def test_boolean_network_fixed_points_and_predicate():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.fixed_points() == [
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
    ]
    assert bn.fixed_points(limit=1) == [{"A": 0, "B": 0}]
    assert bn.is_fixed_point({"A": True, "B": True}) is True
    assert bn.is_fixed_point({"A": 1, "B": 0}) is False
    assert bn.is_fixed_point(bt.logic.ba.Hypercube({"A": 1, "B": 1})) is True

    no_fixed_point = bt.logic.bn.BooleanNetwork({"A": "~A"})
    assert no_fixed_point.fixed_points() == []


def test_boolean_network_reachable_attractors_synchronous_cycle():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    attractors = bn.reachable_attractors(
        {"A": 0},
        update="synchronous",
    )

    assert len(attractors) == 1
    assert isinstance(attractors[0], bt.logic.ba.ConfigurationSet)
    assert _canonical_attractors(attractors) == _canonical_attractors(
        (({"A": 0}, {"A": 1}),)
    )


def test_boolean_network_reachable_attractors_synchronous_partial_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "B"})

    attractors = bn.reachable_attractors(
        {"A": "*"},
        update="synchronous",
    )

    assert _canonical_attractors(attractors) == _canonical_attractors(
        (
            (
                {"A": 0, "B": 0},
                {"A": 1, "B": 0},
            ),
            (
                {"A": 0, "B": 1},
                {"A": 1, "B": 1},
            ),
        )
    )


def test_boolean_network_reachable_attractors_synchronous_reuses_trajectories(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": "A"})
    calls = 0
    next_state_bits = _dynamics._next_state_bits

    def counted_next_state_bits(*args, **kwargs):
        nonlocal calls
        calls += 1
        return next_state_bits(*args, **kwargs)

    monkeypatch.setattr(_dynamics, "_next_state_bits", counted_next_state_bits)
    monkeypatch.setattr(
        _dynamics,
        "_terminal_strongly_connected_components",
        lambda _successors: pytest.fail("synchronous dynamics must not compute SCCs"),
    )

    attractors = bn.reachable_attractors(
        {},
        update="synchronous",
    )

    assert len(attractors) == 1
    assert attractors[0] == bt.logic.ba.ConfigurationSet(
        ("A", "B"),
        [{"A": 1, "B": 1}],
    )
    assert calls == 4


def test_attractor_states_are_compacted_without_incremental_additions(monkeypatch):
    def fail(*_args, **_kwargs):
        pytest.fail("state bitsets must be compacted in bulk")

    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "add", fail)
    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "compress", fail)

    configurations = _dynamics._configuration_set_from_state_bits(
        {0b000, 0b001, 0b100, 0b101, 0b110},
        ("A", "B", "C"),
    )

    assert len(configurations._hypercubes) == 2
    assert {tuple(configuration.items()) for configuration in configurations} == {
        (("A", 0), ("B", 0), ("C", 0)),
        (("A", 1), ("B", 0), ("C", 0)),
        (("A", 0), ("B", 0), ("C", 1)),
        (("A", 1), ("B", 0), ("C", 1)),
        (("A", 0), ("B", 1), ("C", 1)),
    }


def test_boolean_network_reachable_attractors_synchronous_bdd_matches_explicit():
    pytest.importorskip("dd.autoref")

    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "~A",
            "B": "B",
            "C": "A | C",
        }
    )

    explicit = _reachable_attractors_with_backend(
        bn,
        {"A": "*", "B": 1, "C": 0},
        update="synchronous",
        backend="explicit",
    )
    bdd = _reachable_attractors_with_backend(
        bn,
        {"A": "*", "B": 1, "C": 0},
        update="synchronous",
        backend="bdd",
    )

    assert _canonical_attractors(bdd) == _canonical_attractors(explicit)


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

    explicit = _reachable_attractors_with_backend(
        bn,
        {"A": 0, "C": "*"},
        update=cast(Any, update),
        backend="explicit",
    )
    bdd = _reachable_attractors_with_backend(
        bn,
        {"A": 0, "C": "*"},
        update=cast(Any, update),
        backend="bdd",
    )

    assert _canonical_attractors(bdd) == _canonical_attractors(explicit)


def test_synchronous_explicit_attractors_consume_encoded_states(monkeypatch):
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    monkeypatch.setattr(
        bt.logic.ba.ConfigurationSet,
        "__iter__",
        lambda _self: pytest.fail("synchronous dynamics must consume bitsets"),
    )

    attractors = _reachable_attractors_with_backend(
        bn,
        {"A": "*"},
        update="synchronous",
        backend="explicit",
    )

    assert len(attractors) == 3
    assert sorted(len(attractor) for attractor in attractors) == [1, 1, 2]


def test_synchronous_bdd_attractors_use_functional_cycles(monkeypatch):
    pytest.importorskip("dd.autoref")
    components = tuple(f"x{index}" for index in range(8))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )

    monkeypatch.setattr(
        _dynamics._BDDTransitionSystem,
        "_terminal_sccs",
        lambda *_args, **_kwargs: pytest.fail(
            "synchronous dynamics must use functional cycles"
        ),
    )

    attractors = _reachable_attractors_with_backend(
        bn,
        {},
        update="synchronous",
        backend="bdd",
    )

    assert len(attractors) == 128
    assert all(len(attractor) == 2 for attractor in attractors)


def test_synchronous_bdd_attractors_remove_transient_states_functionally(
    monkeypatch,
):
    pytest.importorskip("dd.autoref")
    components = tuple(f"x{index}" for index in range(8))
    bn = bt.logic.bn.BooleanNetwork({component: 0 for component in components})

    monkeypatch.setattr(
        _dynamics._BDDTransitionSystem,
        "_terminal_sccs",
        lambda *_args, **_kwargs: pytest.fail(
            "synchronous dynamics must remove transients functionally"
        ),
    )

    attractors = _reachable_attractors_with_backend(
        bn,
        {},
        update="synchronous",
        backend="bdd",
    )

    assert len(attractors) == 1
    assert attractors[0].enumerate() == ({component: 0 for component in components},)


@pytest.mark.parametrize("update", ["asynchronous", "general"])
def test_bdd_reachable_attractors_keep_transient_states_symbolic(
    monkeypatch,
    update,
):
    pytest.importorskip("dd.autoref")
    components = tuple(f"x{index}" for index in range(8))
    bn = bt.logic.bn.BooleanNetwork(
        {
            **{component: f"~{component}" for component in components},
            "selector": "selector",
        }
    )

    monkeypatch.setattr(
        _dynamics,
        "_compile_bitset_rules",
        lambda *_args, **_kwargs: pytest.fail("BDD attractors must remain symbolic"),
    )
    monkeypatch.setattr(
        _dynamics,
        "_terminal_strongly_connected_components",
        lambda _successors: pytest.fail("BDD attractors must not build an STG"),
    )

    attractors = _reachable_attractors_with_backend(
        bn,
        {component: 0 for component in components},
        update=cast(Any, update),
        backend="bdd",
    )

    assert len(attractors) == 2
    assert all(len(attractor) == 256 for attractor in attractors)
    assert all(len(attractor._hypercubes) == 1 for attractor in attractors)


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_attractors_bdd_reports_missing_dd(
    monkeypatch,
    update,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    def missing_dd(name: str) -> Any:
        if name in {"dd.autoref", "dd.cudd"}:
            raise ImportError("missing dd")
        return __import__(name)

    monkeypatch.setattr(_dynamics, "import_module", missing_dd)

    with pytest.raises(ImportError, match=r"bonesistools\[bdd\]"):
        _reachable_attractors_with_backend(
            bn,
            {"A": 0},
            update=cast(Any, update),
            backend="bdd",
        )


@pytest.mark.parametrize(
    ("n_components", "expected_backend"),
    [(13, "explicit"), (14, "bdd")],
)
def test_asynchronous_attractor_backend_uses_closed_hypercube_bound(
    monkeypatch,
    n_components,
    expected_backend,
):
    components = tuple(f"x{index}" for index in range(n_components))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("bdd") or (),
    )

    bn.reachable_attractors(
        {component: 0 for component in components},
        update="asynchronous",
    )

    assert selected == [expected_backend]


@pytest.mark.parametrize(
    ("n_components", "expected_backend"),
    [(10, "explicit"), (11, "bdd")],
)
def test_general_attractor_backend_uses_closed_hypercube_bound(
    monkeypatch,
    n_components,
    expected_backend,
):
    components = tuple(f"x{index}" for index in range(n_components))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("bdd") or (),
    )

    bn.reachable_attractors(
        {component: 0 for component in components},
        update="general",
    )

    assert selected == [expected_backend]


@pytest.mark.parametrize(
    ("n_components", "expected_backend"),
    [(10, "explicit"), (11, "bdd")],
)
def test_synchronous_attractor_backend_uses_initial_state_bound(
    monkeypatch,
    n_components,
    expected_backend,
):
    components = tuple(f"x{index}" for index in range(n_components))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        _network,
        "_synchronous_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("bdd") or (),
    )

    bn.reachable_attractors({}, update="synchronous")

    assert selected == [expected_backend]


def test_asynchronous_attractor_backend_uses_initial_state_bound(monkeypatch):
    components = tuple(f"x{index}" for index in range(11))
    bn = bt.logic.bn.BooleanNetwork({component: component for component in components})
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("bdd") or (),
    )

    bn.reachable_attractors({}, update="asynchronous")

    assert selected == ["bdd"]


@pytest.mark.parametrize("update", ["asynchronous", "synchronous", "general"])
def test_automatic_attractor_backend_falls_back_without_dd(monkeypatch, update):
    components = tuple(f"x{index}" for index in range(14))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: None)
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_synchronous_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("explicit") or (),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_attractors",
        lambda *_args, **_kwargs: selected.append("bdd") or (),
    )

    bn.reachable_attractors(
        {component: 0 for component in components},
        update=cast(Any, update),
    )

    assert selected == ["explicit"]


@pytest.mark.parametrize(
    ("update", "n_components", "expected_backend"),
    [
        ("asynchronous", 12, "explicit"),
        ("asynchronous", 13, "bdd"),
        ("general", 8, "explicit"),
        ("general", 9, "bdd"),
    ],
)
def test_reachable_configurations_select_backend_from_closed_hypercube_bound(
    monkeypatch,
    update,
    n_components,
    expected_backend,
):
    components = tuple(f"x{index}" for index in range(n_components))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_configurations",
        lambda *_args, **_kwargs: selected.append("explicit") or iter(()),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_configurations",
        lambda *_args, **_kwargs: selected.append("bdd") or iter(()),
    )

    configurations = bn.reachable_configurations(
        {component: 0 for component in components},
        update=cast(Any, update),
    )

    assert tuple(configurations) == ()
    assert selected == [expected_backend]


@pytest.mark.parametrize("update", ["asynchronous", "general"])
def test_reachable_configurations_fall_back_to_explicit_without_dd(
    monkeypatch,
    update,
):
    components = tuple(f"x{index}" for index in range(14))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    selected = []

    monkeypatch.setattr(_network, "find_spec", lambda _name: None)
    monkeypatch.setattr(
        _network,
        "_explicit_reachable_configurations",
        lambda *_args, **_kwargs: selected.append("explicit") or iter(()),
    )
    monkeypatch.setattr(
        _network,
        "_bdd_reachable_configurations",
        lambda *_args, **_kwargs: selected.append("bdd") or iter(()),
    )

    configurations = bn.reachable_configurations(
        {component: 0 for component in components},
        update=cast(Any, update),
    )

    assert tuple(configurations) == ()
    assert selected == ["explicit"]


def test_boolean_network_reachable_attractors_asynchronous_branching():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 1},
        update="asynchronous",
    )

    assert _canonical_attractors(attractors) == _canonical_attractors(
        (
            ({"A": 0, "B": 0},),
            ({"A": 1, "B": 1},),
        )
    )


def test_boolean_network_reachable_attractors_asynchronous_partial_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0},
        update="asynchronous",
    )

    assert _canonical_attractors(attractors) == _canonical_attractors(
        (
            ({"A": 0, "B": 0},),
            ({"A": 1, "B": 1},),
        )
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_attractors_defaults_to_fully_free_state(update):
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    omitted = bn.reachable_attractors(
        update=cast(Any, update),
    )
    explicit = bn.reachable_attractors(
        {},
        update=cast(Any, update),
    )

    assert _canonical_attractors(omitted) == _canonical_attractors(explicit)


def test_boolean_network_reachable_attractors_explores_shared_paths_once(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    calls = 0
    next_state_bits = _dynamics._next_state_bits

    def counted_next_state_bits(*args, **kwargs):
        nonlocal calls
        calls += 1
        return next_state_bits(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_next_state_bits",
        counted_next_state_bits,
    )
    monkeypatch.setattr(
        _dynamics,
        "_terminal_strongly_connected_components",
        lambda _successors: pytest.fail("asynchronous dynamics must use Tarjan"),
    )

    attractors = bn.reachable_attractors(
        {},
        update="asynchronous",
    )

    assert _canonical_attractors(attractors) == _canonical_attractors(
        (
            ({"A": 0, "B": 0},),
            ({"A": 1, "B": 1},),
        )
    )
    assert calls == 4


def test_explicit_reachable_attractors_consumes_initial_states_lazily(monkeypatch):
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})
    explored = False
    next_state_bits = _dynamics._next_state_bits

    class InitialStates:
        def __iter__(self):
            yield {"A": 0}
            assert explored
            yield {"A": 1}

    def counted_next_state_bits(*args, **kwargs):
        nonlocal explored
        explored = True
        return next_state_bits(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_next_state_bits",
        counted_next_state_bits,
    )

    attractors = _dynamics._explicit_reachable_attractors(
        bn,
        cast(Any, InitialStates()),
        update="asynchronous",
    )

    assert len(attractors) == 2
    assert all(len(attractor) == 1 for attractor in attractors)


def test_boolean_network_reachable_attractors_general_uses_subsets_of_updates():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    attractors = bn.reachable_attractors(
        {"A": 0, "B": 1},
        update="general",
    )

    assert _canonical_attractors(attractors) == _canonical_attractors(
        (
            ({"A": 0, "B": 0},),
            ({"A": 1, "B": 1},),
        )
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_boolean_network_reachable_configurations_match_reference_stg(update):
    networks = (
        bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"}),
        bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"}),
        bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"}),
        bt.logic.bn.BooleanNetwork({"A": "A | B", "B": "A & ~B"}),
    )

    for bn in networks:
        for values in product((0, 1), repeat=len(bn)):
            initial_state = dict(zip(bn, values))
            expected = _reference_reachable_configurations(
                bn,
                initial_state,
                update,
            )
            observed = tuple(
                bn.reachable_configurations(
                    initial_state,
                    update=cast(Any, update),
                )
            )

            assert _canonical_configurations(observed) == _canonical_configurations(
                expected
            )
            assert len(observed) == len(_canonical_configurations(observed))


def test_boolean_network_reachable_configurations_distinguish_update_semantics():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    initial_state = {"A": 0, "B": 1}

    synchronous = _canonical_configurations(
        bn.reachable_configurations(initial_state, update="synchronous")
    )
    asynchronous = _canonical_configurations(
        bn.reachable_configurations(initial_state, update="asynchronous")
    )
    general = _canonical_configurations(
        bn.reachable_configurations(initial_state, update="general")
    )

    assert synchronous == _canonical_configurations(
        ({"A": 0, "B": 1}, {"A": 1, "B": 0})
    )
    assert asynchronous == _canonical_configurations(
        ({"A": 0, "B": 1}, {"A": 0, "B": 0}, {"A": 1, "B": 1})
    )
    assert general == _canonical_configurations(
        (
            {"A": 0, "B": 0},
            {"A": 0, "B": 1},
            {"A": 1, "B": 0},
            {"A": 1, "B": 1},
        )
    )


@pytest.mark.parametrize("update", ["asynchronous", "general"])
def test_bdd_reachable_configurations_match_explicit_traversal(update):
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "A",
            "C": "~C",
        }
    )
    initial_state = bt.logic.ba.Hypercube({"A": 0, "B": 1, "C": 0})

    explicit = _dynamics._explicit_reachable_configurations(
        bn,
        initial_state,
        update=cast(Any, update),
    )
    bdd = _dynamics._bdd_reachable_configurations(
        bn,
        initial_state,
        update=cast(Any, update),
    )

    assert _canonical_configurations(bdd) == _canonical_configurations(explicit)


def test_explicit_reachable_configurations_yield_initial_state_lazily(monkeypatch):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    monkeypatch.setattr(
        _dynamics,
        "_next_state_bits",
        lambda *_args, **_kwargs: pytest.fail(
            "successors must be computed after yielding the current state"
        ),
    )

    configurations = bn.reachable_configurations(
        {"A": 0},
        update="asynchronous",
    )
    assert next(configurations) == {"A": 0}


@pytest.mark.parametrize("update", ["asynchronous", "general"])
def test_explicit_reachable_configurations_stream_successors(monkeypatch, update):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "~B"})

    monkeypatch.setattr(
        _dynamics,
        "_successor_state_bits",
        lambda *_args, **_kwargs: pytest.fail(
            "reachable configurations must not materialize successor tuples"
        ),
    )

    configurations = tuple(
        bn.reachable_configurations(
            {"A": 0, "B": 0},
            update=cast(Any, update),
        )
    )
    assert len(configurations) == 4


def test_boolean_network_transition_distinguishes_update_semantics():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "~B"})
    initial = {"A": 0, "B": 0}

    assert bn.transition(
        initial,
        {"A": 1, "B": 1},
        update="synchronous",
    )
    assert not bn.transition(
        initial,
        {"A": 1, "B": 1},
        update="asynchronous",
    )
    assert bn.transition(
        initial,
        {"A": 1, "B": 0},
        update="asynchronous",
    )
    assert bn.transition(
        initial,
        {"A": 1, "B": 1},
        update="general",
    )


def test_boolean_network_transition_handles_identity_by_update_semantics():
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})
    state = {"A": 0}

    assert bn.transition(state, state, update="synchronous")
    assert not bn.transition(state, state, update="asynchronous")
    assert not bn.transition(state, state, update="general")


def test_boolean_network_transition_concrete_states_do_not_require_bdd(monkeypatch):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "~B"})

    def unexpected_bdd_import():
        raise AssertionError("concrete transitions should use bitsets")

    monkeypatch.setattr(_dynamics, "_import_dd_backend", unexpected_bdd_import)

    for update in ("synchronous", "asynchronous", "general"):
        assert bn.transition(
            {"A": 0, "B": 0},
            {"A": 1, "B": int(update != "asynchronous")},
            update=cast(Any, update),
        )

    assert bn.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="synchronous",
        quantifier="exists",
    )
    assert not bn.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="synchronous",
        quantifier="universal",
    )
    assert bn.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="asynchronous",
        quantifier="exists",
    )
    assert not bn.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="asynchronous",
        quantifier="universal",
    )
    assert bn.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="general",
        quantifier="universal",
    )


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_bdd_transition_quantifiers_match_explicit_relation(update):
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"})
    components = tuple(bn)
    states = tuple(
        dict(zip(components, values))
        for values in product((0, 1), repeat=len(components))
    )
    hypercubes = tuple(
        {
            component: value
            for component, value in zip(components, values)
            if value is not None
        }
        for values in product((None, 0, 1), repeat=len(components))
    )
    successors = {
        tuple(state.items()): _explicit_successors(bn, state, update)
        for state in states
    }
    transition_system = _dynamics._BDDTransitionSystem(
        bn,
        update=cast(Any, update),
    )

    for initial_hypercube, target_hypercube in product(hypercubes, repeat=2):
        initial_states = tuple(
            state for state in states if _matches_hypercube(state, initial_hypercube)
        )
        target_states = tuple(
            state for state in states if _matches_hypercube(state, target_hypercube)
        )
        relation = {
            tuple(initial.items()): {
                tuple(target.items()) for target in successors[tuple(initial.items())]
            }
            for initial in initial_states
        }
        expected = {
            "exists": any(
                tuple(target.items()) in relation[tuple(initial.items())]
                for initial in initial_states
                for target in target_states
            ),
            "robust": all(
                any(
                    tuple(target.items()) in relation[tuple(initial.items())]
                    for target in target_states
                )
                for initial in initial_states
            ),
            "universal": all(
                tuple(target.items()) in relation[tuple(initial.items())]
                for initial in initial_states
                for target in target_states
            ),
        }

        for quantifier, expected_result in expected.items():
            assert (
                transition_system.transition(
                    initial_hypercube,
                    target_hypercube,
                    quantifier=cast(Any, quantifier),
                )
                is expected_result
            )
            if len(initial_hypercube) == len(components):
                assert (
                    bn.transition(
                        initial_hypercube,
                        target_hypercube,
                        update=cast(Any, update),
                        quantifier=cast(Any, quantifier),
                    )
                    is expected_result
                )


def test_boolean_network_transition_quantifies_partial_states():
    pytest.importorskip("dd.autoref")
    oscillator = bt.logic.bn.BooleanNetwork({"A": "~A"})

    assert oscillator.transition(
        {},
        {"A": 1},
        update="asynchronous",
        quantifier="exists",
    )
    assert not oscillator.transition(
        {},
        {"A": 1},
        update="asynchronous",
        quantifier="robust",
    )
    assert oscillator.transition(
        {},
        {},
        update="asynchronous",
        quantifier="robust",
    )
    assert not oscillator.transition(
        {},
        {},
        update="asynchronous",
        quantifier="universal",
    )

    general = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "~B"})
    assert general.transition(
        {"A": 0, "B": 0},
        {"A": 1},
        update="general",
        quantifier="universal",
    )


def test_boolean_network_reachability_asynchronous_uses_bdd_backward_closure():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    initial = {"A": 0, "B": 1}

    assert bn.reachability(
        initial,
        {"A": 0, "B": 1},
        update="asynchronous",
    )
    assert bn.reachability(
        initial,
        {"A": 0, "B": 0},
        update="asynchronous",
    )
    assert bn.reachability(
        initial,
        {"A": 1, "B": 1},
        update="asynchronous",
    )
    assert not bn.reachability(
        initial,
        {"A": 1, "B": 0},
        update="asynchronous",
    )


def test_boolean_network_reachability_synchronous_follows_unique_trajectory():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})
    initial = {"A": 0, "B": 0}

    assert bn.reachability(initial, {"A": 1, "B": 0}, update="synchronous")
    assert bn.reachability(initial, {"A": 0, "B": 1}, update="synchronous")
    assert not bn.reachability(initial, {"A": 1, "B": 1}, update="synchronous")
    assert bn.reachability(initial, {"A": 1, "B": 1}, update="asynchronous")


def test_boolean_network_reachability_synchronous_quantifies_partial_states():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    assert bn.reachability(
        {},
        {"A": 1},
        update="synchronous",
        quantifier="exists",
    )
    assert not bn.reachability(
        {},
        {"A": 1},
        update="synchronous",
        quantifier="robust",
    )


def test_boolean_network_reachability_synchronous_avoids_bdd_for_concrete_state(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})

    def unexpected_bdd_import():
        raise AssertionError("the concrete synchronous trajectory should be direct")

    monkeypatch.setattr(_dynamics, "_import_dd_backend", unexpected_bdd_import)

    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
        update="synchronous",
    )

    oscillator = bt.logic.bn.BooleanNetwork({"A": "~A"})
    assert oscillator.reachability(
        {"A": 0},
        {},
        update="synchronous",
        quantifier="universal",
    )


def test_boolean_network_reachability_general_updates_multiple_components():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"})
    initial = {"A": 0, "B": 0}
    target = {"A": 1, "B": 1}

    assert bn.reachability(initial, target, update="general")
    assert not bn.reachability(initial, target, update="asynchronous")


def test_boolean_network_reachability_general_quantifies_partial_states():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"})

    assert bn.reachability(
        {},
        {"A": 1, "B": 1},
        update="general",
        quantifier="exists",
    )
    assert not bn.reachability(
        {},
        {"A": 1, "B": 1},
        update="general",
        quantifier="robust",
    )


def test_general_bdd_transition_system_reuses_compiled_partitions(monkeypatch):
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"})
    original = _dynamics._bdd_general_transition_partitions
    calls = 0

    def counted_transition_partitions(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_bdd_general_transition_partitions",
        counted_transition_partitions,
    )

    transition_system = _dynamics._BDDTransitionSystem(
        bn,
        update="general",
    )

    assert transition_system.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 0},
        quantifier="robust",
    )
    assert transition_system.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
        quantifier="robust",
    )
    assert calls == 1


def test_synchronous_bdd_transition_system_reuses_compiled_partitions(monkeypatch):
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "~B", "B": "~A"})
    original = _dynamics._bdd_synchronous_transition_partitions
    calls = 0

    def counted_transition_partitions(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_bdd_synchronous_transition_partitions",
        counted_transition_partitions,
    )

    transition_system = _dynamics._BDDTransitionSystem(
        bn,
        update="synchronous",
    )

    assert transition_system.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
        quantifier="robust",
    )
    assert transition_system.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 0},
        quantifier="robust",
    )
    assert calls == 1


def test_bdd_transition_system_reuses_compiled_partitions(monkeypatch):
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})
    original = _dynamics._bdd_asynchronous_transition_partitions
    calls = 0

    def counted_transition_partitions(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        _dynamics,
        "_bdd_asynchronous_transition_partitions",
        counted_transition_partitions,
    )

    transition_system = _dynamics._BDDTransitionSystem(
        bn,
        update="asynchronous",
    )

    assert transition_system.reachability(
        {"A": 0, "B": 1},
        {"A": 0, "B": 0},
        quantifier="robust",
    )
    assert transition_system.reachability(
        {"A": 0, "B": 1},
        {"A": 1, "B": 1},
        quantifier="robust",
    )
    assert calls == 1


def test_boolean_network_reachability_asynchronous_quantifies_partial_states():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    assert bn.reachability(
        {},
        {"A": 1},
        update="asynchronous",
        quantifier="exists",
    )
    assert not bn.reachability(
        {},
        {"A": 1},
        update="asynchronous",
        quantifier="robust",
    )


@pytest.mark.parametrize(
    ("update", "initial_state"),
    [
        ("asynchronous", {"A": 0}),
        ("general", {"A": 0}),
        ("synchronous", {"A": "*"}),
    ],
)
def test_boolean_network_reachability_reports_missing_dd(
    monkeypatch,
    update,
    initial_state,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    def missing_dd(name: str) -> Any:
        if name in {"dd.autoref", "dd.cudd"}:
            raise ImportError("missing dd")
        return __import__(name)

    monkeypatch.setattr(_dynamics, "import_module", missing_dd)

    with pytest.raises(ImportError, match=r"bonesistools\[bdd\]"):
        bn.reachability(
            initial_state,
            {"A": 1},
            update=cast(Any, update),
        )


def test_boolean_network_transition_reports_missing_dd_for_partial_states(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    def missing_dd(name: str) -> Any:
        if name in {"dd.autoref", "dd.cudd"}:
            raise ImportError("missing dd")
        return __import__(name)

    monkeypatch.setattr(_dynamics, "import_module", missing_dd)

    with pytest.raises(ImportError, match=r"bonesistools\[bdd\]"):
        bn.transition(
            {},
            {"A": 1},
            update="asynchronous",
        )


@pytest.mark.parametrize(
    "update",
    ["synchronous", "asynchronous", "general", "most-permissive"],
)
def test_boolean_network_reachability_universal_quantification(update):
    if update != "most-permissive":
        pytest.importorskip("dd.autoref")

    identity = bt.logic.bn.BooleanNetwork({"A": "A"})
    assert identity.reachability(
        {},
        {},
        update=cast(Any, update),
        quantifier="robust",
    )
    assert not identity.reachability(
        {},
        {},
        update=cast(Any, update),
        quantifier="universal",
    )

    oscillator = bt.logic.bn.BooleanNetwork({"A": "~A"})
    assert oscillator.reachability(
        {},
        {},
        update=cast(Any, update),
        quantifier="universal",
    )


def test_boolean_network_reachability_validate_options_and_states():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachable_configurations(
            {"A": 0, "B": 0},
            update=cast(Any, "parallel"),
        )

    with pytest.raises(ValueError, match="initial_state must define fixed values"):
        list(bn.reachable_configurations({"A": 0}))

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachability(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            update=cast(Any, "parallel"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'quantifier'"):
        bn.reachability(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            quantifier=cast(Any, "all"),
        )

    with pytest.raises(ValueError, match="unknown components in initial_state"):
        bn.reachability({"A": 0, "B": 0, "C": 0}, {"A": 1, "B": 0})

    with pytest.raises(ValueError, match="unknown components in target_state"):
        bn.reachability({"A": 0, "B": 0}, {"A": 1, "B": 0, "C": 0})


def test_boolean_network_transition_validates_options_and_states():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.transition(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            update=cast(Any, "parallel"),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'quantifier'"):
        bn.transition(
            {"A": 0, "B": 0},
            {"A": 1, "B": 0},
            quantifier=cast(Any, "all"),
        )

    with pytest.raises(ValueError, match="unknown components in initial_state"):
        bn.transition({"A": 0, "B": 0, "C": 0}, {"A": 1, "B": 0})

    with pytest.raises(ValueError, match="unknown components in target_state"):
        bn.transition({"A": 0, "B": 0}, {"A": 1, "B": 0, "C": 0})


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

    assert "backend" not in signature(bn.reachable_attractors).parameters

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        bn.reachable_attractors({"A": 0}, update=cast(Any, "parallel"))

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
