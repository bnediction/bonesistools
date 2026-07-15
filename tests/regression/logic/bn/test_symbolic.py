#!/usr/bin/env python

from inspect import signature
from itertools import combinations, product
from random import Random
from typing import Any, cast

import networkx as nx
import pytest

import bonesistools as bt


def _configuration_key(configuration, components):
    return tuple(configuration[component] for component in components)


def _configuration_keys(configurations, components):
    return {
        _configuration_key(configuration, components)
        for configuration in configurations
    }


def _all_configurations(components):
    return tuple(
        dict(zip(components, values))
        for values in product((0, 1), repeat=len(components))
    )


def _successors(network, configuration, update):
    next_configuration = network.next_configuration(configuration)
    unstable = tuple(
        component
        for component in network
        if configuration[component] != next_configuration[component]
    )

    if update == "synchronous":
        return (next_configuration,)

    updated_subsets = (
        ((component,) for component in unstable)
        if update == "asynchronous"
        else (
            subset
            for size in range(1, len(unstable) + 1)
            for subset in combinations(unstable, size)
        )
    )
    successors = []
    for updated in updated_subsets:
        successor = dict(configuration)
        for component in updated:
            successor[component] = next_configuration[component]
        successors.append(successor)
    return tuple(successors)


def _reachable_keys(network, initial_keys, update, within_keys):
    components = tuple(network)
    scheduled = set(initial_keys) & set(within_keys)
    pending = list(scheduled)

    while pending:
        configuration_key = pending.pop()
        configuration = dict(zip(components, configuration_key))
        for successor in _successors(network, configuration, update):
            successor_key = _configuration_key(successor, components)
            if successor_key in within_keys and successor_key not in scheduled:
                scheduled.add(successor_key)
                pending.append(successor_key)

    return scheduled


def _coreachable_keys(network, target_keys, update, within_keys):
    components = tuple(network)
    graph = nx.DiGraph()
    graph.add_nodes_from(within_keys)
    for source_key in within_keys:
        source = dict(zip(components, source_key))
        graph.add_edges_from(
            (source_key, target_key)
            for target_key in (
                _configuration_key(successor, components)
                for successor in _successors(network, source, update)
            )
            if target_key in within_keys
        )

    reachable = set(target_keys) & set(within_keys)
    for target_key in tuple(reachable):
        reachable.update(nx.ancestors(graph, target_key))
    return reachable


def _terminal_scc_keys(network, region_keys, update):
    components = tuple(network)
    graph = nx.DiGraph()
    graph.add_nodes_from(region_keys)
    for source_key in region_keys:
        source = dict(zip(components, source_key))
        graph.add_edges_from(
            (source_key, target_key)
            for target_key in (
                _configuration_key(successor, components)
                for successor in _successors(network, source, update)
            )
            if target_key in region_keys
        )

    terminal = []
    for component in nx.strongly_connected_components(graph):
        if not any(
            target not in component
            for source in component
            for target in graph.successors(source)
        ):
            terminal.append(frozenset(component))
    return frozenset(terminal)


def _symbolic_scc_keys(symbolic_sccs, components):
    return frozenset(
        frozenset(_configuration_keys(component, components))
        for component in symbolic_sccs
    )


def test_symbolic_classes_are_public_and_use_asynchronous_updates_by_default():
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})
    system = bn.symbolic()

    assert isinstance(system, bt.logic.bn.SymbolicTransitionSystem)
    assert isinstance(system.universe(), bt.logic.bn.SymbolicConfigurationSet)
    assert system.update == "asynchronous"
    assert signature(bn.symbolic).parameters["update"].default == "asynchronous"
    assert (
        signature(bt.logic.bn.SymbolicTransitionSystem).parameters["update"].default
        == "asynchronous"
    )


def test_symbolic_configuration_set_constructor_is_managed_by_system():
    with pytest.raises(TypeError, match="must be created"):
        bt.logic.bn.SymbolicConfigurationSet()


def test_symbolic_system_constructs_empty_universe_and_hypercube_union():
    system = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"}).symbolic()

    assert system.empty().is_empty()
    assert system.empty().count() == 0
    assert system.universe().is_universe()
    assert system.universe().count() == 4

    states = system.configurations(({"A": 0, "B": 0}, {"A": 1, "B": 1}))
    assert _configuration_keys(states, system.components) == {(0, 0), (1, 1)}


def test_symbolic_system_encodes_configuration_set_without_enumerating_it(
    monkeypatch,
):
    system = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"}).symbolic()
    materialized = bt.logic.ba.ConfigurationSet(
        ("B", "A"),
        ({"A": 0},),
    )

    def fail_iteration(self):
        raise AssertionError("concrete configurations were enumerated")

    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "__iter__", fail_iteration)
    states = system.configurations(materialized)

    assert states.count() == 2
    assert states.contains({"A": 0})


def test_symbolic_system_rejects_mismatched_configuration_set_components():
    system = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"}).symbolic()
    configurations = bt.logic.ba.ConfigurationSet(("A", "C"), ({"A": 0},))

    with pytest.raises(ValueError, match="components do not match"):
        system.configurations(configurations)


def test_symbolic_system_supports_empty_network_and_atypical_component_names():
    empty_system = bt.logic.bn.BooleanNetwork({}).symbolic()
    assert empty_system.universe().enumerate() == ({},)
    assert empty_system.universe().count() == 1

    system = bt.logic.bn.BooleanNetwork({"A-B": 0, "x,1": 1}).symbolic()
    states = system.configurations({"A-B": 1})
    assert states.count() == 2
    assert states.contains({"A-B": 1})


def test_symbolic_configuration_set_algebra_matches_python_sets():
    components = ("A", "B", "C")
    system = bt.logic.bn.BooleanNetwork(
        {component: component for component in components}
    ).symbolic()
    configurations = _all_configurations(components)
    universe = {
        _configuration_key(configuration, components)
        for configuration in configurations
    }
    left_keys = {key for key in universe if sum(key) % 2 == 0}
    right_keys = {key for key in universe if key[1] == 1}
    left = system.configurations(dict(zip(components, key)) for key in left_keys)
    right = system.configurations(dict(zip(components, key)) for key in right_keys)

    assert _configuration_keys(left | right, components) == left_keys | right_keys
    assert _configuration_keys(left & right, components) == left_keys & right_keys
    assert _configuration_keys(left - right, components) == left_keys - right_keys
    assert _configuration_keys(left ^ right, components) == left_keys ^ right_keys
    assert _configuration_keys(~left, components) == universe - left_keys

    assert left == system.configurations(
        dict(zip(components, key)) for key in left_keys
    )
    assert left <= system.universe()
    assert left < system.universe()
    assert system.universe() >= left
    assert system.universe() > left
    assert not left < left


def test_symbolic_configuration_set_inspection_and_pick():
    system = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"}).symbolic()
    states = system.configurations({"A": 0})

    assert states.contains({"A": 0, "B": 1})
    assert states.contains({"A": 0})
    assert not states.contains({"B": 0})
    assert {"A": 0, "B": 1} in states
    assert {"A": 0} in states
    assert {"B": 0} not in states
    assert {"C": 0} not in states
    assert {"A": 2} not in states
    assert object() not in states
    assert states.intersects({"B": 1})
    assert not states.intersects({"A": 1})
    assert states.pick() in states.enumerate()

    with pytest.raises(ValueError, match="unknown components"):
        states.contains({"C": 0})
    with pytest.raises(ValueError):
        states.contains({"A": 2})

    with pytest.raises(ValueError, match="empty symbolic set"):
        system.empty().pick()


def test_symbolic_configuration_set_conversion_uses_direct_encoded_path(
    monkeypatch,
):
    system = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B", "C": "C"}).symbolic()
    states = system.configurations(({"A": 0, "B": 0}, {"A": 1, "C": 1}))

    def fail(*args, **kwargs):
        raise AssertionError("incremental ConfigurationSet construction was used")

    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "add", fail)
    monkeypatch.setattr(bt.logic.ba.ConfigurationSet, "compress", fail)
    materialized = states.configurations()

    assert materialized.count() == states.count()
    assert _configuration_keys(materialized, system.components) == _configuration_keys(
        states,
        system.components,
    )


def test_symbolic_sets_from_different_systems_are_incompatible():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})
    first_system = bn.symbolic(update="asynchronous")
    second_system = bn.symbolic(update="asynchronous")
    first = first_system.universe()
    second = second_system.universe()

    assert first != second
    with pytest.raises(ValueError, match="different transition systems"):
        _ = first | second
    with pytest.raises(ValueError, match="different transition systems"):
        _ = first <= second
    with pytest.raises(ValueError, match="different transition systems"):
        first.intersects(second)
    with pytest.raises(ValueError, match="different transition systems"):
        first_system.post(second)
    with pytest.raises(ValueError, match="different transition systems"):
        first.reachable(within=second)
    with pytest.raises(TypeError, match="cannot be implicitly converted"):
        second_system.configurations(first)


def test_symbolic_set_keeps_transition_system_alive():
    system = bt.logic.bn.BooleanNetwork({"A": "~A"}).symbolic()
    states = system.universe()

    del system

    assert states.count() == 2
    assert states.system.update == "asynchronous"


def test_symbolic_system_captures_network_at_compilation():
    network = bt.logic.bn.BooleanNetwork({"A": 0})
    system = network.symbolic(update="synchronous")
    network["A"] = 1
    exposed = system.to_boolean_network()
    exposed["A"] = 1

    successors = system.configurations({"A": 0}).post()

    assert successors.enumerate() == ({"A": 0},)
    assert system.to_boolean_network().rules == {"A": "0"}
    assert system.to_boolean_network() is not exposed


def test_symbolic_configuration_set_repr_is_compact():
    states = (
        bt.logic.bn.BooleanNetwork({"A": "A", "B": "B", "C": "C"}).symbolic().universe()
    )

    assert repr(states) == (
        "SymbolicConfigurationSet(components=3, update='asynchronous')"
    )


@pytest.mark.parametrize("update", ["asynchronous", "general"])
def test_non_reflexive_symbolic_updates_give_fixed_points_no_direct_neighbors(
    update,
):
    system = bt.logic.bn.BooleanNetwork({"A": "A"}).symbolic(update=update)
    fixed = system.configurations({"A": 1})

    assert fixed.post().is_empty()
    assert fixed.pre().is_empty()


def test_synchronous_symbolic_updates_keep_fixed_point_self_loops():
    system = bt.logic.bn.BooleanNetwork({"A": "A"}).symbolic(update="synchronous")
    fixed = system.configurations({"A": 1})

    assert fixed.post() == fixed
    assert fixed.pre() == fixed


def test_symbolic_terminal_sccs_are_relative_to_the_induced_region():
    system = bt.logic.bn.BooleanNetwork({"A": 1}).symbolic(update="general")
    transient_region = system.configurations({"A": 0})

    relative = transient_region.terminal_sccs()
    global_sccs = system.terminal_sccs()

    assert tuple(component.enumerate() for component in relative) == (({"A": 0},),)
    assert tuple(component.enumerate() for component in global_sccs) == (({"A": 1},),)


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_symbolic_dynamics_match_explicit_state_graph(update):
    network = bt.logic.bn.BooleanNetwork(
        {
            "A": "~B",
            "B": "A | C",
            "C": "~A",
        }
    )
    system = network.symbolic(update=update)
    components = system.components
    configurations = _all_configurations(components)
    universe_keys = {
        _configuration_key(configuration, components)
        for configuration in configurations
    }
    initial_keys = {(0, 0, 0), (0, 0, 1)}
    target_keys = {(1, 1, 0)}
    region_keys = {key for key in universe_keys if key[1] == 0 or key[2] == 0}
    initial = system.configurations(dict(zip(components, key)) for key in initial_keys)
    target = system.configurations(dict(zip(components, key)) for key in target_keys)
    region = system.configurations(dict(zip(components, key)) for key in region_keys)

    expected_post = set()
    for initial_key in initial_keys:
        configuration = dict(zip(components, initial_key))
        expected_post.update(
            _configuration_key(successor, components)
            for successor in _successors(network, configuration, update)
        )
    assert _configuration_keys(initial.post(), components) == expected_post

    expected_pre = {
        source_key
        for source_key in universe_keys
        if any(
            _configuration_key(successor, components) in target_keys
            for successor in _successors(
                network,
                dict(zip(components, source_key)),
                update,
            )
        )
    }
    assert _configuration_keys(target.pre(), components) == expected_pre

    assert _configuration_keys(
        initial.reachable(within=region),
        components,
    ) == _reachable_keys(network, initial_keys, update, region_keys)
    assert _configuration_keys(
        target.coreachable(within=region),
        components,
    ) == _coreachable_keys(network, target_keys, update, region_keys)

    assert _symbolic_scc_keys(
        system.terminal_sccs(),
        components,
    ) == _terminal_scc_keys(
        network,
        universe_keys,
        update,
    )
    assert _symbolic_scc_keys(
        region.terminal_sccs(),
        components,
    ) == _terminal_scc_keys(
        network,
        region_keys,
        update,
    )


def test_symbolic_dynamics_match_random_explicit_state_graphs():
    rng = Random(0)
    components = ("A", "B", "C")
    rules = (
        "0",
        "1",
        "A",
        "~A",
        "B",
        "~B",
        "C",
        "~C",
        "A & B",
        "A | C",
        "(A & ~B) | (~A & B)",
    )
    universe_keys = set(product((0, 1), repeat=len(components)))

    for _ in range(10):
        network = bt.logic.bn.BooleanNetwork(
            {component: rng.choice(rules) for component in components}
        )
        initial_keys = set(rng.sample(tuple(universe_keys), 2))
        target_keys = set(rng.sample(tuple(universe_keys), 2))
        region_keys = set(rng.sample(tuple(universe_keys), 6))

        for update in ("synchronous", "asynchronous", "general"):
            system = network.symbolic(update=update)
            initial = system.configurations(
                dict(zip(components, key)) for key in initial_keys
            )
            target = system.configurations(
                dict(zip(components, key)) for key in target_keys
            )
            region = system.configurations(
                dict(zip(components, key)) for key in region_keys
            )

            assert _configuration_keys(
                initial.reachable(within=region),
                components,
            ) == _reachable_keys(network, initial_keys, update, region_keys)
            assert _configuration_keys(
                target.coreachable(within=region),
                components,
            ) == _coreachable_keys(network, target_keys, update, region_keys)
            assert _symbolic_scc_keys(
                system.terminal_sccs(),
                components,
            ) == _terminal_scc_keys(network, universe_keys, update)


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_small_symbolic_universe_uses_explicit_terminal_sccs(update, monkeypatch):
    network = bt.logic.bn.BooleanNetwork({"A": "~A", "B": "A"})
    system = network.symbolic(update=update)
    original = type(system)._explicit_terminal_scc_nodes
    calls = 0

    def tracked(self, states):
        nonlocal calls
        calls += 1
        return original(self, states)

    monkeypatch.setattr(type(system), "_explicit_terminal_scc_nodes", tracked)
    expected = system.terminal_sccs()
    observed = system.terminal_sccs(system.universe())
    delegated = system.universe().terminal_sccs()

    assert observed == expected
    assert delegated == expected
    assert calls == 3


@pytest.mark.parametrize(
    ("update", "algorithm"),
    [
        ("synchronous", "_synchronous_terminal_cycles"),
        ("asynchronous", "_transition_guided_reduction"),
        ("general", "_terminal_sccs"),
    ],
)
def test_large_symbolic_universe_uses_semantics_specific_terminal_sccs(
    update,
    algorithm,
    monkeypatch,
):
    network = bt.logic.bn.BooleanNetwork(
        {f"x{index}": f"~x{index}" for index in range(8)}
    )
    system = network.symbolic(update=update)
    original = getattr(type(system), algorithm)
    calls = 0

    def tracked(self, states):
        nonlocal calls
        calls += 1
        return original(self, states)

    monkeypatch.setattr(type(system), algorithm, tracked)

    system.terminal_sccs()

    assert calls == 1


@pytest.mark.parametrize(
    ("update", "specialized_algorithm"),
    [
        ("synchronous", "_synchronous_terminal_sccs"),
        ("asynchronous", "_transition_guided_reduction"),
    ],
)
def test_large_non_closed_region_uses_induced_symbolic_sccs(
    update,
    specialized_algorithm,
    monkeypatch,
):
    network = bt.logic.bn.BooleanNetwork(
        {f"x{index}": f"~x{index}" for index in range(9)}
    )
    system = network.symbolic(update=update)
    region = system.configurations({"x0": 0})
    calls = 0

    def generic(_self, _states):
        nonlocal calls
        calls += 1
        return ()

    def specialized(*_args, **_kwargs):
        raise AssertionError("non-closed regions require induced SCC semantics")

    monkeypatch.setattr(type(system), "_terminal_sccs", generic)
    monkeypatch.setattr(type(system), specialized_algorithm, specialized)

    assert region.terminal_sccs() == ()
    assert calls == 1


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_large_closed_reachable_region_matches_boolean_network_attractors(update):
    components = ("input",) + tuple(f"x{index}" for index in range(8))
    network = bt.logic.bn.BooleanNetwork(
        {
            "input": "input",
            **{component: f"~{component}" for component in components[1:]},
        }
    )
    system = network.symbolic(update=update)
    reachable = system.configurations({"input": 0}).reachable()

    expected = network.attractors({"input": 0}, update=update)
    observed = reachable.terminal_sccs()

    assert reachable.count() == 256
    assert _symbolic_scc_keys(observed, components) == _symbolic_scc_keys(
        expected,
        components,
    )


def test_non_unate_synchronous_system_uses_symbolic_terminal_scc_search(
    monkeypatch,
):
    components = tuple(f"x{index}" for index in range(20))
    rule = " | ".join(
        f"({left} & {right}) | (~{left} & ~{right})"
        for left, right in zip(components[::2], components[1::2])
    )
    system = bt.logic.bn.BooleanNetwork(
        {component: rule for component in components}
    ).symbolic(update="synchronous")
    cycle_calls = 0
    symbolic_calls = 0

    def cycles(_self, _states):
        nonlocal cycle_calls
        cycle_calls += 1
        return ()

    def symbolic(_self, _states):
        nonlocal symbolic_calls
        symbolic_calls += 1
        return ()

    monkeypatch.setattr(type(system), "_synchronous_terminal_cycles", cycles)
    monkeypatch.setattr(type(system), "_terminal_sccs", symbolic)

    assert system.terminal_sccs() == ()
    assert cycle_calls == 0
    assert symbolic_calls == 1


def test_large_synchronous_system_uses_symbolic_terminal_scc_search(monkeypatch):
    components = tuple(f"x{index}" for index in range(90))
    system = bt.logic.bn.BooleanNetwork(
        {component: component for component in components}
    ).symbolic(update="synchronous")
    cycle_calls = 0
    symbolic_calls = 0

    def cycles(_self, _states):
        nonlocal cycle_calls
        cycle_calls += 1
        return ()

    def symbolic(_self, _states):
        nonlocal symbolic_calls
        symbolic_calls += 1
        return ()

    monkeypatch.setattr(type(system), "_synchronous_terminal_cycles", cycles)
    monkeypatch.setattr(type(system), "_terminal_sccs", symbolic)

    assert system.terminal_sccs() == ()
    assert cycle_calls == 0
    assert symbolic_calls == 1


def test_symbolic_api_does_not_expose_backend_objects_or_instance_dicts():
    system = bt.logic.bn.BooleanNetwork({"A": "A"}).symbolic()
    states = system.universe()

    assert not hasattr(system, "__dict__")
    assert not hasattr(states, "__dict__")
    assert not hasattr(system, "bdd")
    assert not hasattr(states, "node")


def test_symbolic_system_rejects_most_permissive_updates():
    network = bt.logic.bn.BooleanNetwork({"A": "A"})

    with pytest.raises(ValueError, match="invalid argument value for 'update'"):
        network.symbolic(update=cast(Any, "most-permissive"))
