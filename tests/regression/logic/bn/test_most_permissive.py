#!/usr/bin/env python

from typing import Any, cast

import pytest

import bonesistools as bt
from bonesistools.logic.boolean_network import _most_permissive


def _configuration_from_bits(bits):
    return {f"x{index}": int(bit) for index, bit in enumerate(bits, start=1)}


def _configuration_bits(configuration):
    return "".join(str(configuration[f"x{index}"]) for index in range(1, 4))


def _canonical_attractors(attractors):
    return frozenset(
        frozenset(tuple(sorted(configuration.items())) for configuration in attractor)
        for attractor in attractors
    )


def _assert_most_permissive_stg(bn, expected_targets):
    states = tuple(expected_targets)

    for source in states:
        enumerated = [
            _configuration_bits(configuration)
            for configuration in bn.reachable_configurations(
                _configuration_from_bits(source),
                update="most-permissive",
            )
        ]
        enumerated_targets = {target for target in enumerated if target != source}

        assert len(enumerated) == len(set(enumerated))
        assert set(enumerated) == expected_targets[source] | {source}
        assert enumerated_targets == expected_targets[source]

        for target in states:
            if target == source:
                continue

            assert bn.reachability(
                _configuration_from_bits(source),
                _configuration_from_bits(target),
                update="most-permissive",
                quantifier="exists",
            ) is (target in expected_targets[source])


def test_boolean_network_trap_spaces_minimal_fixed_points():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.trap_spaces() == (
        bt.logic.ba.Hypercube({"A": 0, "B": 0}),
        bt.logic.ba.Hypercube({"A": 1, "B": 1}),
    )


def test_boolean_network_trap_spaces_universal_when_no_fixed_subspace():
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})

    assert bn.trap_spaces() == (bt.logic.ba.Hypercube({}),)


def test_boolean_network_trap_spaces_match_unrestricted_mp_attractors():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "A",
            "B": "A | ~B",
        }
    )

    trap_spaces = bn.trap_spaces()
    attractors = bn.attractors(update="most-permissive")

    assert tuple(attractor.hypercubes() for attractor in attractors) == tuple(
        (trap_space,) for trap_space in trap_spaces
    )


def test_boolean_network_trap_spaces_rejects_invalid_kind():
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    with pytest.raises(ValueError, match="invalid argument value for 'kind'"):
        bn.trap_spaces(kind=cast(Any, "maximal"))


def test_boolean_network_trap_spaces_rejects_open_network():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="undefined components"):
        bn.trap_spaces()


def test_boolean_network_principal_trap_space_preserves_fixed_point():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.principal_trap_space({"A": 0, "B": 0}) == bt.logic.ba.Hypercube(
        {"A": 0, "B": 0}
    )


def test_boolean_network_principal_trap_space_relaxes_unstable_components():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": 1,
            "x2": "x1",
            "x3": "(~x1 & x2) | x3",
        }
    )

    assert bn.principal_trap_space(
        {"x1": 0, "x2": 0, "x3": 1}
    ) == bt.logic.ba.Hypercube({"x3": 1})


def test_boolean_network_principal_trap_space_rejects_invalid_configuration():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    with pytest.raises(ValueError, match="configuration must define fixed values"):
        bn.principal_trap_space({"A": 0})


def test_boolean_network_principal_trap_space_rejects_open_network():
    bn = bt.logic.bn.BooleanNetwork({"A": "B"}, check=False)

    with pytest.raises(ValueError, match="undefined components"):
        bn.principal_trap_space({"A": 0})


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
    initial_configuration = {"x1": 0, "x2": 0, "x3": 1}

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_configuration,
        relaxed_components=("x1",),
    ) == bt.logic.ba.Hypercube({"x2": 0, "x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_configuration,
        relaxed_components=("x1", "x2"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_configuration,
        relaxed_components=("x1", "x2", "x3"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    initial_configuration = {"x1": 0, "x2": 1, "x3": 1}

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_configuration,
        relaxed_components=("x1", "x2", "x3"),
    ) == bt.logic.ba.Hypercube({"x3": 1})

    assert _most_permissive._smallest_closed_hypercube(
        bn,
        initial_configuration,
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
    initial_configuration = {"x1": 0, "x2": 0, "x3": 1}

    assert bn.reachability(
        initial_configuration,
        {"x1": 1, "x2": 0, "x3": 1},
        update="most-permissive",
    )
    assert bn.reachability(
        initial_configuration,
        {"x1": 1, "x2": 1, "x3": 1},
        update="most-permissive",
    )
    assert not bn.reachability(
        initial_configuration,
        {"x1": 0, "x2": 0, "x3": 0},
        update="most-permissive",
    )
    assert not bn.reachability(
        initial_configuration,
        {"x1": 0, "x2": 1, "x3": 1},
        update="most-permissive",
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


def test_boolean_network_reachable_configurations_iterates_regions_lazily(
    monkeypatch,
):
    bn = bt.logic.bn.BooleanNetwork({"A": "~A"})
    initial_configuration = bt.logic.ba.Hypercube({"A": 0})
    first_region = _most_permissive._MPTransitionRegion(
        closure_components=("A",),
        free_components=(),
        closed_hypercube=initial_configuration,
        irreversible_components=(),
    )

    def guarded_regions(*args, **kwargs):
        yield first_region
        raise AssertionError("the next region was evaluated eagerly")

    monkeypatch.setattr(
        _most_permissive,
        "_iter_most_permissive_transition_regions",
        guarded_regions,
    )

    configurations = bn.reachable_configurations(
        {"A": 0},
        update="most-permissive",
    )

    assert next(configurations) == {"A": 0}


def test_boolean_network_reachable_configurations_iterates_combinations_lazily(
    monkeypatch,
):
    components = tuple(f"x{index}" for index in range(12))
    bn = bt.logic.bn.BooleanNetwork(
        {component: f"~{component}" for component in components}
    )
    initial_configuration = {component: 0 for component in components}
    space = _most_permissive._most_permissive_transition_space(
        bn, initial_configuration
    )
    original_combinations = _most_permissive.combinations

    def guarded_combinations(values, length):
        if length > 0:
            raise AssertionError("later combinations were evaluated eagerly")
        return original_combinations(values, length)

    monkeypatch.setattr(
        _most_permissive,
        "combinations",
        guarded_combinations,
    )

    assert next(iter(space)) == initial_configuration


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


def test_boolean_network_attractors_most_permissive():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": 1,
            "B": "A",
            "C": "(~A & B) | C",
        }
    )

    attractors = bn.attractors(
        {"A": 0, "B": 0, "C": 0},
        update="most-permissive",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"A": 1, "B": 1, "C": 0},),
        ({"A": 1, "B": 1, "C": 1},),
    )


def test_boolean_network_attractors_most_permissive_partial_state():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "AA": "AA",
            "BB": "AA | BB",
            "CC": "~AA | CC",
        }
    )

    attractors = bn.attractors(
        {"AA": 0},
        update="most-permissive",
    )

    assert tuple(attractor.enumerate() for attractor in attractors) == (
        ({"AA": 0, "BB": 0, "CC": 1},),
        ({"AA": 0, "BB": 1, "CC": 1},),
    )


def test_boolean_network_attractors_evaluates_non_unate_rules_exactly():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "~A",
            "B": "~B",
            "always": "A | ~A",
            "never": "A & ~A",
            "xor": "(A & ~B) | (~A & B)",
        }
    )

    attractors = bn.attractors(update="most-permissive")

    assert len(attractors) == 1
    assert attractors[0].hypercubes() == (
        bt.logic.ba.Hypercube({"always": 1, "never": 0}),
    )


def test_boolean_network_attractors_preserves_incomparable_minima():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "A",
            "B": "A | ~B",
        }
    )

    attractors = bn.attractors(update="most-permissive")

    assert tuple(attractor.hypercubes() for attractor in attractors) == (
        (bt.logic.ba.Hypercube({"A": 0}),),
        (bt.logic.ba.Hypercube({"A": 1, "B": 1}),),
    )


def test_boolean_network_attractors_most_permissive_defaults_to_free_state():
    bn = bt.logic.bn.BooleanNetwork({"A": "A", "B": "B"})

    omitted = bn.attractors(update="most-permissive")
    explicit = bn.attractors({}, update="most-permissive")

    assert tuple(attractor.enumerate() for attractor in omitted) == tuple(
        attractor.enumerate() for attractor in explicit
    )


def test_boolean_network_reachability_most_permissive_tracks_irreversibles():
    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": "A"})

    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 0},
        update="most-permissive",
    )
    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 0},
        update="most-permissive",
    )
    assert bn.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
        update="most-permissive",
    )
    assert not bn.reachability(
        {"A": 0, "B": 0},
        {"A": 0, "B": 1},
        update="most-permissive",
    )


def test_boolean_network_reachability_most_permissive_allows_reversible_space():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 0, "B": 0},
        update="most-permissive",
    )
    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 1, "B": 1},
        update="most-permissive",
    )
    assert bn.reachability(
        {"A": 0, "B": 1},
        {"A": 1, "B": 0},
        update="most-permissive",
    )


def test_boolean_network_reachability_most_permissive_rejects_self_support():
    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": "A"})

    assert not bn.reachability(
        {"A": 0, "B": 0},
        {"A": 1, "B": 1},
        update="most-permissive",
        quantifier="exists",
    )


def test_boolean_network_reachability_accepts_partial_states_existentially():
    bn = bt.logic.bn.BooleanNetwork({"A": "A"})

    assert bn.reachability(
        {},
        {"A": 1},
        update="most-permissive",
        quantifier="exists",
    )
    assert not bn.reachability(
        {},
        {"A": 1},
        update="most-permissive",
        quantifier="robust",
    )


def test_boolean_network_reachability_defaults_to_robust_quantification():
    bn = bt.logic.bn.BooleanNetwork({"A": 1})

    assert bn.reachability({}, {"A": 1}, update="most-permissive")


def test_boolean_network_reachability_evaluates_non_unate_rules_exactly():
    bn = bt.logic.bn.BooleanNetwork(
        {
            "A": "~A",
            "always": "A | ~A",
            "never": "A & ~A",
        }
    )
    initial_configuration = {"A": 0, "always": 1, "never": 0}

    assert bn.reachability(
        initial_configuration,
        {"A": 1, "always": 1, "never": 0},
        update="most-permissive",
        quantifier="exists",
    )
    assert not bn.reachability(
        initial_configuration,
        {"always": 0},
        update="most-permissive",
        quantifier="exists",
    )
    assert not bn.reachability(
        initial_configuration,
        {"never": 1},
        update="most-permissive",
        quantifier="exists",
    )


@pytest.mark.parametrize("quantifier", ["exists", "robust", "universal"])
def test_boolean_network_transition_equals_most_permissive_reachability(
    quantifier,
):
    bn = bt.logic.bn.BooleanNetwork({"A": 1, "B": "A"})

    for initial_configuration, target_configuration in (
        ({"A": 0, "B": 0}, {"A": 1, "B": 1}),
        ({"A": 0, "B": 0}, {"A": 0, "B": 1}),
        ({}, {"A": 1}),
    ):
        transition = bn.transition(
            initial_configuration,
            target_configuration,
            update="most-permissive",
            quantifier=cast(Any, quantifier),
        )
        reachability = bn.reachability(
            initial_configuration,
            target_configuration,
            update="most-permissive",
            quantifier=cast(Any, quantifier),
        )

        assert transition is reachability


def test_i3_feed_forward_loop_is_reachable_only_with_most_permissive_dynamics():
    pytest.importorskip("dd.autoref")
    bn = bt.logic.bn.BooleanNetwork(
        {
            "x1": 1,
            "x2": "x1",
            "x3": "~x1 & x2",
        }
    )
    initial_configuration = {"x1": 0, "x2": 0, "x3": 0}
    target_configuration = {"x1": 1, "x2": 1, "x3": 1}

    for update in ("asynchronous", "synchronous", "general"):
        assert not bn.transition(
            initial_configuration,
            target_configuration,
            update=cast(Any, update),
        )
        assert not bn.reachability(
            initial_configuration,
            target_configuration,
            update=cast(Any, update),
        )

    assert bn.transition(
        initial_configuration,
        target_configuration,
        update="most-permissive",
    )
    assert bn.reachability(
        initial_configuration,
        target_configuration,
        update="most-permissive",
    )
    assert target_configuration in tuple(
        bn.reachable_configurations(
            initial_configuration,
            update="most-permissive",
        )
    )
