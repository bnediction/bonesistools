#!/usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from itertools import product
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Mapping,
    Set,
    Tuple,
    Union,
)

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import Configuration, HypercubeLike
from ._asynchronous import (
    _asynchronous_successor_configuration_bits,
    _tarjan_asynchronous_terminal_sccs,
)
from ._bdd import _bdd_equivalence, _bdd_rule
from ._general import (
    _general_successor_configuration_bits,
)
from ._synchronous import (
    _synchronous_successor_configuration_bits,
)

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import _BDDConfigurationSetNode, _BDDManager

_ConfigurationBits = int
_CompiledRule = Callable[[_ConfigurationBits], int]
_FIXED_POINTS_BITSET_COMPONENT_LIMIT = 12


@dataclass
class _ForwardReduction:
    """Forward reachability task for one transition label."""

    transition_index: int
    forward: _BDDConfigurationSetNode


@dataclass
class _ExtendedComponentReduction:
    """Backward extended-component task within a forward set."""

    transition_index: int
    forward: _BDDConfigurationSetNode
    extended_component: _BDDConfigurationSetNode


@dataclass
class _ForwardBasinReduction:
    """Backward basin task used to discard states before a forward set."""

    forward: _BDDConfigurationSetNode
    basin: _BDDConfigurationSetNode


@dataclass
class _BottomBasinReduction:
    """Backward basin task used to discard states outside a bottom set."""

    bottom: _BDDConfigurationSetNode
    basin: _BDDConfigurationSetNode


_TransitionGuidedReduction = Union[
    _ForwardReduction,
    _ExtendedComponentReduction,
    _ForwardBasinReduction,
    _BottomBasinReduction,
]


def _transition_guided_reduction_size(
    reduction: _TransitionGuidedReduction,
) -> int:
    """Return the BDD size currently driving one reduction task."""

    if isinstance(reduction, _ForwardReduction):
        return len(reduction.forward)
    if isinstance(reduction, _ExtendedComponentReduction):
        return len(reduction.extended_component)
    return len(reduction.basin)


def _restrict_transition_guided_reduction(
    reduction: _TransitionGuidedReduction,
    states: _BDDConfigurationSetNode,
) -> None:
    """Restrict all state sets held by a reduction task in place."""

    if isinstance(reduction, _ForwardReduction):
        reduction.forward &= states
    elif isinstance(reduction, _ExtendedComponentReduction):
        reduction.forward &= states
        reduction.extended_component &= states
    elif isinstance(reduction, _ForwardBasinReduction):
        reduction.forward &= states
        reduction.basin &= states
    else:
        reduction.bottom &= states
        reduction.basin &= states


def _automatic_reachable_attractors_backend(
    network: "BooleanNetwork",
    initial: HypercubeLike,
    initial_configurations: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Literal["explicit", "bdd"]:
    """Select a conservative backend from bounds on the reachable state space."""

    if update == "synchronous":
        return "explicit" if len(initial_configurations) <= 1_024 else "bdd"

    if update == "asynchronous" and len(initial_configurations) >= 2_048:
        return "bdd"

    n_free = _reachable_configuration_dimension_upper_bound(network, initial)
    explicit_dimension_limit = 13 if update == "asynchronous" else 10
    return "explicit" if n_free <= explicit_dimension_limit else "bdd"


def _automatic_reachable_configurations_backend(
    network: "BooleanNetwork",
    initial: Hypercube,
    *,
    update: Literal["asynchronous", "general"],
) -> Literal["explicit", "bdd"]:
    """Select a traversal from a conservative reachable-state upper bound."""

    n_free = _reachable_configuration_dimension_upper_bound(network, initial)
    explicit_dimension_limit = 12 if update == "asynchronous" else 8
    return "explicit" if n_free <= explicit_dimension_limit else "bdd"


def _reachable_configuration_dimension_upper_bound(
    network: "BooleanNetwork",
    initial: HypercubeLike,
) -> int:
    """Bound reachable configurations by the smallest closed hypercube."""

    from ._most_permissive import _smallest_closed_hypercube

    closed_hypercube = _smallest_closed_hypercube(
        network,
        initial,
        relaxed_components=network.keys(),
    )
    n_fixed = sum(value.is_fixed for value in closed_hypercube.values())
    return len(network) - n_fixed


def _fixed_points(
    network: "BooleanNetwork",
    *,
    limit: int,
) -> List[Configuration]:
    """Return fixed points with bitsets or a direct symbolic constraint."""

    components = tuple(network.keys())
    if len(components) <= _FIXED_POINTS_BITSET_COMPONENT_LIMIT:
        return _bitset_fixed_points(network, components, limit=limit)

    return _bdd_fixed_points(network, components, limit=limit)


def _bitset_fixed_points(
    network: "BooleanNetwork",
    components: Tuple[str, ...],
    *,
    limit: int,
) -> List[Configuration]:
    """Enumerate a small state space with compiled bitset rules."""

    compiled_rules = _compile_bitset_rules(network, components)
    fixed_points = []

    for values in product((0, 1), repeat=len(components)):
        configuration = sum(value << index for index, value in enumerate(values))
        if _next_configuration_bits(compiled_rules, configuration) != configuration:
            continue

        fixed_points.append(dict(zip(components, values)))
        if limit and len(fixed_points) >= limit:
            break

    return fixed_points


def _bdd_fixed_points(
    network: "BooleanNetwork",
    components: Tuple[str, ...],
    *,
    limit: int,
) -> List[Configuration]:
    """Solve all local fixed-point equations as one BDD constraint."""

    bdd = _import_dd_backend().BDD()
    variables = tuple(f"x{index}" for index in range(len(components)))
    bdd.declare(*variables)
    component_variables = dict(zip(components, variables))

    fixed_points = bdd.true
    for component, variable_name in zip(components, variables):
        rule = _bdd_rule(
            network,
            network[component],
            bdd=bdd,
            variables=component_variables,
        )
        fixed_points &= _bdd_equivalence(bdd.var(variable_name), rule)

    configurations = []
    for values in _iter_bdd_assignments(
        bdd,
        fixed_points,
        variables,
    ):
        configurations.append(dict(zip(components, values)))
        if limit and len(configurations) >= limit:
            break

    return configurations


def _iter_bdd_assignments(
    bdd: _BDDManager,
    node: _BDDConfigurationSetNode,
    variables: Tuple[str, ...],
) -> Iterator[Tuple[int, ...]]:
    """Yield satisfying assignments in lexicographic Boolean order."""

    pending: List[Tuple[_BDDConfigurationSetNode, int, Tuple[int, ...]]] = [
        (node, 0, ())
    ]
    while pending:
        current, index, values = pending.pop()
        if current == bdd.false:
            continue
        if index == len(variables):
            yield values
            continue

        variable = variables[index]
        high = bdd.let({variable: True}, current)
        low = bdd.let({variable: False}, current)
        pending.append((high, index + 1, values + (1,)))
        pending.append((low, index + 1, values + (0,)))


def _bdd_reachable_attractors(
    network: "BooleanNetwork",
    initial_configurations: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Compute reachable attractors through a BDD reachability set."""

    from ._symbolic import SymbolicTransitionSystem

    transition_system = SymbolicTransitionSystem._from_validated_network(
        network,
        update=update,
    )
    return transition_system._reachable_attractors(initial_configurations)


def _bdd_reachable_configurations(
    network: "BooleanNetwork",
    initial_configuration: Hypercube,
    *,
    update: Literal["asynchronous", "general"],
) -> Iterator[Configuration]:
    """Enumerate a symbolic reachable-state closure."""

    from ._symbolic import SymbolicTransitionSystem

    transition_system = SymbolicTransitionSystem._from_validated_network(
        network,
        update=update,
    )
    return transition_system._reachable_configurations(initial_configuration)


def _finite_state_transition(
    network: "BooleanNetwork",
    source: HypercubeLike,
    target: HypercubeLike,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test a finite-state transition directly or through a symbolic BDD."""

    components = tuple(network.keys())
    source_hypercube = _validate_hypercube(
        components,
        source,
        name="source",
    )
    target_hypercube = _validate_hypercube(
        components,
        target,
        name="target",
    )
    source_is_concrete = (
        source_hypercube.components == frozenset(components)
        and source_hypercube.is_fully_specified
    )
    if not source_is_concrete:
        return _bdd_transition(
            network,
            source_hypercube,
            target_hypercube,
            update=update,
            quantifier=quantifier,
        )

    source_bits = _configuration_bits_from_hypercube(source_hypercube, components)
    if update == "synchronous":
        return _synchronous_transition_from_configuration(
            network,
            source_bits,
            target_hypercube,
            components,
            quantifier=quantifier,
        )

    return _componentwise_transition_from_configuration(
        network,
        source_bits,
        target_hypercube,
        components,
        update=update,
        quantifier=quantifier,
    )


def _synchronous_transition_from_configuration(
    network: "BooleanNetwork",
    source: _ConfigurationBits,
    target: Hypercube,
    components: Tuple[str, ...],
    *,
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test a synchronous transition from one concrete configuration."""

    next_configuration = _next_configuration_bits(
        _compile_bitset_rules(network, components),
        source,
    )
    if not _configuration_bits_match_hypercube(
        next_configuration,
        target,
        components,
    ):
        return False

    if quantifier != "universal":
        return True

    return target.components == frozenset(components) and target.is_fully_specified


def _componentwise_transition_from_configuration(
    network: "BooleanNetwork",
    source: _ConfigurationBits,
    target: Hypercube,
    components: Tuple[str, ...],
    *,
    update: Literal["asynchronous", "general"],
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test an asynchronous or general transition from one configuration."""

    fixed_mask, fixed_values = _hypercube_fixed_masks(target, components)
    all_components_mask = (1 << len(components)) - 1
    free_mask = all_components_mask ^ fixed_mask
    required_changes = fixed_mask & (source ^ fixed_values)

    if update == "asynchronous":
        if quantifier == "universal" and free_mask:
            return False

        if required_changes:
            if required_changes & (required_changes - 1):
                return False
            candidates = required_changes
        else:
            candidates = free_mask

        return _any_candidate_is_unstable(
            network,
            source,
            components,
            candidates,
        )

    if quantifier == "universal":
        if not required_changes:
            return False
        candidates = required_changes | free_mask
        return _all_candidates_are_unstable(
            network,
            source,
            components,
            candidates,
        )

    if required_changes:
        return _all_candidates_are_unstable(
            network,
            source,
            components,
            required_changes,
        )

    return _any_candidate_is_unstable(
        network,
        source,
        components,
        free_mask,
    )


def _hypercube_fixed_masks(
    hypercube: Hypercube,
    components: Tuple[str, ...],
) -> Tuple[_ConfigurationBits, _ConfigurationBits]:
    """Encode the fixed dimensions and values of one hypercube."""

    fixed_mask = 0
    value_mask = 0
    for index, component in enumerate(components):
        if component not in hypercube or not hypercube[component].is_fixed:
            continue

        bit = 1 << index
        fixed_mask |= bit
        if hypercube[component].value:
            value_mask |= bit

    return fixed_mask, value_mask


def _any_candidate_is_unstable(
    network: "BooleanNetwork",
    state: _ConfigurationBits,
    components: Tuple[str, ...],
    candidates: _ConfigurationBits,
) -> bool:
    """Return whether any candidate component is unstable in one state."""

    component_indices = {component: index for index, component in enumerate(components)}
    while candidates:
        bit = candidates & -candidates
        index = bit.bit_length() - 1
        component = components[index]
        rule = _compile_bitset_rule(
            network,
            network[component],
            component_indices,
        )
        if bool(rule(state)) != bool(state & bit):
            return True
        candidates ^= bit

    return False


def _all_candidates_are_unstable(
    network: "BooleanNetwork",
    state: _ConfigurationBits,
    components: Tuple[str, ...],
    candidates: _ConfigurationBits,
) -> bool:
    """Return whether every candidate component is unstable in one state."""

    component_indices = {component: index for index, component in enumerate(components)}
    while candidates:
        bit = candidates & -candidates
        index = bit.bit_length() - 1
        component = components[index]
        rule = _compile_bitset_rule(
            network,
            network[component],
            component_indices,
        )
        if bool(rule(state)) == bool(state & bit):
            return False
        candidates ^= bit

    return True


def _bdd_transition(
    network: "BooleanNetwork",
    source: HypercubeLike,
    target: HypercubeLike,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test a one-step transition through a symbolic BDD relation."""

    from ._symbolic import SymbolicTransitionSystem

    transition_system = SymbolicTransitionSystem._from_validated_network(
        network,
        update=update,
    )
    return transition_system._transition(
        source,
        target,
        quantifier=quantifier,
    )


def _bdd_reachability(
    network: "BooleanNetwork",
    source: HypercubeLike,
    target: HypercubeLike,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test reachability through a symbolic BDD backward closure."""

    from ._symbolic import SymbolicTransitionSystem

    transition_system = SymbolicTransitionSystem._from_validated_network(
        network,
        update=update,
    )
    return transition_system._reachability(
        source,
        target,
        quantifier=quantifier,
    )


def _synchronous_reachability(
    network: "BooleanNetwork",
    source: HypercubeLike,
    target: HypercubeLike,
    *,
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Use a direct trajectory when the synchronous source is concrete."""

    components = tuple(network.keys())
    source_hypercube = _validate_hypercube(
        components,
        source,
        name="source",
    )
    if (
        source_hypercube.components != frozenset(components)
        or not source_hypercube.is_fully_specified
    ):
        return _bdd_reachability(
            network,
            source_hypercube,
            target,
            update="synchronous",
            quantifier=quantifier,
        )

    target_hypercube = _validate_hypercube(
        components,
        target,
        name="target",
    )
    compiled_rules = _compile_bitset_rules(network, components)
    state = 0
    for index, component in enumerate(components):
        if source_hypercube[component].value:
            state |= 1 << index

    visited = set()
    reached_target_configurations = set()
    n_fixed_target_components = sum(
        component in target_hypercube and target_hypercube[component].is_fixed
        for component in components
    )
    n_target_configurations = 1 << (len(components) - n_fixed_target_components)

    while state not in visited:
        if _configuration_bits_match_hypercube(state, target_hypercube, components):
            if quantifier != "universal":
                return True

            reached_target_configurations.add(state)
            if len(reached_target_configurations) == n_target_configurations:
                return True

        visited.add(state)
        state = _next_configuration_bits(compiled_rules, state)

    return False


def _synchronous_reachable_attractors(
    network: "BooleanNetwork",
    initial_configurations: ConfigurationSet,
) -> Tuple[ConfigurationSet, ...]:
    """Follow deterministic trajectories and return their terminal cycles."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    attractor_by_state: Dict[_ConfigurationBits, FrozenSet[_ConfigurationBits]] = {}
    attractors: Set[FrozenSet[_ConfigurationBits]] = set()

    for state in initial_configurations._iter_configuration_bits():
        if state in attractor_by_state:
            attractors.add(attractor_by_state[state])
            continue

        path: List[_ConfigurationBits] = []
        path_positions: Dict[_ConfigurationBits, int] = {}
        while state not in attractor_by_state and state not in path_positions:
            path_positions[state] = len(path)
            path.append(state)
            state = _next_configuration_bits(compiled_rules, state)

        if state in attractor_by_state:
            attractor = attractor_by_state[state]
        else:
            attractor = frozenset(path[path_positions[state] :])

        path_positions.clear()
        for path_state in path:
            attractor_by_state[path_state] = attractor

        attractors.add(attractor)

    return tuple(
        _configuration_set_from_configuration_bits(attractor, components)
        for attractor in attractors
    )


def _explicit_reachable_attractors(
    network: "BooleanNetwork",
    initial_configurations: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Enumerate terminal SCCs reachable in the explicit state-transition graph."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    if update == "asynchronous":
        terminal_components = _tarjan_asynchronous_terminal_sccs(
            (
                _configuration_bits_from_mapping(initial_configuration, components)
                for initial_configuration in initial_configurations
            ),
            lambda state: state ^ _next_configuration_bits(compiled_rules, state),
        )
        return tuple(
            _configuration_set_from_configuration_bits(component, components)
            for component in terminal_components
        )

    successors: Dict[_ConfigurationBits, Tuple[_ConfigurationBits, ...]] = {}
    for initial_configuration in initial_configurations:
        initial_configuration_bits = _configuration_bits_from_mapping(
            initial_configuration,
            components,
        )
        if initial_configuration_bits in successors:
            continue

        pending = [initial_configuration_bits]
        scheduled = {initial_configuration_bits}
        while pending:
            configuration_bits = pending.pop()
            scheduled.discard(configuration_bits)
            if configuration_bits in successors:
                continue

            state_successors = _successor_configuration_bits(
                compiled_rules,
                configuration_bits,
                update=update,
                n_components=len(components),
            )
            successors[configuration_bits] = state_successors

            for successor in state_successors:
                if successor not in successors and successor not in scheduled:
                    pending.append(successor)
                    scheduled.add(successor)

    terminal_components = _terminal_strongly_connected_components(successors)

    return tuple(
        _configuration_set_from_configuration_bits(component, components)
        for component in terminal_components
    )


def _synchronous_reachable_configurations(
    network: "BooleanNetwork",
    initial_configuration: Hypercube,
) -> Iterator[Configuration]:
    """Iterate over the unique synchronous trajectory until its first cycle."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    initial_configuration_bits = _configuration_bits_from_hypercube(
        initial_configuration,
        components,
    )

    def iterate() -> Iterator[Configuration]:
        state = initial_configuration_bits
        visited: Set[_ConfigurationBits] = set()
        while state not in visited:
            visited.add(state)
            yield {
                component: (state >> index) & 1
                for index, component in enumerate(components)
            }
            state = _next_configuration_bits(compiled_rules, state)

    return iterate()


def _explicit_reachable_configurations(
    network: "BooleanNetwork",
    initial_configuration: Hypercube,
    *,
    update: Literal["asynchronous", "general"],
) -> Iterator[Configuration]:
    """Traverse a finite reachable-state graph lazily with bitsets."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    initial_configuration_bits = _configuration_bits_from_hypercube(
        initial_configuration,
        components,
    )

    def iterate() -> Iterator[Configuration]:
        pending = [initial_configuration_bits]
        scheduled = {initial_configuration_bits}

        while pending:
            state = pending.pop()
            yield {
                component: (state >> index) & 1
                for index, component in enumerate(components)
            }

            next_state = _next_configuration_bits(compiled_rules, state)
            unstable_mask = state ^ next_state
            if update == "asynchronous":
                while unstable_mask:
                    bit = unstable_mask & -unstable_mask
                    successor = state ^ bit
                    unstable_mask ^= bit
                    if successor not in scheduled:
                        scheduled.add(successor)
                        pending.append(successor)
            else:
                updated_mask = unstable_mask
                while updated_mask:
                    successor = state ^ updated_mask
                    updated_mask = (updated_mask - 1) & unstable_mask
                    if successor not in scheduled:
                        scheduled.add(successor)
                        pending.append(successor)

    return iterate()


def _import_dd_backend() -> Any:
    """Import the fastest available BDD backend."""

    try:
        return import_module("dd.cudd")

    except (ImportError, OSError):
        pass

    return import_module("dd.autoref")


def _validate_hypercube(
    components: Tuple[str, ...],
    configuration: HypercubeLike,
    *,
    name: str,
) -> Hypercube:
    """Validate a possibly partial Boolean configuration."""

    if not isinstance(configuration, Mapping):
        raise TypeError(f"{name} must be a mapping")

    hypercube = Hypercube(configuration)
    unknown_components = sorted(set(hypercube) - set(components))
    if unknown_components:
        raise ValueError(f"unknown components in {name}: {unknown_components}")

    return hypercube


def _validate_complete_hypercube(
    components: Tuple[str, ...],
    configuration: HypercubeLike,
    *,
    name: str,
) -> Hypercube:
    """Validate a complete Boolean configuration encoded as a hypercube."""

    hypercube = _validate_hypercube(components, configuration, name=name)
    invalid_components = [
        component
        for component in components
        if component not in hypercube or not hypercube[component].is_fixed
    ]
    if invalid_components:
        raise ValueError(f"{name} must define fixed values for all network components")

    return hypercube


def _compile_bitset_rules(
    network: "BooleanNetwork",
    components: Tuple[str, ...],
) -> Tuple[_CompiledRule, ...]:
    """Compile Boolean rules into bitset evaluators."""

    component_indices = {component: index for index, component in enumerate(components)}

    return tuple(
        _compile_bitset_rule(network, rule, component_indices)
        for rule in network.values()
    )


def _compile_bitset_rule(
    network: "BooleanNetwork",
    rule: Any,
    component_indices: Mapping[str, int],
) -> _CompiledRule:
    """Compile one Boolean expression into a bitset evaluator."""

    source = _bitset_rule_source(network, rule, component_indices)
    # The generated source contains only integer indices and fixed operators.
    return eval(  # noqa: S307
        f"lambda state: {source}",
        {"__builtins__": {}},
    )


def _bitset_rule_source(
    network: "BooleanNetwork",
    rule: Any,
    component_indices: Mapping[str, int],
) -> str:
    """Translate one parsed Boolean expression into safe Python operations."""

    if rule is network.ba.TRUE:
        return "1"

    if rule is network.ba.FALSE:
        return "0"

    if isinstance(rule, network.ba.Symbol):
        index = component_indices[str(rule)]
        return f"((state >> {index}) & 1)"

    if isinstance(rule, network.ba.NOT):
        child = _bitset_rule_source(network, rule.args[0], component_indices)
        return f"(1 - {child})"

    if isinstance(rule, network.ba.AND):
        children = (
            _bitset_rule_source(network, child, component_indices)
            for child in rule.args
        )
        return f"({' and '.join(children)})"

    if isinstance(rule, network.ba.OR):
        children = (
            _bitset_rule_source(network, child, component_indices)
            for child in rule.args
        )
        return f"({' or '.join(children)})"

    raise TypeError(f"unsupported Boolean expression type: {type(rule)}")


def _configuration_bits_from_mapping(
    configuration: Mapping[str, int],
    components: Tuple[str, ...],
) -> _ConfigurationBits:
    """Encode a concrete configuration as an integer bitset."""

    configuration_bits = 0
    for index, component in enumerate(components):
        if configuration[component]:
            configuration_bits |= 1 << index

    return configuration_bits


def _configuration_bits_from_hypercube(
    hypercube: Hypercube,
    components: Tuple[str, ...],
) -> _ConfigurationBits:
    """Encode a fully specified hypercube as an integer bitset."""

    configuration_bits = 0
    for index, component in enumerate(components):
        if hypercube[component].value:
            configuration_bits |= 1 << index

    return configuration_bits


def _successor_configuration_bits(
    compiled_rules: Tuple[_CompiledRule, ...],
    configuration: _ConfigurationBits,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    n_components: int,
) -> Tuple[_ConfigurationBits, ...]:
    """Return successor bitsets under a finite-state update semantics."""

    next_configuration = _next_configuration_bits(compiled_rules, configuration)

    if update == "synchronous":
        return _synchronous_successor_configuration_bits(
            configuration,
            next_configuration,
        )

    unstable_mask = configuration ^ next_configuration
    if unstable_mask == 0:
        return ()

    if update == "asynchronous":
        return _asynchronous_successor_configuration_bits(
            configuration,
            unstable_mask,
        )

    return _general_successor_configuration_bits(configuration, unstable_mask)


def _next_configuration_bits(
    compiled_rules: Tuple[_CompiledRule, ...],
    configuration: _ConfigurationBits,
) -> _ConfigurationBits:
    """Evaluate all compiled rules on a bitset configuration."""

    next_configuration = 0
    for index, rule in enumerate(compiled_rules):
        if rule(configuration):
            next_configuration |= 1 << index

    return next_configuration


def _configuration_bits_match_hypercube(
    configuration: _ConfigurationBits,
    hypercube: Hypercube,
    components: Tuple[str, ...],
) -> bool:
    """Return whether a configuration bitset belongs to a hypercube."""

    for index, component in enumerate(components):
        if component not in hypercube:
            continue

        value = hypercube[component]
        if value.is_fixed and bool((configuration >> index) & 1) != bool(value.value):
            return False

    return True


def _configuration_set_from_configuration_bits(
    configuration_bits: Iterable[_ConfigurationBits],
    components: Tuple[str, ...],
) -> ConfigurationSet:
    """Convert explicit configuration bitsets into a ConfigurationSet."""

    hypercubes = _encoded_hypercubes_from_configuration_bits(
        frozenset(configuration_bits),
        n_components=len(components),
    )
    return ConfigurationSet._from_encoded_hypercubes(components, hypercubes)


def _encoded_hypercubes_from_configuration_bits(
    configurations: FrozenSet[_ConfigurationBits],
    *,
    n_components: int,
) -> FrozenSet[Tuple[_ConfigurationBits, _ConfigurationBits]]:
    """Compress explicit configurations into exact disjoint hypercubes."""

    all_mask = (1 << n_components) - 1
    hypercubes = {
        (all_mask, configuration & all_mask) for configuration in configurations
    }

    for index in reversed(range(n_components)):
        bit = 1 << index
        branches: Dict[Tuple[_ConfigurationBits, _ConfigurationBits], int] = {}
        for fixed_mask, value_mask in hypercubes:
            reduced = fixed_mask ^ bit, value_mask & ~bit
            branch = 2 if value_mask & bit else 1
            branches[reduced] = branches.get(reduced, 0) | branch

        hypercubes = set()
        for (fixed_mask, value_mask), branches_present in branches.items():
            if branches_present == 3:
                hypercubes.add((fixed_mask, value_mask))
                continue

            fixed_mask |= bit
            if branches_present == 2:
                value_mask |= bit
            hypercubes.add((fixed_mask, value_mask))

    return frozenset(hypercubes)


def _terminal_strongly_connected_components(
    successors: Mapping[_ConfigurationBits, Tuple[_ConfigurationBits, ...]],
) -> Tuple[Tuple[_ConfigurationBits, ...], ...]:
    """Return SCCs with no outgoing transition to another SCC."""

    components = _kosaraju_strongly_connected_components(successors)
    component_of = {
        state: component_index
        for component_index, component in enumerate(components)
        for state in component
    }
    terminal_components = []

    for component_index, component in enumerate(components):
        if all(
            component_of[successor] == component_index
            for state in component
            for successor in successors[state]
        ):
            terminal_components.append(tuple(component))

    return tuple(terminal_components)


def _kosaraju_strongly_connected_components(
    successors: Mapping[_ConfigurationBits, Tuple[_ConfigurationBits, ...]],
) -> Tuple[Tuple[_ConfigurationBits, ...], ...]:
    """Compute SCCs from a successor mapping with Kosaraju's algorithm."""

    visited: Set[_ConfigurationBits] = set()
    finishing_order: List[_ConfigurationBits] = []

    for start in successors:
        if start in visited:
            continue

        stack = [(start, False)]
        while stack:
            state, expanded = stack.pop()
            if expanded:
                finishing_order.append(state)
                continue

            if state in visited:
                continue

            visited.add(state)
            stack.append((state, True))
            for successor in successors[state]:
                if successor not in visited:
                    stack.append((successor, False))

    predecessors: Dict[_ConfigurationBits, List[_ConfigurationBits]] = {
        state: [] for state in successors
    }
    for state, state_successors in successors.items():
        for successor in state_successors:
            predecessors[successor].append(state)

    assigned: Set[_ConfigurationBits] = set()
    components: List[Tuple[_ConfigurationBits, ...]] = []
    for start in reversed(finishing_order):
        if start in assigned:
            continue

        component = []
        stack = [start]
        assigned.add(start)

        while stack:
            state = stack.pop()
            component.append(state)
            for predecessor in predecessors[state]:
                if predecessor not in assigned:
                    assigned.add(predecessor)
                    stack.append(predecessor)

        components.append(tuple(component))

    return tuple(components)
