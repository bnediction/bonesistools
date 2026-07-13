#!/usr/bin/env python

from __future__ import annotations

from importlib import import_module
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Mapping,
    Set,
    Tuple,
)

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import HypercubeLike
from ._asynchronous import (
    _asynchronous_successor_state_bits,
    _bdd_asynchronous_transition_partitions,
    _tarjan_asynchronous_terminal_sccs,
)
from ._bdd import _bdd_equivalence
from ._general import (
    _bdd_general_transition_partitions,
    _general_successor_state_bits,
)
from ._synchronous import (
    _bdd_synchronous_transition_partitions,
    _synchronous_successor_state_bits,
)

if TYPE_CHECKING:
    from ._network import BooleanNetwork

_StateBits = int
_CompiledRule = Callable[[_StateBits], int]
_SYMBOLIC_SCC_THRESHOLD = 128


class _BDDTransitionSystem:
    """Compiled symbolic transition system backed by binary decision diagrams."""

    def __init__(
        self,
        network: "BooleanNetwork",
        *,
        update: Literal["asynchronous", "synchronous", "general"],
    ) -> None:

        bdd_module = _import_dd_backend()
        self._network = network
        self._bdd = bdd_module.BDD()
        self._and_exists = getattr(bdd_module, "and_exists", None)
        self._update: Literal["asynchronous", "synchronous", "general"] = update
        self._components = tuple(network.keys())
        self._current_variables = tuple(
            f"x{index}" for index in range(len(self._components))
        )
        self._next_variables = tuple(
            f"y{index}" for index in range(len(self._components))
        )
        ordered_variables = tuple(
            variable
            for pair in zip(self._current_variables, self._next_variables)
            for variable in pair
        )
        self._bdd.declare(*ordered_variables)
        self._disjunctive_transitions: Tuple[Any, ...] = ()
        self._transition_clusters: Tuple[Any, ...] = ()
        self._forward_quantification: Tuple[Tuple[str, ...], ...] = ()
        if update in {"synchronous", "general"}:
            transition_partitions = (
                _bdd_synchronous_transition_partitions
                if update == "synchronous"
                else _bdd_general_transition_partitions
            )
            (
                self._transition_clusters,
                self._forward_quantification,
            ) = transition_partitions(
                network,
                bdd=self._bdd,
                components=self._components,
                current_vars=self._current_variables,
                next_vars=self._next_variables,
            )
        else:
            self._disjunctive_transitions = _bdd_asynchronous_transition_partitions(
                network,
                bdd=self._bdd,
                components=self._components,
                current_vars=self._current_variables,
                next_vars=self._next_variables,
            )
        self._rename_current_to_next = dict(
            zip(self._current_variables, self._next_variables)
        )
        self._rename_next_to_current = dict(
            zip(self._next_variables, self._current_variables)
        )

    @property
    def components(self) -> Tuple[str, ...]:
        """Return components in BDD variable order."""

        return self._components

    def reachability(
        self,
        initial_state: HypercubeLike,
        target_state: HypercubeLike,
        *,
        quantifier: Literal["exists", "robust", "universal"],
    ) -> bool:
        """Test whether initial configurations can reach a target subspace."""

        initial_hypercube = _validate_hypercube(
            self._components,
            initial_state,
            name="initial_state",
        )
        target_hypercube = _validate_hypercube(
            self._components,
            target_state,
            name="target_state",
        )
        initial = self._encode_hypercube(initial_hypercube)
        target = self._encode_hypercube(target_hypercube)
        reachable_from_initial = self._forward_reachable_states(initial)
        reachable_target = reachable_from_initial & target

        if reachable_target == self._bdd.false:
            return False

        if quantifier == "universal":
            return self._all_initial_states_reach_all_target_states(
                initial,
                target_hypercube,
                reachable_from_initial,
            )

        initial_is_concrete = (
            initial_hypercube.components == frozenset(self._components)
            and initial_hypercube.is_fully_specified
        )
        if quantifier == "exists" or initial_is_concrete:
            return True

        states_reaching_target = reachable_target

        while True:
            if self._matches_reachability_quantifier(
                initial,
                states_reaching_target,
                quantifier=quantifier,
            ):
                return True

            predecessors = self._predecessors(states_reaching_target)
            predecessors &= reachable_from_initial
            updated = states_reaching_target | predecessors
            if updated == states_reaching_target:
                return False

            states_reaching_target = updated

    def reachable_attractors(
        self,
        initial_states: ConfigurationSet,
    ) -> Tuple[ConfigurationSet, ...]:
        """Return reachable terminal SCCs without decoding transient states."""

        initial = self._encode_configurations(initial_states)
        reachable = self._forward_reachable_states(initial)
        n_reachable = self._bdd.count(
            reachable,
            nvars=len(self._components),
        )
        if n_reachable <= _SYMBOLIC_SCC_THRESHOLD:
            return self._explicit_terminal_sccs(reachable)

        if self._update == "synchronous":
            return tuple(
                self._decode_configuration_set(cycle)
                for cycle in self._synchronous_terminal_cycles(reachable)
            )

        return tuple(
            self._decode_configuration_set(component)
            for component in self._terminal_sccs(reachable)
        )

    def _encode_hypercube(self, hypercube: Hypercube) -> Any:
        """Encode one validated Boolean hypercube."""

        return self._encode_hypercube_variables(
            hypercube,
            self._current_variables,
        )

    def _encode_hypercube_variables(
        self,
        hypercube: Hypercube,
        variables: Tuple[str, ...],
    ) -> Any:
        """Encode one hypercube over a selected BDD variable family."""

        encoded = self._bdd.true
        for component, variable_name in zip(self._components, variables):
            if component not in hypercube:
                continue

            value = hypercube[component]
            if value.is_fixed:
                variable = self._bdd.var(variable_name)
                encoded &= variable if value.value else ~variable

        return encoded

    def _encode_configurations(self, configurations: ConfigurationSet) -> Any:
        """Encode a ConfigurationSet as a BDD over current variables."""

        states = self._bdd.false
        for fixed_mask, value_mask in configurations._iter_encoded_hypercubes():
            encoded = self._bdd.true
            for index, variable_name in enumerate(self._current_variables):
                bit = 1 << index
                if not fixed_mask & bit:
                    continue

                variable = self._bdd.var(variable_name)
                encoded &= variable if value_mask & bit else ~variable
            states |= encoded

        return states

    def _forward_reachable_states(self, initial: Any, within: Any = None) -> Any:
        """Compute the symbolic forward reachability fixed point."""

        reachable = initial
        frontier = initial

        while frontier != self._bdd.false:
            successors = self._successors(frontier)
            if within is not None:
                successors &= within
            new_frontier = successors & ~reachable
            if new_frontier == self._bdd.false:
                break

            reachable |= new_frontier
            frontier = new_frontier

        return reachable

    def _backward_reachable_states(self, initial: Any, within: Any) -> Any:
        """Compute a symbolic backward reachability fixed point within a region."""

        reachable = initial
        frontier = initial

        while frontier != self._bdd.false:
            predecessors = self._predecessors(frontier) & within
            new_frontier = predecessors & ~reachable
            if new_frontier == self._bdd.false:
                break

            reachable |= new_frontier
            frontier = new_frontier

        return reachable

    def _synchronous_terminal_cycles(self, states: Any) -> Tuple[Any, ...]:
        """Extract cycles from a deterministic synchronous transition system."""

        remaining = self._synchronous_recurrent_states(states)
        cycles = []

        while remaining != self._bdd.false:
            seed = self._pick_state(remaining)
            cycle = seed
            successor = self._successors(seed)

            while successor & cycle == self._bdd.false:
                cycle |= successor
                successor = self._successors(successor)

            cycles.append(cycle)
            remaining &= ~cycle

        return tuple(cycles)

    def _synchronous_recurrent_states(self, states: Any) -> Any:
        """Remove transient states until only synchronous cycles remain."""

        recurrent = states

        while True:
            updated = recurrent & self._successors(recurrent)
            if updated == recurrent:
                return recurrent
            recurrent = updated

    def _terminal_sccs(self, states: Any) -> Tuple[Any, ...]:
        """Extract terminal SCCs while keeping state regions symbolic."""

        remaining = states
        terminal_components = []
        while remaining != self._bdd.false:
            component = self._terminal_scc(remaining)
            terminal_components.append(component)
            basin = self._backward_reachable_states(component, remaining)
            remaining &= ~basin

        return tuple(terminal_components)

    def _terminal_scc(self, states: Any) -> Any:
        """Follow the symbolic SCC graph until reaching one terminal SCC."""

        region = states
        seed = self._pick_state(region)

        while True:
            forward = self._forward_reachable_states(seed, within=region)
            backward = self._backward_reachable_states(seed, forward)
            component = forward & backward
            successors = self._successors(component) & region & ~component
            if successors == self._bdd.false:
                return component

            region = forward & ~backward
            seed = self._pick_state(successors)

    def _explicit_terminal_sccs(self, states: Any) -> Tuple[ConfigurationSet, ...]:
        """Decode a small reachable set and compute its terminal SCCs explicitly."""

        state_bits = self._decode_state_bits(states)
        compiled_rules = _compile_bitset_rules(self._network, self._components)
        successors = {
            state: _successor_state_bits(
                compiled_rules,
                state,
                update=self._update,
                n_components=len(self._components),
            )
            for state in state_bits
        }
        return tuple(
            _configuration_set_from_state_bits(component, self._components)
            for component in _terminal_strongly_connected_components(successors)
        )

    def _successors(self, configurations: Any) -> Any:
        """Return symbolic successors of a configuration set."""

        if self._update == "asynchronous":
            successors = self._bdd.false
            for transition in self._disjunctive_transitions:
                successors |= self._relational_product(
                    configurations,
                    transition,
                    self._current_variables,
                )
        else:
            successors = configurations
            for cluster, quantified_variables in zip(
                self._transition_clusters,
                self._forward_quantification,
            ):
                if quantified_variables:
                    successors = self._relational_product(
                        successors,
                        cluster,
                        quantified_variables,
                    )
                else:
                    successors &= cluster
        return self._bdd.let(
            self._rename_next_to_current,
            successors,
        )

    def _predecessors(self, configurations: Any) -> Any:
        """Return symbolic predecessors of a configuration set."""

        configurations_next = self._bdd.let(
            self._rename_current_to_next,
            configurations,
        )
        if self._update == "asynchronous":
            predecessors = self._bdd.false
            for transition in self._disjunctive_transitions:
                predecessors |= self._relational_product(
                    transition,
                    configurations_next,
                    self._next_variables,
                )
            return predecessors

        predecessors = configurations_next
        for cluster, next_variable in zip(
            self._transition_clusters,
            self._next_variables,
        ):
            predecessors = self._relational_product(
                predecessors,
                cluster,
                (next_variable,),
            )
        return predecessors

    def _relational_product(
        self,
        left: Any,
        right: Any,
        quantified_variables: Tuple[str, ...],
    ) -> Any:
        """Conjoin two BDDs while existentially quantifying variables."""

        if self._and_exists is not None:
            return self._and_exists(left, right, quantified_variables)

        return self._bdd.exist(quantified_variables, left & right)

    def _all_initial_states_reach_all_target_states(
        self,
        initial: Any,
        target_hypercube: Hypercube,
        reachable_from_initial: Any,
    ) -> bool:
        """Test universal reachability without enumerating target states."""

        target_variables = tuple(f"z{index}" for index in range(len(self._components)))
        undeclared_variables = tuple(
            variable for variable in target_variables if variable not in self._bdd.vars
        )
        if undeclared_variables:
            self._bdd.declare(*undeclared_variables)

        target_states = self._encode_hypercube_variables(
            target_hypercube,
            target_variables,
        )
        current_matches_target = target_states
        for current_variable, target_variable in zip(
            self._current_variables,
            target_variables,
        ):
            current_matches_target &= _bdd_equivalence(
                self._bdd.var(current_variable),
                self._bdd.var(target_variable),
            )

        required_pairs = initial & target_states
        states_reaching_targets = current_matches_target & reachable_from_initial

        while True:
            if required_pairs & ~states_reaching_targets == self._bdd.false:
                return True

            predecessors = self._predecessors(states_reaching_targets)
            predecessors &= reachable_from_initial
            predecessors &= target_states
            updated = states_reaching_targets | predecessors
            if updated == states_reaching_targets:
                return False

            states_reaching_targets = updated

    def _matches_reachability_quantifier(
        self,
        initial: Any,
        states_reaching_target: Any,
        *,
        quantifier: Literal["exists", "robust"],
    ) -> bool:
        """Test whether the current backward closure answers the query."""

        if quantifier == "exists":
            return initial & states_reaching_target != self._bdd.false

        return initial & ~states_reaching_target == self._bdd.false

    def _pick_state(self, states: Any) -> Any:
        """Select one concrete state from a non-empty symbolic state set."""

        assignment = next(
            self._bdd.pick_iter(
                states,
                care_vars=self._current_variables,
            )
        )
        selected = self._bdd.true
        for variable_name in self._current_variables:
            variable = self._bdd.var(variable_name)
            selected &= variable if assignment[variable_name] else ~variable

        return selected

    def _decode_configuration_set(self, states: Any) -> ConfigurationSet:
        """Decode a symbolic state set into exact disjoint hypercubes."""

        hypercubes = []
        for assignment in self._bdd.pick_iter(states):
            fixed_mask = 0
            value_mask = 0
            for index, variable_name in enumerate(self._current_variables):
                if variable_name not in assignment:
                    continue

                bit = 1 << index
                fixed_mask |= bit
                if assignment[variable_name]:
                    value_mask |= bit
            hypercubes.append((fixed_mask, value_mask))

        return ConfigurationSet._from_encoded_hypercubes(
            self._components,
            hypercubes,
        )

    def _decode_state_bits(self, states: Any) -> Tuple[_StateBits, ...]:
        """Decode a symbolic state set into explicit bitsets."""

        decoded = []
        for assignment in self._bdd.pick_iter(
            states,
            care_vars=self._current_variables,
        ):
            state_bits = 0
            for index, variable_name in enumerate(self._current_variables):
                if assignment[variable_name]:
                    state_bits |= 1 << index
            decoded.append(state_bits)

        return tuple(decoded)


def _bdd_reachable_attractors(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Compute reachable attractors through a BDD reachability set."""

    transition_system = _BDDTransitionSystem(
        network,
        update=update,
    )
    return transition_system.reachable_attractors(initial_states)


def _bdd_reachability(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test reachability through a symbolic BDD backward closure."""

    transition_system = _BDDTransitionSystem(
        network,
        update=update,
    )
    return transition_system.reachability(
        initial_state,
        target_state,
        quantifier=quantifier,
    )


def _synchronous_reachability(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
    *,
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Use a direct trajectory when the synchronous initial state is concrete."""

    components = tuple(network.keys())
    initial_hypercube = _validate_hypercube(
        components,
        initial_state,
        name="initial_state",
    )
    if (
        initial_hypercube.components != frozenset(components)
        or not initial_hypercube.is_fully_specified
    ):
        return _bdd_reachability(
            network,
            initial_hypercube,
            target_state,
            update="synchronous",
            quantifier=quantifier,
        )

    target_hypercube = _validate_hypercube(
        components,
        target_state,
        name="target_state",
    )
    compiled_rules = _compile_bitset_rules(network, components)
    state = 0
    for index, component in enumerate(components):
        if initial_hypercube[component].value:
            state |= 1 << index

    visited = set()
    reached_target_states = set()
    n_fixed_target_components = sum(
        component in target_hypercube and target_hypercube[component].is_fixed
        for component in components
    )
    n_target_states = 1 << (len(components) - n_fixed_target_components)

    while state not in visited:
        if _state_bits_match_hypercube(state, target_hypercube, components):
            if quantifier != "universal":
                return True

            reached_target_states.add(state)
            if len(reached_target_states) == n_target_states:
                return True

        visited.add(state)
        state = _next_state_bits(compiled_rules, state)

    return False


def _synchronous_reachable_attractors(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
) -> Tuple[ConfigurationSet, ...]:
    """Follow deterministic trajectories and return their terminal cycles."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    attractor_by_state: Dict[_StateBits, FrozenSet[_StateBits]] = {}
    attractors: Set[FrozenSet[_StateBits]] = set()

    for state in initial_states._iter_state_bits():
        if state in attractor_by_state:
            attractors.add(attractor_by_state[state])
            continue

        path: List[_StateBits] = []
        path_positions: Dict[_StateBits, int] = {}
        while state not in attractor_by_state and state not in path_positions:
            path_positions[state] = len(path)
            path.append(state)
            state = _next_state_bits(compiled_rules, state)

        if state in attractor_by_state:
            attractor = attractor_by_state[state]
        else:
            attractor = frozenset(path[path_positions[state] :])

        path_positions.clear()
        for path_state in path:
            attractor_by_state[path_state] = attractor

        attractors.add(attractor)

    return tuple(
        _configuration_set_from_state_bits(attractor, components)
        for attractor in attractors
    )


def _explicit_reachable_attractors(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Enumerate terminal SCCs reachable in the explicit state-transition graph."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    if update == "asynchronous":
        terminal_components = _tarjan_asynchronous_terminal_sccs(
            (
                _state_bits_from_mapping(initial_state, components)
                for initial_state in initial_states
            ),
            lambda state: state ^ _next_state_bits(compiled_rules, state),
        )
        return tuple(
            _configuration_set_from_state_bits(component, components)
            for component in terminal_components
        )

    successors: Dict[_StateBits, Tuple[_StateBits, ...]] = {}
    for initial_state in initial_states:
        initial_state_bits = _state_bits_from_mapping(initial_state, components)
        if initial_state_bits in successors:
            continue

        pending = [initial_state_bits]
        scheduled = {initial_state_bits}
        while pending:
            state_bits = pending.pop()
            scheduled.discard(state_bits)
            if state_bits in successors:
                continue

            state_successors = _successor_state_bits(
                compiled_rules,
                state_bits,
                update=update,
                n_components=len(components),
            )
            successors[state_bits] = state_successors

            for successor in state_successors:
                if successor not in successors and successor not in scheduled:
                    pending.append(successor)
                    scheduled.add(successor)

    terminal_components = _terminal_strongly_connected_components(successors)

    return tuple(
        _configuration_set_from_state_bits(component, components)
        for component in terminal_components
    )


def _import_dd_backend() -> Any:
    """Import the fastest available optional BDD backend."""

    try:
        return import_module("dd.cudd")

    except (ImportError, OSError):
        pass

    try:
        return import_module("dd.autoref")

    except ImportError as error:
        raise ImportError(
            "BDD-based Boolean dynamics require the optional dependency `dd`; "
            "install it with `pip install 'bonesistools[bdd]'`."
        ) from error


def _validate_hypercube(
    components: Tuple[str, ...],
    state: HypercubeLike,
    *,
    name: str,
) -> Hypercube:
    """Validate a possibly partial Boolean state."""

    hypercube = Hypercube(state)
    unknown_components = sorted(set(hypercube) - set(components))
    if unknown_components:
        raise ValueError(f"unknown components in {name}: {unknown_components}")

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


def _state_bits_from_mapping(
    state: Mapping[str, int],
    components: Tuple[str, ...],
) -> _StateBits:
    """Encode a concrete state mapping as an integer bitset."""

    state_bits = 0
    for index, component in enumerate(components):
        if state[component]:
            state_bits |= 1 << index

    return state_bits


def _successor_state_bits(
    compiled_rules: Tuple[_CompiledRule, ...],
    state: _StateBits,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    n_components: int,
) -> Tuple[_StateBits, ...]:
    """Return successor bitsets under a finite-state update semantics."""

    next_state = _next_state_bits(compiled_rules, state)

    if update == "synchronous":
        return _synchronous_successor_state_bits(state, next_state)

    unstable_mask = state ^ next_state
    if unstable_mask == 0:
        return ()

    if update == "asynchronous":
        return _asynchronous_successor_state_bits(state, unstable_mask)

    return _general_successor_state_bits(state, unstable_mask)


def _next_state_bits(
    compiled_rules: Tuple[_CompiledRule, ...],
    state: _StateBits,
) -> _StateBits:
    """Evaluate all compiled rules on a bitset state."""

    next_state = 0
    for index, rule in enumerate(compiled_rules):
        if rule(state):
            next_state |= 1 << index

    return next_state


def _state_bits_match_hypercube(
    state: _StateBits,
    hypercube: Hypercube,
    components: Tuple[str, ...],
) -> bool:
    """Return whether a state bitset belongs to a hypercube."""

    for index, component in enumerate(components):
        if component not in hypercube:
            continue

        value = hypercube[component]
        if value.is_fixed and bool((state >> index) & 1) != bool(value.value):
            return False

    return True


def _configuration_set_from_state_bits(
    state_bits: Iterable[_StateBits],
    components: Tuple[str, ...],
) -> ConfigurationSet:
    """Convert explicit state bitsets into a compact ConfigurationSet."""

    hypercubes = _encoded_hypercubes_from_state_bits(
        frozenset(state_bits),
        n_components=len(components),
    )
    return ConfigurationSet._from_encoded_hypercubes(components, hypercubes)


def _encoded_hypercubes_from_state_bits(
    states: FrozenSet[_StateBits],
    *,
    n_components: int,
) -> FrozenSet[Tuple[_StateBits, _StateBits]]:
    """Compress explicit states into exact disjoint encoded hypercubes."""

    all_mask = (1 << n_components) - 1
    hypercubes = {(all_mask, state & all_mask) for state in states}

    for index in reversed(range(n_components)):
        bit = 1 << index
        branches: Dict[Tuple[_StateBits, _StateBits], int] = {}
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
    successors: Mapping[_StateBits, Tuple[_StateBits, ...]],
) -> Tuple[Tuple[_StateBits, ...], ...]:
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
    successors: Mapping[_StateBits, Tuple[_StateBits, ...]],
) -> Tuple[Tuple[_StateBits, ...], ...]:
    """Compute SCCs from a successor mapping with Kosaraju's algorithm."""

    visited: Set[_StateBits] = set()
    finishing_order: List[_StateBits] = []

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

    predecessors: Dict[_StateBits, List[_StateBits]] = {
        state: [] for state in successors
    }
    for state, state_successors in successors.items():
        for successor in state_successors:
            predecessors[successor].append(state)

    assigned: Set[_StateBits] = set()
    components: List[Tuple[_StateBits, ...]] = []
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
