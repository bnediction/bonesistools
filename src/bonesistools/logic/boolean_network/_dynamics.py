#!/usr/bin/env python

from __future__ import annotations

from importlib import import_module
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
)

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import HypercubeLike
from ._asynchronous import (
    _asynchronous_successor_state_bits,
    _bdd_asynchronous_transition_relation,
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


class _BDDTransitionSystem:
    """Compiled symbolic transition system backed by binary decision diagrams."""

    def __init__(
        self,
        network: "BooleanNetwork",
        *,
        update: Literal["asynchronous", "synchronous", "general"],
    ) -> None:

        bdd_module = _import_dd_backend()
        self._bdd = bdd_module.BDD()
        self._and_exists = getattr(bdd_module, "and_exists", None)
        self._update = update
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
        self._transition = None
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
            self._transition = _bdd_transition_relation(
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

    def reachable_state_bits(
        self,
        initial_states: ConfigurationSet,
    ) -> Tuple[_StateBits, ...]:
        """Return bitsets for states reachable from initial configurations."""

        initial = self._encode_configurations(initial_states)
        reachable = self._forward_reachable_states(initial)
        return self._decode_state_bits(reachable)

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
        for hypercube in configurations._as_hypercubes():
            states |= self._encode_hypercube(hypercube)

        return states

    def _forward_reachable_states(self, initial: Any) -> Any:
        """Compute the symbolic forward reachability fixed point."""

        reachable = initial
        frontier = initial

        while frontier != self._bdd.false:
            successors = self._successors(frontier)
            new_frontier = successors & ~reachable
            if new_frontier == self._bdd.false:
                break

            reachable |= new_frontier
            frontier = new_frontier

        return reachable

    def _successors(self, configurations: Any) -> Any:
        """Return symbolic successors of a configuration set."""

        if self._update in {"synchronous", "general"}:
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
        else:
            successors = self._relational_product(
                configurations,
                self._transition,
                self._current_variables,
            )

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
        if self._update in {"synchronous", "general"}:
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

        return self._relational_product(
            self._transition,
            configurations_next,
            self._next_variables,
        )

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

    def _decode_state_bits(self, states: Any) -> Tuple[_StateBits, ...]:
        """Decode a BDD state set into sorted bitsets."""

        decoded = []
        for assignment in self._bdd.pick_iter(
            states,
            care_vars=self._current_variables,
        ):
            state_bits = 0
            for index, current_variable in enumerate(self._current_variables):
                if assignment.get(current_variable, False):
                    state_bits |= 1 << index
            decoded.append(state_bits)

        return tuple(sorted(decoded))


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
    components = transition_system.components
    reachable_state_bits = transition_system.reachable_state_bits(initial_states)

    compiled_rules = _compile_bitset_rules(network, components)
    successors = {
        state: _successor_state_bits(
            compiled_rules,
            state,
            update=update,
            n_components=len(components),
        )
        for state in reachable_state_bits
    }
    terminal_components = _terminal_strongly_connected_components(successors)

    return tuple(
        _configuration_set_from_state_bits(component, components)
        for component in sorted(terminal_components, key=lambda states: min(states))
    )


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


def _explicit_reachable_attractors(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Enumerate terminal SCCs reachable in the explicit state-transition graph."""

    components = tuple(network.keys())
    compiled_rules = _compile_bitset_rules(network, components)
    successors: Dict[_StateBits, Tuple[_StateBits, ...]] = {}
    pending = [
        _state_bits_from_mapping(initial_state, components)
        for initial_state in initial_states
    ]
    scheduled = set(pending)

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
        for component in sorted(terminal_components, key=lambda states: min(states))
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


def _bdd_transition_relation(
    network: "BooleanNetwork",
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Any:
    """Build the asynchronous transition relation as a BDD."""

    return _bdd_asynchronous_transition_relation(
        network,
        bdd=bdd,
        components=components,
        current_vars=current_vars,
        next_vars=next_vars,
    )


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

    if rule is network.ba.TRUE:
        return lambda state: 1

    if rule is network.ba.FALSE:
        return lambda state: 0

    if isinstance(rule, network.ba.Symbol):
        index = component_indices[str(rule)]
        return lambda state, index=index: (state >> index) & 1

    if isinstance(rule, network.ba.NOT):
        child = _compile_bitset_rule(network, rule.args[0], component_indices)
        return lambda state, child=child: 1 - child(state)

    if isinstance(rule, network.ba.AND):
        children = tuple(
            _compile_bitset_rule(network, child, component_indices)
            for child in rule.args
        )
        return lambda state, children=children: int(
            all(child(state) for child in children)
        )

    if isinstance(rule, network.ba.OR):
        children = tuple(
            _compile_bitset_rule(network, child, component_indices)
            for child in rule.args
        )
        return lambda state, children=children: int(
            any(child(state) for child in children)
        )

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

    states = tuple(sorted(state_bits))
    hypercube = _complete_hypercube_from_state_bits(
        states,
        n_components=len(components),
    )
    if hypercube is not None:
        return ConfigurationSet(
            components,
            [_partial_state_mapping_from_bits(hypercube, components)],
        )

    configurations = ConfigurationSet(components)
    for state in states:
        configurations.add(_state_mapping_from_bits(state, components))

    configurations.compress()

    return configurations


def _complete_hypercube_from_state_bits(
    states: Tuple[_StateBits, ...],
    *,
    n_components: int,
) -> Optional[Tuple[_StateBits, _StateBits]]:
    """Return fixed bits when states cover one complete hypercube."""

    if not states:
        return None

    all_mask = (1 << n_components) - 1
    common_ones = all_mask
    common_zeros = all_mask

    for state in states:
        common_ones &= state
        common_zeros &= ~state & all_mask

    fixed_mask = common_ones | common_zeros
    expected_size = 1 << (n_components - fixed_mask.bit_count())
    if len(states) != expected_size:
        return None

    fixed_values = common_ones & fixed_mask
    if all((state & fixed_mask) == fixed_values for state in states):
        return fixed_mask, fixed_values

    return None


def _partial_state_mapping_from_bits(
    hypercube: Tuple[_StateBits, _StateBits],
    components: Tuple[str, ...],
) -> Dict[str, int]:
    """Decode fixed bits into a partial configuration mapping."""

    fixed_mask, value_mask = hypercube
    return {
        component: 1 if value_mask & (1 << index) else 0
        for index, component in enumerate(components)
        if fixed_mask & (1 << index)
    }


def _state_mapping_from_bits(
    state: _StateBits,
    components: Tuple[str, ...],
) -> Dict[str, int]:
    """Decode a bitset state into a component-value mapping."""

    return {
        component: (state >> index) & 1 for index, component in enumerate(components)
    }


def _terminal_strongly_connected_components(
    successors: Mapping[_StateBits, Tuple[_StateBits, ...]],
) -> Tuple[Tuple[_StateBits, ...], ...]:
    """Return SCCs with no outgoing transition to another SCC."""

    components = _strongly_connected_components(successors)
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
            terminal_components.append(tuple(sorted(component)))

    return tuple(terminal_components)


def _strongly_connected_components(
    successors: Mapping[_StateBits, Tuple[_StateBits, ...]],
) -> Tuple[Tuple[_StateBits, ...], ...]:
    """Compute SCCs of a directed graph encoded as a successor mapping."""

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
