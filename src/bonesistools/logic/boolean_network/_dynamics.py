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
    Optional,
    Set,
    Tuple,
    cast,
)

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import HypercubeLike
from ._most_permissive import (
    _asp_string,
    _boolean_function_asp_facts,
    _boolean_function_evaluation_asp_rules,
    _validate_hypercube,
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
        self._general_updates = None
        self._general_clusters: Tuple[Any, ...] = ()
        self._general_forward_quantification: Tuple[Tuple[str, ...], ...] = ()
        if update == "general":
            (
                self._general_updates,
                self._general_clusters,
                self._general_forward_quantification,
            ) = _bdd_general_transition_partitions(
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
                update=update,
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
        quantifier: Literal["exists", "robust"],
    ) -> bool:
        """Test whether initial configurations can reach a target subspace."""

        initial_hypercube = _validate_hypercube(
            self._components,
            initial_state,
            name="initial_state",
        )
        initial = self._encode_hypercube(initial_hypercube)
        target = self._encode_configuration(target_state, name="target_state")
        reachable_from_initial = self._forward_reachable_states(initial)
        reachable_target = reachable_from_initial & target

        if reachable_target == self._bdd.false:
            return False

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

    def _encode_configuration(
        self,
        configuration: HypercubeLike,
        *,
        name: str,
    ) -> Any:
        """Encode one possibly partial Boolean configuration."""

        hypercube = _validate_hypercube(
            self._components,
            configuration,
            name=name,
        )
        return self._encode_hypercube(hypercube)

    def _encode_hypercube(self, hypercube: Hypercube) -> Any:
        """Encode one validated Boolean hypercube."""

        return self._encode_configurations(
            ConfigurationSet(self._components, [hypercube])
        )

    def _encode_configurations(self, configurations: ConfigurationSet) -> Any:
        """Encode a ConfigurationSet as a BDD over current variables."""

        states = self._bdd.false
        for hypercube in configurations._as_hypercubes():
            state = self._bdd.true
            for component, current_variable in zip(
                self._components,
                self._current_variables,
            ):
                if component not in hypercube:
                    continue

                value = hypercube[component]
                if value.is_fixed:
                    variable = self._bdd.var(current_variable)
                    state &= variable if value.value else ~variable

            states |= state

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

        if self._update == "general":
            successors = configurations & self._general_updates
            for cluster, quantified_variables in zip(
                self._general_clusters,
                self._general_forward_quantification,
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
        if self._update == "general":
            predecessors = configurations_next & self._general_updates
            for cluster, next_variable in zip(
                self._general_clusters,
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


def _reachable_attractors_with_bdd_backend(
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


def _reachability_with_bdd_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
    *,
    update: Literal["asynchronous", "general"],
    quantifier: Literal["exists", "robust"],
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


def _reachable_attractors_with_explicit_backend(
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
            "`BooleanNetwork.reachable_attractors(..., backend='bdd')` requires "
            "the optional dependency `dd` to be installed."
        ) from error


def _bdd_transition_relation(
    network: "BooleanNetwork",
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
    update: Literal["asynchronous", "synchronous"],
) -> Any:
    """Build a transition relation as a BDD."""

    component_vars = dict(zip(components, current_vars))
    current_values = tuple(bdd.var(variable) for variable in current_vars)
    next_values = tuple(bdd.var(variable) for variable in next_vars)
    rules = tuple(
        _bdd_rule(network, network[component], bdd=bdd, variables=component_vars)
        for component in components
    )

    if update == "synchronous":
        relation = bdd.true
        for next_value, rule in zip(next_values, rules):
            relation &= _bdd_equivalence(next_value, rule)

        return relation

    relation = bdd.false
    unchanged_prefixes = [bdd.true]
    for current_value, next_value in zip(current_values, next_values):
        unchanged_prefixes.append(
            unchanged_prefixes[-1] & _bdd_equivalence(next_value, current_value)
        )

    unchanged_suffixes = [bdd.true] * (len(components) + 1)
    for index in range(len(components) - 1, -1, -1):
        unchanged_suffixes[index] = unchanged_suffixes[index + 1] & (
            _bdd_equivalence(next_values[index], current_values[index])
        )

    for updated_index, (current_value, next_value, rule) in enumerate(
        zip(current_values, next_values, rules)
    ):
        transition = _bdd_exclusive_or(current_value, rule)
        transition &= _bdd_equivalence(next_value, rule)
        transition &= unchanged_prefixes[updated_index]
        transition &= unchanged_suffixes[updated_index + 1]

        relation |= transition

    return relation


def _bdd_general_transition_partitions(
    network: "BooleanNetwork",
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Tuple[Any, Tuple[Any, ...], Tuple[Tuple[str, ...], ...]]:
    """Build conjunctive partitions for the general transition relation."""

    component_vars = dict(zip(components, current_vars))
    rules = tuple(
        _bdd_rule(network, network[component], bdd=bdd, variables=component_vars)
        for component in components
    )
    clusters = []
    updates = bdd.false
    for current_var, next_var, rule in zip(current_vars, next_vars, rules):
        current_value = bdd.var(current_var)
        next_value = bdd.var(next_var)
        updated = _bdd_exclusive_or(current_value, rule)
        updated &= _bdd_equivalence(next_value, rule)
        unchanged = _bdd_equivalence(next_value, current_value)
        clusters.append(unchanged | updated)
        updates |= updated

    current_variables = set(current_vars)
    last_cluster = {variable: -1 for variable in current_vars}
    for index, cluster in enumerate(clusters):
        for variable in bdd.support(cluster) & current_variables:
            last_cluster[variable] = index

    forward_quantification = tuple(
        tuple(variable for variable in current_vars if last_cluster[variable] == index)
        for index in range(len(clusters))
    )

    return updates, tuple(clusters), forward_quantification


def _bdd_equivalence(left: Any, right: Any) -> Any:
    """Return a BDD encoding logical equivalence."""

    return (left & right) | (~left & ~right)


def _bdd_exclusive_or(left: Any, right: Any) -> Any:
    """Return a BDD encoding exclusive disjunction."""

    return (left & ~right) | (~left & right)


def _bdd_rule(
    network: "BooleanNetwork",
    rule: Any,
    *,
    bdd: Any,
    variables: Mapping[str, str],
) -> Any:
    """Convert a Boolean rule into a BDD over current-state variables."""

    if rule is network.ba.TRUE:
        return bdd.true

    if rule is network.ba.FALSE:
        return bdd.false

    if isinstance(rule, network.ba.Symbol):
        return bdd.var(variables[str(rule)])

    if isinstance(rule, network.ba.NOT):
        return ~_bdd_rule(network, rule.args[0], bdd=bdd, variables=variables)

    if isinstance(rule, network.ba.AND):
        value = bdd.true
        for child in rule.args:
            value &= _bdd_rule(network, child, bdd=bdd, variables=variables)
        return value

    if isinstance(rule, network.ba.OR):
        value = bdd.false
        for child in rule.args:
            value |= _bdd_rule(network, child, bdd=bdd, variables=variables)
        return value

    raise TypeError(f"unsupported Boolean expression type: {type(rule)}")


def _reachable_attractors_with_most_permissive_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
) -> Tuple[ConfigurationSet, ...]:
    """Return minimal trap spaces reachable under most-permissive semantics."""

    components = tuple(network.keys())
    initial_hypercube = _validate_initial_hypercube(components, initial_state)
    clingo = network._import_clingo()
    base_program = _most_permissive_asp_base_program(
        network,
        clingo=clingo,
        components=components,
        initial_state=initial_hypercube,
    )
    trapspaces = _solve_reachable_trapspaces_literals(
        clingo,
        base_program=base_program,
    )

    return tuple(
        ConfigurationSet(components, [_hypercube_from_fixed_literals(literals)])
        for literals in trapspaces
    )


def _validate_initial_hypercube(
    components: Tuple[str, ...],
    initial_state: HypercubeLike,
) -> Hypercube:
    """Validate an initial state without enumerating compatible configurations."""

    configurations = ConfigurationSet(components, [initial_state])
    if len(configurations._hypercubes) != 1:
        raise ValueError("invalid initial state representation")

    return configurations._as_hypercubes()[0]


def _most_permissive_asp_base_program(
    network: "BooleanNetwork",
    *,
    clingo: object,
    components: Tuple[str, ...],
    initial_state: Hypercube,
) -> str:
    """Build the ASP program for reachable minimal trap spaces."""

    facts = ["opposite(0,1).", "opposite(1,0)."]

    for component in components:
        facts.append(f"node({_asp_string(clingo, component)}).")

    for component, value in initial_state.items():
        if value.is_fixed:
            facts.append(
                "initial_fixed("
                f"{_asp_string(clingo, component)}, {cast(int, value.value)}"
                ")."
            )

    facts.extend(_boolean_function_asp_facts(network, clingo=clingo))

    return "\n".join(
        [
            *facts,
            _boolean_function_evaluation_asp_rules(),
            """
            timepoint(initial).
            timepoint(final).

            mp_state(initial,N,V) :- initial_fixed(N,V).
            1 { mp_state(initial,N,0); mp_state(initial,N,1) } 1 :-
                node(N), not initial_fixed(N,0), not initial_fixed(N,1).

            1 { final_reach(N,0); final_reach(N,1) } 2 :- node(N).
            mp_reach(final,N,V) :- final_reach(N,V).

            final_reach(N,V) :- mp_eval(final,N,V).

            mp_reach(initial,N,V) :- mp_state(initial,N,V).
            mp_state(final,N,V) :- final_reach(N,V).

            mp_ext(N,V) :- mp_eval(initial,N,V), mp_state(final,N,V).
            { mp_ext(N,V) } :-
                mp_eval(initial,N,V),
                opposite(V,W),
                not mp_state(final,N,V),
                mp_state(final,N,W).
            mp_reach(initial,N,V) :- mp_ext(N,V).

            :- mp_state(final,N,V), not mp_reach(initial,N,V).
            :- mp_state(final,N,V),
               opposite(V,W),
               mp_ext(N,W),
               not mp_ext(N,V).

            #show final_reach/2.
            """,
        ]
    )


def _solve_reachable_trapspaces_literals(
    clingo: object,
    *,
    base_program: str,
) -> Tuple[FrozenSet[Tuple[str, int]], ...]:
    """Enumerate inclusion-minimal reachable trap spaces with domain recursion."""

    control = clingo.Control(  # pyright: ignore[reportAttributeAccessIssue]
        ["--warn=none", "--single-shot"],
        logger=lambda _code, _message: None,
    )
    control.configuration.solve.models = 0
    control.configuration.solve.project = 1
    control.configuration.solve.enum_mode = "domRec"
    control.configuration.solver[0].heuristic = "Domain"
    control.configuration.solver[0].dom_mod = "5,16"
    control.add("base", [], base_program)
    control.ground([("base", [])])

    trapspaces = set()
    with control.solve(yield_=True) as handle:
        for model in handle:
            reachable_values = frozenset(
                (atom.arguments[0].string, atom.arguments[1].number)
                for atom in model.symbols(shown=True)
            )
            trapspaces.add(
                frozenset(
                    (component, value)
                    for component, value in reachable_values
                    if (component, 1 - value) not in reachable_values
                )
            )

    return tuple(sorted(trapspaces, key=_fixed_literals_sort_key))


def _hypercube_from_fixed_literals(
    literals: FrozenSet[Tuple[str, int]],
) -> Hypercube:
    """Build a hypercube from fixed component-value literals."""

    return Hypercube(dict(sorted(literals)))


def _fixed_literals_sort_key(
    literals: FrozenSet[Tuple[str, int]],
) -> Tuple[int, Tuple[Tuple[str, int], ...]]:
    """Return a deterministic sort key for fixed-literal sets."""

    return len(literals), tuple(sorted(literals))


def _hypercube_sort_key(
    hypercube: Hypercube,
) -> Tuple[int, Tuple[Tuple[str, int], ...]]:
    """Return a deterministic sort key for hypercubes."""

    literals = tuple(
        sorted(
            (component, cast(int, value.value))
            for component, value in hypercube.items()
            if value.is_fixed
        )
    )

    return len(literals), literals


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
        if next_state == state:
            return ()

        return (next_state,)

    unstable_mask = state ^ next_state
    if unstable_mask == 0:
        return ()

    if update == "asynchronous":
        return tuple(state ^ bit for bit in _iter_set_bits(unstable_mask))

    return tuple(
        state ^ updated_mask for updated_mask in _iter_nonzero_submasks(unstable_mask)
    )


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


def _iter_set_bits(mask: int) -> Iterable[int]:
    """Yield one-bit masks from low to high bit."""

    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def _iter_nonzero_submasks(mask: int) -> Iterable[int]:
    """Yield all non-zero submasks of a bit mask."""

    submask = mask
    while submask:
        yield submask
        submask = (submask - 1) & mask


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
