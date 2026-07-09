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
from ..boolean_algebra import (
    ConfigurationSet,
    Hypercube,
    dnf_implicants,
    prime_implicants,
)
from ..boolean_algebra._structure import Implicants
from ..boolean_algebra._typing import HypercubeLike

if TYPE_CHECKING:
    from ._network import BooleanNetwork

_StateBits = int
_CompiledRule = Callable[[_StateBits], int]
ImplicantsByValue = Dict[Literal[0, 1], Implicants]


def _reachable_attractors_with_bdd_backend(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Compute reachable attractors through a BDD reachability set."""

    bdd_module = _import_dd_autoref()
    bdd = bdd_module.BDD()
    components = tuple(network.keys())
    current_vars = tuple(f"x{index}" for index in range(len(components)))
    next_vars = tuple(f"y{index}" for index in range(len(components)))

    ordered_vars = tuple(
        variable for pair in zip(current_vars, next_vars) for variable in pair
    )
    bdd.declare(*ordered_vars)

    transition = _bdd_transition_relation(
        network,
        bdd=bdd,
        components=components,
        current_vars=current_vars,
        next_vars=next_vars,
        update=update,
    )
    initial = _bdd_initial_states(
        bdd=bdd,
        components=components,
        current_vars=current_vars,
        initial_states=initial_states,
    )
    reachable = _bdd_reachable_states(
        bdd=bdd,
        initial=initial,
        transition=transition,
        current_vars=current_vars,
        next_vars=next_vars,
    )
    reachable_state_bits = _bdd_state_bits(
        bdd=bdd,
        states=reachable,
        current_vars=current_vars,
    )

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


def _import_dd_autoref() -> Any:
    """Import the optional BDD backend."""

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
    update: Literal["asynchronous", "synchronous", "general"],
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

    if update == "asynchronous":
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

    relation = bdd.true
    updates = bdd.false
    for current_value, next_value, rule in zip(current_values, next_values, rules):
        updated = _bdd_exclusive_or(current_value, rule)
        updated &= _bdd_equivalence(next_value, rule)
        unchanged = _bdd_equivalence(next_value, current_value)
        relation &= unchanged | updated
        updates |= updated

    return relation & updates


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


def _bdd_initial_states(
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    initial_states: ConfigurationSet,
) -> Any:
    """Encode an initial ConfigurationSet as a BDD."""

    states = bdd.false
    for hypercube in initial_states._as_hypercubes():
        state = bdd.true
        for component, current_var in zip(components, current_vars):
            if component not in hypercube:
                continue

            value = hypercube[component]
            if value.is_fixed:
                variable = bdd.var(current_var)
                state &= variable if value.value else ~variable

        states |= state

    return states


def _bdd_reachable_states(
    *,
    bdd: Any,
    initial: Any,
    transition: Any,
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Any:
    """Compute the reachable-state set as a BDD fixpoint."""

    reachable = initial
    frontier = initial
    rename_next_to_current = dict(zip(next_vars, current_vars))

    while frontier != bdd.false:
        successors = bdd.exist(current_vars, frontier & transition)
        successors = bdd.let(rename_next_to_current, successors)
        new_frontier = successors & ~reachable
        if new_frontier == bdd.false:
            break

        reachable |= new_frontier
        frontier = new_frontier

    return reachable


def _bdd_state_bits(
    *,
    bdd: Any,
    states: Any,
    current_vars: Tuple[str, ...],
) -> Tuple[_StateBits, ...]:
    """Decode a BDD state set into sorted bitsets."""

    decoded = []
    for assignment in bdd.pick_iter(states, care_vars=current_vars):
        state_bits = 0
        for index, current_var in enumerate(current_vars):
            if assignment.get(current_var, False):
                state_bits |= 1 << index
        decoded.append(state_bits)

    return tuple(sorted(decoded))


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
        components=components,
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

    clause_id = 0
    for target, rule in network.items():
        target_symbol = _asp_string(clingo, target)
        for implicant in dnf_implicants(rule, value=1, ba=network.ba):
            facts.append(f"positive_clause({target_symbol}, {clause_id}).")
            for source, source_value in implicant.items():
                if not source_value.is_fixed:
                    raise ValueError(
                        "invalid DNF implicant with free component " f"{source!r}"
                    )

                facts.append(
                    "positive_clause_literal("
                    f"{target_symbol}, {clause_id}, "
                    f"{_asp_string(clingo, source)}, "
                    f"{cast(int, source_value.value)}"
                    ")."
                )
            clause_id += 1

    return "\n".join(
        [
            *facts,
            """
            timepoint(initial).
            timepoint(final).

            mp_state(initial,N,V) :- initial_fixed(N,V).
            1 { mp_state(initial,N,0); mp_state(initial,N,1) } 1 :-
                node(N), not initial_fixed(N,0), not initial_fixed(N,1).

            1 { final_reach(N,0); final_reach(N,1) } 2 :- node(N).
            mp_reach(final,N,V) :- final_reach(N,V).

            has_positive_clause(N) :- positive_clause(N,C).

            positive_clause_satisfied(T,N,C) :-
                timepoint(T),
                positive_clause(N,C),
                mp_reach(T,L,W) : positive_clause_literal(N,C,L,W).
            positive_clause_falsified(T,N,C) :-
                timepoint(T),
                positive_clause_literal(N,C,L,W),
                opposite(W,V),
                mp_reach(T,L,V).

            mp_eval(T,N,1) :- positive_clause_satisfied(T,N,C).
            mp_eval(T,N,0) :- timepoint(T), node(N), not has_positive_clause(N).
            mp_eval(T,N,0) :-
                timepoint(T),
                has_positive_clause(N),
                positive_clause_falsified(T,N,C) : positive_clause(N,C).

            final_reach(N,V) :- mp_eval(final,N,V).
            fixed(N,V) :-
                final_reach(N,V),
                opposite(V,W),
                not final_reach(N,W).

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

            #maximize { 1,N,V : fixed(N,V) }.
            #show fixed/2.
            """,
        ]
    )


def _solve_reachable_trapspace_literals(
    clingo: object,
    *,
    base_program: str,
    blockers: List[FrozenSet[Tuple[str, int]]],
    components: Tuple[str, ...],
) -> Optional[FrozenSet[Tuple[str, int]]]:
    """Solve one reachable minimal trap-space candidate."""

    control = clingo.Control(  # pyright: ignore[reportAttributeAccessIssue]
        ["--opt-mode=opt", "--opt-strategy=usc", "--warn=none"]
    )
    control.add(
        "base",
        [],
        "\n".join(
            [
                base_program,
                _trapspace_blocker_program(
                    clingo,
                    blockers=blockers,
                    components=components,
                ),
            ]
        ),
    )
    control.ground([("base", [])])

    optimal_literals = None
    with control.solve(yield_=True) as handle:
        for model in handle:
            optimal_literals = frozenset(
                (atom.arguments[0].string, atom.arguments[1].number)
                for atom in model.symbols(shown=True)
            )

    return optimal_literals


def _solve_reachable_trapspaces_literals(
    clingo: object,
    *,
    base_program: str,
    components: Tuple[str, ...],
) -> Tuple[FrozenSet[Tuple[str, int]], ...]:
    """Enumerate reachable minimal trap spaces with incremental blockers."""

    control = clingo.Control(  # pyright: ignore[reportAttributeAccessIssue]
        ["--opt-mode=opt", "--opt-strategy=usc", "--warn=none"]
    )
    control.add("base", [], base_program)
    control.ground([("base", [])])

    blocker_id = 0
    trapspaces: List[FrozenSet[Tuple[str, int]]] = []
    seen: Set[FrozenSet[Tuple[str, int]]] = set()

    while True:
        optimal_literals = None
        with control.solve(yield_=True) as handle:
            for model in handle:
                optimal_literals = frozenset(
                    (atom.arguments[0].string, atom.arguments[1].number)
                    for atom in model.symbols(shown=True)
                )

        if optimal_literals is None:
            break

        if optimal_literals not in seen:
            trapspaces.append(optimal_literals)
            seen.add(optimal_literals)

        blocker_part = f"blocker_{blocker_id}"
        control.add(
            blocker_part,
            [],
            _incremental_trapspace_blocker_program(
                clingo,
                blocker_id=blocker_id,
                literals=optimal_literals,
                components=components,
            ),
        )
        control.ground([(blocker_part, [])])
        blocker_id += 1

    return tuple(sorted(trapspaces, key=_fixed_literals_sort_key))


def _incremental_trapspace_blocker_program(
    clingo: object,
    *,
    blocker_id: int,
    literals: FrozenSet[Tuple[str, int]],
    components: Tuple[str, ...],
) -> str:
    """Return one blocker constraint for an already enumerated trap space."""

    outside = f"outside_{blocker_id}"
    lines = []

    for component in components:
        for value in (0, 1):
            if (component, value) not in literals:
                lines.append(
                    f"{outside} :- fixed({_asp_string(clingo, component)}, {value})."
                )

    lines.append(f":- not {outside}.")
    return "\n".join(lines)


def _trapspace_blocker_program(
    clingo: object,
    *,
    blockers: List[FrozenSet[Tuple[str, int]]],
    components: Tuple[str, ...],
) -> str:
    """Return ASP constraints excluding already enumerated trap spaces."""

    lines = []
    for blocker_id, literals in enumerate(blockers):
        lines.append(f"blocker({blocker_id}).")

        for component in components:
            for value in (0, 1):
                if (component, value) not in literals:
                    lines.append(
                        "outside_literal("
                        f"{blocker_id}, {_asp_string(clingo, component)}, {value}"
                        ")."
                    )

    if blockers:
        lines.extend(
            [
                "outside(B) :- outside_literal(B,N,V), fixed(N,V).",
                ":- blocker(B), not outside(B).",
            ]
        )

    return "\n".join(lines)


def _asp_string(clingo: object, value: str) -> str:
    """Return a clingo string literal."""

    return str(clingo.String(value))  # pyright: ignore[reportAttributeAccessIssue]


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


def _smallest_closed_hypercube(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    *,
    components_to_relax: Iterable[str],
    depth: Optional[int] = None,
) -> Hypercube:
    """Return the smallest closed hypercube obtained by relaxing components in K."""

    hypercube = Hypercube(initial_state)
    fixed = {
        component
        for component in components_to_relax
        if component in hypercube and hypercube[component].is_fixed
    }
    implicants: Dict[str, ImplicantsByValue] = {}

    iteration = 0
    while fixed and (depth is None or iteration < depth):
        changed = False
        current_hypercube = hypercube.copy() if depth is not None else hypercube

        for component in tuple(fixed):
            current_value = hypercube[component]
            if not current_value.is_fixed:
                fixed.remove(component)
                continue

            current_fixed_value = cast(Literal[0, 1], current_value.value)
            target_value = cast(Literal[0, 1], 1 - current_fixed_value)
            if _rule_can_take_value_in_hypercube(
                network,
                component,
                value=target_value,
                hypercube=current_hypercube,
                implicants=implicants,
            ):
                hypercube.pop(component, None)
                fixed.remove(component)
                changed = True

        if not changed:
            break

        iteration += 1

    return hypercube


def _rule_can_take_value_in_hypercube(
    network: "BooleanNetwork",
    component: str,
    *,
    value: Literal[0, 1],
    hypercube: Hypercube,
    implicants: Dict[str, ImplicantsByValue],
) -> bool:
    """Test whether a rule can take a value inside a hypercube."""

    if component not in implicants:
        implicants[component] = {}

    if value not in implicants[component]:
        implicants[component][value] = prime_implicants(
            network[component],
            value=value,
            backend="asp",
            ba=network.ba,
        )

    return any(
        _implicant_intersects_hypercube(implicant, hypercube)
        for implicant in implicants[component][value]
    )


def _implicant_intersects_hypercube(
    implicant: Hypercube,
    hypercube: Hypercube,
) -> bool:
    """Test whether an implicant has at least one state in a hypercube."""

    for component, value in implicant.items():
        if component in hypercube:
            hypercube_value = hypercube[component]
            if hypercube_value.is_fixed and hypercube_value != value:
                return False

    return True


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
        state ^ updated_mask
        for updated_mask in _iter_nonzero_submasks(unstable_mask)
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
