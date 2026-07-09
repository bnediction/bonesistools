#!/usr/bin/env python

from __future__ import annotations

from itertools import combinations
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
    HypercubeLike,
    Implicants,
    dnf_implicants,
    prime_implicants,
)

if TYPE_CHECKING:
    from ._network import BooleanNetwork

_StateBits = int
_CompiledRule = Callable[[_StateBits], int]
ImplicantsByValue = Dict[Literal[0, 1], Implicants]


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

    while pending:
        state_bits = pending.pop()
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
            if successor not in successors:
                pending.append(successor)

    terminal_components = _terminal_strongly_connected_components(successors)

    return tuple(
        _configuration_set_from_state_bits(component, components)
        for component in sorted(terminal_components, key=lambda states: min(states))
    )


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
    blockers: List[FrozenSet[Tuple[str, int]]] = []
    trapspaces: Set[FrozenSet[Tuple[str, int]]] = set()

    while True:
        literals = _solve_reachable_trapspace_literals(
            clingo,
            base_program=base_program,
            blockers=blockers,
            components=components,
        )

        if literals is None:
            break

        if literals in trapspaces:
            blockers.append(literals)
            continue

        trapspaces.add(literals)
        blockers.append(literals)

    return tuple(
        ConfigurationSet(components, [_hypercube_from_fixed_literals(literals)])
        for literals in sorted(trapspaces, key=_fixed_literals_sort_key)
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

    implicant_id = 0
    for target, rule in network.items():
        target_symbol = _asp_string(clingo, target)
        for value in (0, 1):
            for implicant in dnf_implicants(rule, value=value, ba=network.ba):
                facts.append(f"implicant({target_symbol}, {value}, {implicant_id}).")
                for source, source_value in implicant.items():
                    if not source_value.is_fixed:
                        raise ValueError(
                            "invalid DNF implicant with free component " f"{source!r}"
                        )

                    facts.append(
                        "implicant_literal("
                        f"{target_symbol}, {value}, {implicant_id}, "
                        f"{_asp_string(clingo, source)}, "
                        f"{cast(int, source_value.value)}"
                        ")."
                    )
                implicant_id += 1

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

            mp_eval(T,N,V) :-
                timepoint(T),
                implicant(N,V,P),
                mp_reach(T,L,W) : implicant_literal(N,V,P,L,W).

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

    n_unstable = unstable_mask.bit_count()
    return tuple(
        state ^ updated_mask
        for size in range(1, n_unstable + 1)
        for updated_mask in _updated_bit_masks(unstable_mask, size, n_components)
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

    bit = 1
    while bit <= mask:
        if mask & bit:
            yield bit
        bit <<= 1


def _updated_bit_masks(
    unstable_mask: int,
    size: int,
    n_components: int,
) -> Iterable[int]:
    """Yield masks selecting a fixed number of unstable components."""

    unstable_bits = [
        1 << index for index in range(n_components) if unstable_mask & (1 << index)
    ]

    for selected_bits in combinations(unstable_bits, size):
        updated_mask = 0
        for bit in selected_bits:
            updated_mask |= bit
        yield updated_mask


def _configuration_set_from_state_bits(
    state_bits: Iterable[_StateBits],
    components: Tuple[str, ...],
) -> ConfigurationSet:
    """Convert explicit state bitsets into a compact ConfigurationSet."""

    configurations = ConfigurationSet(components)
    for state in sorted(state_bits):
        configurations.add(_state_mapping_from_bits(state, components))

    configurations.compress()

    return configurations


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

    index = 0
    stack: List[_StateBits] = []
    on_stack: Set[_StateBits] = set()
    indices: Dict[_StateBits, int] = {}
    lowlinks: Dict[_StateBits, int] = {}
    components: List[Tuple[_StateBits, ...]] = []

    def visit(state: _StateBits) -> None:
        nonlocal index

        indices[state] = index
        lowlinks[state] = index
        index += 1
        stack.append(state)
        on_stack.add(state)

        for successor in successors[state]:
            if successor not in indices:
                visit(successor)
                lowlinks[state] = min(lowlinks[state], lowlinks[successor])
            elif successor in on_stack:
                lowlinks[state] = min(lowlinks[state], indices[successor])

        if lowlinks[state] != indices[state]:
            return None

        component = []
        while True:
            successor = stack.pop()
            on_stack.remove(successor)
            component.append(successor)
            if successor == state:
                break

        components.append(tuple(component))

        return None

    for state in successors:
        if state not in indices:
            visit(state)

    return tuple(components)
