#!/usr/bin/env python

from __future__ import annotations

from itertools import combinations
from typing import (
    TYPE_CHECKING,
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

_StateKey = Tuple[int, ...]
ImplicantsByValue = Dict[Literal[0, 1], Implicants]


def _reachable_attractors_with_explicit_backend(
    network: "BooleanNetwork",
    initial_states: ConfigurationSet,
    *,
    update: Literal["asynchronous", "synchronous", "general"],
) -> Tuple[ConfigurationSet, ...]:
    """Enumerate terminal SCCs reachable in the explicit state-transition graph."""

    components = tuple(network.keys())
    successors: Dict[_StateKey, Tuple[_StateKey, ...]] = {}
    pending = [
        tuple(initial_state[component] for component in components)
        for initial_state in initial_states
    ]

    while pending:
        state_key = pending.pop()
        if state_key in successors:
            continue

        state = dict(zip(components, state_key))
        state_successors = _successor_state_keys(
            network,
            state,
            update=update,
            components=components,
        )
        successors[state_key] = state_successors

        for successor in state_successors:
            if successor not in successors:
                pending.append(successor)

    terminal_components = _terminal_strongly_connected_components(successors)

    return tuple(
        _configuration_set_from_state_keys(component, components)
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

    return configurations._hypercubes[0]


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
                            "invalid DNF implicant with free component "
                            f"{source!r}"
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
    literals: FrozenSet[Tuple[str, int]]
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


def _successor_state_keys(
    network: "BooleanNetwork",
    state: Dict[str, int],
    *,
    update: Literal["asynchronous", "synchronous", "general"],
    components: Tuple[str, ...],
) -> Tuple[_StateKey, ...]:
    """Return successor state keys under a finite-state update semantics."""

    next_state = network.next_configuration(state)

    if update == "synchronous":
        if next_state == state:
            return ()

        return (tuple(next_state[component] for component in components),)

    unstable = tuple(
        component
        for component in components
        if next_state[component] != state[component]
    )

    if not unstable:
        return ()

    if update == "asynchronous":
        return tuple(
            _updated_state_key(
                state,
                next_state,
                updated_components=(component,),
                components=components,
            )
            for component in unstable
        )

    return tuple(
        _updated_state_key(
            state,
            next_state,
            updated_components=updated_components,
            components=components,
        )
        for size in range(1, len(unstable) + 1)
        for updated_components in combinations(unstable, size)
    )


def _updated_state_key(
    state: Mapping[str, int],
    next_state: Mapping[str, int],
    *,
    updated_components: Tuple[str, ...],
    components: Tuple[str, ...],
) -> _StateKey:
    """Return the key obtained by updating selected components."""

    updated = dict(state)
    for component in updated_components:
        updated[component] = next_state[component]

    return tuple(updated[component] for component in components)


def _configuration_set_from_state_keys(
    state_keys: Iterable[_StateKey],
    components: Tuple[str, ...],
) -> ConfigurationSet:
    """Convert explicit state keys into a compact ConfigurationSet."""

    configurations = ConfigurationSet(components)
    for state_key in sorted(state_keys):
        configurations.add(dict(zip(components, state_key)))

    configurations.compress()

    return configurations


def _terminal_strongly_connected_components(
    successors: Mapping[_StateKey, Tuple[_StateKey, ...]]
) -> Tuple[Tuple[_StateKey, ...], ...]:
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
    successors: Mapping[_StateKey, Tuple[_StateKey, ...]]
) -> Tuple[Tuple[_StateKey, ...], ...]:
    """Compute SCCs of a directed graph encoded as a successor mapping."""

    index = 0
    stack: List[_StateKey] = []
    on_stack: Set[_StateKey] = set()
    indices: Dict[_StateKey, int] = {}
    lowlinks: Dict[_StateKey, int] = {}
    components: List[Tuple[_StateKey, ...]] = []

    def visit(state: _StateKey) -> None:
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
