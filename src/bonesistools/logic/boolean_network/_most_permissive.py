#!/usr/bin/env python

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import combinations
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    cast,
)

import clingo

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._robdd import ROBDD
from ..boolean_algebra._structure import _dnf_clauses
from ..boolean_algebra._typing import Configuration, HypercubeLike
from ._dynamics import _validate_complete_hypercube, _validate_hypercube

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _most_permissive_reachability(
    network: "BooleanNetwork",
    source: HypercubeLike,
    target: HypercubeLike,
    *,
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test one-step most-permissive reachability with direct ASP search."""

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

    control = _most_permissive_reachability_asp_control(
        network,
        clingo=clingo,
        components=components,
    )

    if quantifier == "exists":
        return _solve_most_permissive_reachability_asp(
            control,
            clingo=clingo,
            source=source_hypercube,
            target=target_hypercube,
        )

    source_configurations = ConfigurationSet(components, [source_hypercube])
    if quantifier == "robust":
        return all(
            _solve_most_permissive_reachability_asp(
                control,
                clingo=clingo,
                source=Hypercube(configuration),
                target=target_hypercube,
            )
            for configuration in source_configurations
        )

    target_configurations = ConfigurationSet(components, [target_hypercube])
    return all(
        _solve_most_permissive_reachability_asp(
            control,
            clingo=clingo,
            source=Hypercube(source_configuration),
            target=Hypercube(target_configuration),
        )
        for source_configuration in source_configurations
        for target_configuration in target_configurations
    )


def _most_permissive_reachable_attractors(
    network: "BooleanNetwork",
    initial: HypercubeLike,
) -> Tuple[ConfigurationSet, ...]:
    """Return minimal trap spaces reachable under most-permissive semantics."""

    components = tuple(network.keys())
    component_indices = {component: index for index, component in enumerate(components)}
    trap_spaces = _most_permissive_reachable_trap_space_literals(network, initial)

    return tuple(
        _configuration_set_from_fixed_literals(
            components,
            component_indices,
            trap_space,
        )
        for trap_space in trap_spaces
    )


def _most_permissive_reachable_trap_spaces(
    network: "BooleanNetwork",
    initial: HypercubeLike,
) -> Tuple[Hypercube, ...]:
    """Return reachable minimal trap spaces using domain recursion."""

    return tuple(
        _hypercube_from_fixed_literals(literals)
        for literals in _most_permissive_reachable_trap_space_literals(
            network,
            initial,
        )
    )


def _most_permissive_reachable_trap_space_literals(
    network: "BooleanNetwork",
    initial: HypercubeLike,
) -> Tuple[FrozenSet[Tuple[str, int]], ...]:
    """Return fixed literals of reachable minimal trap spaces."""

    components = tuple(network.keys())
    initial_hypercube = _validate_initial_hypercube(components, initial)
    if initial_hypercube:
        base_program = _reachable_trap_spaces_asp_program(
            network,
            clingo=clingo,
            components=components,
            initial_configuration=initial_hypercube,
        )
    else:
        base_program = _trap_spaces_asp_program(
            network,
            clingo=clingo,
            components=components,
        )
    project_models = bool(initial_hypercube) and len(initial_hypercube) < len(
        components
    )
    return _solve_minimal_trap_space_literals(
        clingo,
        base_program=base_program,
        project=project_models,
    )


def _most_permissive_reachability_asp_control(
    network: "BooleanNetwork",
    *,
    clingo: Any,
    components: Tuple[str, ...],
) -> Any:
    """Compile a reusable most-permissive reachability solver."""

    control = clingo.Control(  # pyright: ignore[reportAttributeAccessIssue]
        ["--models=1", "--warn=none"],
        logger=lambda _code, _message: None,
    )
    control.add(
        "base",
        [],
        _most_permissive_reachability_asp_program(
            network,
            clingo=clingo,
            components=components,
        ),
    )
    control.ground([("base", [])])

    return control


def _solve_most_permissive_reachability_asp(
    control: Any,
    *,
    clingo: Any,
    source: Hypercube,
    target: Hypercube,
) -> bool:
    """Solve one existential query with a compiled reachability program."""

    assumptions = [
        *_configuration_asp_assumptions(
            clingo,
            configuration_name="initial",
            configuration=source,
        ),
        *_configuration_asp_assumptions(
            clingo,
            configuration_name="target",
            configuration=target,
        ),
    ]

    return bool(control.solve(assumptions=assumptions).satisfiable)


def _configuration_asp_assumptions(
    clingo: Any,
    *,
    configuration_name: str,
    configuration: Hypercube,
) -> List[Tuple[Any, bool]]:
    """Return assumptions fixing one partial configuration in an ASP query."""

    assumptions = []
    for component, value in configuration.items():
        if not value.is_fixed:
            continue

        atom = clingo.Function(
            "mp_state",
            [
                clingo.Function(configuration_name),
                clingo.String(component),
                clingo.Number(cast(int, value.value)),
            ],
        )
        assumptions.append((atom, True))

    return assumptions


def _most_permissive_reachable_configurations(
    network: "BooleanNetwork",
    initial: HypercubeLike,
) -> Iterator[Configuration]:
    """Iterate over one-step most-permissive reachable configurations."""

    components = tuple(network.keys())
    initial_hypercube = _validate_complete_hypercube(
        components,
        initial,
        name="initial",
    )
    rule_bdds: Dict[str, ROBDD] = {}

    def iterate() -> Iterator[Configuration]:
        for region in _iter_most_permissive_transition_regions(
            network,
            initial_hypercube,
            components=components,
            rule_bdds=rule_bdds,
        ):
            yield from region.iter_configurations(
                initial_hypercube,
                components,
            )

    return iterate()


@dataclass(frozen=True)
class _MPTransitionRegion:
    """Compact encoding of one most-permissive transition region."""

    closure_components: Tuple[str, ...]
    free_components: Tuple[str, ...]
    closed_hypercube: Hypercube
    irreversible_components: Tuple[str, ...]

    def iter_configurations(
        self,
        initial: Hypercube,
        components: Tuple[str, ...],
    ) -> Iterator[Configuration]:
        """Iterate lazily over configurations encoded by the region."""

        optional_components = tuple(
            component
            for component in self.free_components
            if component not in self.irreversible_components
        )
        initial_values = {
            component: cast(int, initial[component].value) for component in components
        }
        required_configuration = dict(initial_values)
        for component in self.irreversible_components:
            required_configuration[component] = 1 - initial_values[component]

        for length in range(len(optional_components) + 1):
            for selected_components in combinations(optional_components, length):
                configuration = dict(required_configuration)
                for component in selected_components:
                    configuration[component] = 1 - initial_values[component]
                yield configuration

    def contains(
        self,
        initial: Hypercube,
        target: Hypercube,
        components: Tuple[str, ...],
    ) -> bool:
        """Test whether the region intersects a target hypercube."""

        for component, target_value in target.items():
            if not target_value.is_fixed:
                continue

            initial_value = cast(int, initial[component].value)
            if component in self.irreversible_components:
                expected_value = 1 - initial_value
            elif component not in self.free_components:
                expected_value = initial_value
            else:
                continue

            if target_value.value != expected_value:
                return False

        return True


class _MostPermissiveTransitionSpace:
    """Compact set of one-step most-permissive reachable configurations."""

    def __init__(
        self,
        *,
        components: Tuple[str, ...],
        initial_configuration: Hypercube,
        regions: Iterable[_MPTransitionRegion],
    ) -> None:

        self._components = components
        self._initial_configuration = initial_configuration
        self._regions = tuple(regions)

    def __iter__(self) -> Iterator[Configuration]:
        """Iterate over represented configurations."""

        for region in self._regions:
            yield from region.iter_configurations(
                self._initial_configuration,
                self._components,
            )

    def contains(
        self,
        target: HypercubeLike,
    ) -> bool:
        """Test whether a target hypercube intersects the represented space."""

        target_hypercube = _validate_hypercube(
            self._components,
            target,
            name="target",
        )

        return any(
            region.contains(
                self._initial_configuration,
                target_hypercube,
                self._components,
            )
            for region in self._regions
        )

    def enumerate(self) -> Tuple[Configuration, ...]:
        """Return all configurations represented by the transition space."""

        return tuple(self)


def _most_permissive_transition_space(
    network: "BooleanNetwork",
    initial: HypercubeLike,
    *,
    rule_bdds: Optional[Dict[str, ROBDD]] = None,
) -> _MostPermissiveTransitionSpace:
    """Return the compact one-step most-permissive transition space."""

    components = tuple(network.keys())
    initial_hypercube = _validate_complete_hypercube(
        components,
        initial,
        name="initial",
    )
    active_rule_bdds = {} if rule_bdds is None else rule_bdds
    regions = _iter_most_permissive_transition_regions(
        network,
        initial_hypercube,
        components=components,
        rule_bdds=active_rule_bdds,
    )

    return _MostPermissiveTransitionSpace(
        components=components,
        initial_configuration=initial_hypercube,
        regions=regions,
    )


def _iter_most_permissive_transition_regions(
    network: "BooleanNetwork",
    initial_configuration: Hypercube,
    *,
    components: Tuple[str, ...],
    rule_bdds: Dict[str, ROBDD],
) -> Iterator[_MPTransitionRegion]:
    """Yield compact most-permissive transition regions."""

    waiting = deque([components])
    scheduled: Set[Tuple[str, ...]] = {components}
    processed: Set[Tuple[str, ...]] = set()

    while waiting:
        closure_components = waiting.popleft()
        scheduled.remove(closure_components)
        if closure_components in processed:
            continue

        processed.add(closure_components)
        closed_hypercube = _smallest_closed_hypercube(
            network,
            initial_configuration,
            relaxed_components=closure_components,
            rule_bdds=rule_bdds,
        )
        free_components = _free_components(components, closed_hypercube)
        irreversible_components = _irreversible_components(
            network,
            initial_configuration,
            closed_hypercube,
            free_components=free_components,
            rule_bdds=rule_bdds,
        )

        yield _MPTransitionRegion(
            closure_components=closure_components,
            free_components=free_components,
            closed_hypercube=closed_hypercube,
            irreversible_components=irreversible_components,
        )

        for next_closure_components in _closure_components_without_irreversibles(
            closure_components,
            irreversible_components,
        ):
            if (
                next_closure_components not in processed
                and next_closure_components not in scheduled
            ):
                waiting.append(next_closure_components)
                scheduled.add(next_closure_components)


def _most_permissive_reachability_asp_program(
    network: "BooleanNetwork",
    *,
    clingo: Any,
    components: Tuple[str, ...],
) -> str:
    """Build an ASP decision program for one-step MP reachability."""

    facts = ["opposite(0,1).", "opposite(1,0)."]
    component_symbols = _asp_component_symbols(clingo, components)

    for component in components:
        facts.append(f"node({component_symbols[component]}).")

    facts.extend(
        _boolean_function_asp_facts(
            network,
            component_symbols=component_symbols,
        )
    )

    return "\n".join(
        [
            *facts,
            _boolean_function_evaluation_asp_rules(),
            """
            timepoint(initial).
            timepoint(target).

            1 { mp_state(initial,N,0); mp_state(initial,N,1) } 1 :- node(N).
            1 { mp_state(target,N,0); mp_state(target,N,1) } 1 :- node(N).

            mp_reach(initial,N,V) :- mp_state(initial,N,V).

            mp_ext(N,V) :- mp_eval(initial,N,V), mp_state(target,N,V).
            { mp_ext(N,V) } :-
                mp_eval(initial,N,V),
                opposite(V,W),
                not mp_state(target,N,V),
                mp_state(target,N,W).

            mp_reach(initial,N,V) :- mp_ext(N,V).

            :- mp_state(target,N,V), not mp_reach(initial,N,V).
            :- mp_state(target,N,V),
               opposite(V,W),
               mp_ext(N,W),
               not mp_ext(N,V).
            """,
        ]
    )


def _validate_initial_hypercube(
    components: Tuple[str, ...],
    initial: HypercubeLike,
) -> Hypercube:
    """Validate an initial hypercube without enumerating configurations."""

    configurations = ConfigurationSet(components, [initial])
    if len(configurations._hypercubes) != 1:
        raise ValueError("invalid initial configuration representation")

    return configurations._as_hypercubes()[0]


def _trap_spaces_asp_program(
    network: "BooleanNetwork",
    *,
    clingo: object,
    components: Tuple[str, ...],
) -> str:
    """Build the ASP program for unrestricted minimal trap spaces."""

    facts = ["opposite(0,1).", "opposite(1,0)."]
    component_symbols = _asp_component_symbols(clingo, components)
    facts.extend(f"node({component_symbols[component]})." for component in components)
    facts.extend(
        _boolean_function_asp_facts(
            network,
            component_symbols=component_symbols,
        )
    )

    return "\n".join(
        [
            *facts,
            _boolean_function_evaluation_asp_rules(),
            _trap_space_asp_rules(),
        ]
    )


def _reachable_trap_spaces_asp_program(
    network: "BooleanNetwork",
    *,
    clingo: object,
    components: Tuple[str, ...],
    initial_configuration: Hypercube,
) -> str:
    """Build the ASP program for reachable minimal trap spaces."""

    facts = ["opposite(0,1).", "opposite(1,0)."]
    component_symbols = _asp_component_symbols(clingo, components)

    for component in components:
        facts.append(f"node({component_symbols[component]}).")

    for component, value in initial_configuration.items():
        if value.is_fixed:
            facts.append(
                "initial_fixed("
                f"{component_symbols[component]}, {cast(int, value.value)}"
                ")."
            )

    facts.extend(
        _boolean_function_asp_facts(
            network,
            component_symbols=component_symbols,
        )
    )

    return "\n".join(
        [
            *facts,
            _boolean_function_evaluation_asp_rules(),
            """
            timepoint(initial).

            mp_state(initial,N,V) :- initial_fixed(N,V).
            1 { mp_state(initial,N,0); mp_state(initial,N,1) } 1 :-
                node(N), not initial_fixed(N,0), not initial_fixed(N,1).

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
            """,
            _trap_space_asp_rules(),
        ]
    )


def _trap_space_asp_rules() -> str:
    """Return closure rules shared by unrestricted and reachable trap spaces."""

    return """
        timepoint(final).

        1 { final_reach(N,0); final_reach(N,1) } 2 :- node(N).
        mp_reach(final,N,V) :- final_reach(N,V).
        final_reach(N,V) :- mp_eval(final,N,V).

        #show final_reach/2.
    """


def _solve_minimal_trap_space_literals(
    clingo: object,
    *,
    base_program: str,
    project: bool,
) -> Tuple[FrozenSet[Tuple[str, int]], ...]:
    """Enumerate inclusion-minimal reachable trap spaces with domain recursion."""

    control = clingo.Control(  # pyright: ignore[reportAttributeAccessIssue]
        ["--warn=none", "--single-shot"],
        logger=lambda _code, _message: None,
    )
    control.configuration.solve.models = 0
    control.configuration.solve.enum_mode = "domRec"
    control.configuration.solver[0].heuristic = "Domain"
    control.configuration.solver[0].dom_mod = f"5,{16 if project else 0}"
    if project:
        control.configuration.solve.project = 1
    control.add("base", [], base_program)
    control.ground([("base", [])])

    trap_spaces = set()
    with control.solve(yield_=True) as handle:
        for model in handle:
            fixed_values: Dict[str, Optional[int]] = {}
            for atom in model.symbols(shown=True):
                component_symbol, value_symbol = atom.arguments
                component = component_symbol.string
                if component in fixed_values:
                    fixed_values[component] = None
                else:
                    fixed_values[component] = value_symbol.number

            trap_spaces.add(
                frozenset(
                    (component, value)
                    for component, value in fixed_values.items()
                    if value is not None
                )
            )

    return tuple(sorted(trap_spaces, key=_fixed_literals_sort_key))


def _hypercube_from_fixed_literals(
    literals: FrozenSet[Tuple[str, int]],
) -> Hypercube:
    """Build a hypercube from fixed component-value literals."""

    return Hypercube(dict(sorted(literals)))


def _configuration_set_from_fixed_literals(
    components: Tuple[str, ...],
    component_indices: Dict[str, int],
    literals: FrozenSet[Tuple[str, int]],
) -> ConfigurationSet:
    """Build a configuration set directly from fixed-literal bitsets."""

    fixed_mask = 0
    value_mask = 0
    for component, value in literals:
        bit = 1 << component_indices[component]
        fixed_mask |= bit
        if value == 1:
            value_mask |= bit

    return ConfigurationSet._from_encoded_hypercubes(
        components,
        [(fixed_mask, value_mask)],
    )


def _fixed_literals_sort_key(
    literals: FrozenSet[Tuple[str, int]],
) -> Tuple[int, Tuple[Tuple[str, int], ...]]:
    """Return a deterministic sort key for fixed-literal sets."""

    return len(literals), tuple(sorted(literals))


def _free_components(
    components: Tuple[str, ...],
    closed_hypercube: Hypercube,
) -> Tuple[str, ...]:
    """Return components free in a closed hypercube."""

    return tuple(
        component
        for component in components
        if component not in closed_hypercube or not closed_hypercube[component].is_fixed
    )


def _irreversible_components(
    network: "BooleanNetwork",
    initial_configuration: Hypercube,
    closed_hypercube: Hypercube,
    *,
    free_components: Tuple[str, ...],
    rule_bdds: Dict[str, ROBDD],
) -> Tuple[str, ...]:
    """Return free components forced to leave their initial value."""

    irreversible_components = []
    for component in free_components:
        initial_value = cast(Literal[0, 1], initial_configuration[component].value)
        if not _rule_can_take_value_in_hypercube(
            network,
            component,
            value=initial_value,
            hypercube=closed_hypercube,
            rule_bdds=rule_bdds,
        ):
            irreversible_components.append(component)

    return tuple(irreversible_components)


def _closure_components_without_irreversibles(
    closure_components: Tuple[str, ...],
    irreversible_components: Tuple[str, ...],
) -> Iterable[Tuple[str, ...]]:
    """Yield closure-component sets obtained by removing irreversibles."""

    for length in range(1, len(irreversible_components) + 1):
        for removed in combinations(irreversible_components, length):
            removed_set = set(removed)
            yield tuple(
                component
                for component in closure_components
                if component not in removed_set
            )


def _smallest_closed_hypercube(
    network: "BooleanNetwork",
    initial_configuration: HypercubeLike,
    *,
    relaxed_components: Iterable[str],
    depth: Optional[int] = None,
    rule_bdds: Optional[Dict[str, ROBDD]] = None,
) -> Hypercube:
    """Return the smallest closed hypercube obtained by relaxing components in K."""

    hypercube = Hypercube(initial_configuration)
    fixed = {
        component
        for component in relaxed_components
        if component in hypercube and hypercube[component].is_fixed
    }
    if rule_bdds is None:
        rule_bdds = {}

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
                rule_bdds=rule_bdds,
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
    rule_bdds: Dict[str, ROBDD],
) -> bool:
    """Test whether a rule can take a value inside a hypercube."""

    if component not in rule_bdds:
        rule_bdds[component] = ROBDD._from_rule(
            network[component],
            ba=network.ba,
        )

    robdd = rule_bdds[component]
    fixed_values = {
        source: cast(int, hypercube[source].value)
        for source in robdd.variables
        if source in hypercube and hypercube[source].is_fixed
    }
    return robdd._can_reach_terminal(fixed_values, value)


def _boolean_function_asp_facts(
    network: "BooleanNetwork",
    *,
    component_symbols: Dict[str, str],
) -> List[str]:
    """Encode unate rules as DNF and non-unate rules as exact ROBDDs."""

    facts = []

    for target, rule in network.items():
        target_symbol = component_symbols[target]
        if _expression_is_unate(network, rule):
            facts.append(f"unate({target_symbol}).")
            facts.extend(
                _positive_dnf_asp_facts(
                    network,
                    target=target,
                    component_symbols=component_symbols,
                )
            )
        else:
            facts.extend(
                _non_unate_bdd_asp_facts(
                    network,
                    target=target,
                    component_symbols=component_symbols,
                )
            )

    return facts


def _boolean_function_evaluation_asp_rules() -> str:
    """Return shared MP evaluation rules for DNF and ROBDD encodings."""

    return """
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
        mp_eval(T,N,0) :-
            timepoint(T), unate(N), not has_positive_clause(N).
        mp_eval(T,N,0) :-
            timepoint(T),
            unate(N),
            has_positive_clause(N),
            positive_clause_falsified(T,N,C) : positive_clause(N,C).

        bdd_eval(T,N,0,0) :- timepoint(T), non_unate(N).
        bdd_eval(T,N,1,1) :- timepoint(T), non_unate(N).
        bdd_eval(T,N,B,V) :-
            bdd_node(N,B,L,LO,HI),
            mp_reach(T,L,0),
            bdd_eval(T,N,LO,V).
        bdd_eval(T,N,B,V) :-
            bdd_node(N,B,L,LO,HI),
            mp_reach(T,L,1),
            bdd_eval(T,N,HI,V).
        mp_eval(T,N,V) :-
            non_unate(N),
            bdd_root(N,B),
            bdd_eval(T,N,B,V).
    """


def _asp_string(clingo: object, value: str) -> str:
    """Return a clingo string literal."""

    return str(clingo.String(value))  # pyright: ignore[reportAttributeAccessIssue]


def _asp_component_symbols(
    clingo: object,
    components: Tuple[str, ...],
) -> Dict[str, str]:
    """Encode component names once for reuse across ASP facts."""

    return {component: _asp_string(clingo, component) for component in components}


def _expression_is_unate(network: "BooleanNetwork", rule: Any) -> bool:
    """Return whether every symbol occurs with only one polarity."""

    polarities: Dict[str, Set[int]] = {}

    for literal in rule.literalize().get_literals():
        if isinstance(literal, network.ba.NOT):
            component = str(literal.args[0].obj)
            value = 0
        else:
            component = str(literal.obj)
            value = 1

        polarities.setdefault(component, set()).add(value)

    return all(len(values) == 1 for values in polarities.values())


def _positive_dnf_asp_facts(
    network: "BooleanNetwork",
    *,
    target: str,
    component_symbols: Dict[str, str],
) -> List[str]:
    """Encode positive DNF clauses for one unate function."""

    target_symbol = component_symbols[target]
    facts = []
    clauses = sorted(
        (
            tuple(sorted(clause))
            for clause in _dnf_clauses(
                network.ba,
                network[target].literalize(),
            )
        ),
        key=lambda clause: (len(clause), clause),
    )

    for clause_id, clause in enumerate(clauses):
        facts.append(f"positive_clause({target_symbol}, {clause_id}).")
        for source, source_value in clause:
            facts.append(
                "positive_clause_literal("
                f"{target_symbol}, {clause_id}, "
                f"{component_symbols[source]}, "
                f"{source_value}"
                ")."
            )

    return facts


def _non_unate_bdd_asp_facts(
    network: "BooleanNetwork",
    *,
    target: str,
    component_symbols: Dict[str, str],
) -> List[str]:
    """Encode one non-unate Boolean function as a reduced ordered BDD."""

    rule = network[target]
    robdd = ROBDD._from_rule(rule, ba=network.ba)
    target_symbol = component_symbols[target]
    facts = [
        f"non_unate({target_symbol}).",
        f"bdd_root({target_symbol}, {robdd._root_id}).",
    ]

    for node, source, low, high in robdd._iter_nodes():
        facts.append(
            "bdd_node("
            f"{target_symbol}, {node}, {component_symbols[source]}, "
            f"{low}, {high}"
            ")."
        )

    return facts
