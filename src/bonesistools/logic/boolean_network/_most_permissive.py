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

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube, dnf_implicants
from ..boolean_algebra._robdd import ROBDD
from ..boolean_algebra._typing import HypercubeLike
from ._dynamics import _validate_complete_hypercube, _validate_hypercube

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _most_permissive_reachability(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
    *,
    quantifier: Literal["exists", "robust", "universal"],
) -> bool:
    """Test one-step most-permissive reachability with direct ASP search."""

    components = tuple(network.keys())
    initial_hypercube = _validate_hypercube(
        components,
        initial_state,
        name="initial_state",
    )
    target_hypercube = _validate_hypercube(
        components,
        target_state,
        name="target_state",
    )

    clingo = network._import_clingo()
    control = _most_permissive_reachability_asp_control(
        network,
        clingo=clingo,
        components=components,
    )

    if quantifier == "exists":
        return _solve_most_permissive_reachability_asp(
            control,
            clingo=clingo,
            initial_state=initial_hypercube,
            target_state=target_hypercube,
        )

    initial_configurations = ConfigurationSet(components, [initial_hypercube])
    if quantifier == "robust":
        return all(
            _solve_most_permissive_reachability_asp(
                control,
                clingo=clingo,
                initial_state=Hypercube(configuration),
                target_state=target_hypercube,
            )
            for configuration in initial_configurations
        )

    target_configurations = ConfigurationSet(components, [target_hypercube])
    return all(
        _solve_most_permissive_reachability_asp(
            control,
            clingo=clingo,
            initial_state=Hypercube(initial_configuration),
            target_state=Hypercube(target_configuration),
        )
        for initial_configuration in initial_configurations
        for target_configuration in target_configurations
    )


def _most_permissive_reachable_attractors(
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
    initial_state: Hypercube,
    target_state: Hypercube,
) -> bool:
    """Solve one existential query with a compiled reachability program."""

    assumptions = [
        *_state_asp_assumptions(
            clingo,
            state_name="initial",
            state=initial_state,
        ),
        *_state_asp_assumptions(
            clingo,
            state_name="target",
            state=target_state,
        ),
    ]

    return bool(control.solve(assumptions=assumptions).satisfiable)


def _state_asp_assumptions(
    clingo: Any,
    *,
    state_name: str,
    state: Hypercube,
) -> List[Tuple[Any, bool]]:
    """Return assumptions fixing one partial state in an ASP query."""

    assumptions = []
    for component, value in state.items():
        if not value.is_fixed:
            continue

        atom = clingo.Function(
            "mp_state",
            [
                clingo.Function(state_name),
                clingo.String(component),
                clingo.Number(cast(int, value.value)),
            ],
        )
        assumptions.append((atom, True))

    return assumptions


def _most_permissive_reachable_configurations(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
) -> Iterator[Dict[str, int]]:
    """Iterate over one-step most-permissive reachable configurations."""

    components = tuple(network.keys())
    initial_hypercube = _validate_complete_hypercube(
        components,
        initial_state,
        name="initial_state",
    )
    rule_bdds: Dict[str, ROBDD] = {}

    def iterate() -> Iterator[Dict[str, int]]:
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
        initial_state: Hypercube,
        components: Tuple[str, ...],
    ) -> Iterator[Dict[str, int]]:
        """Iterate lazily over configurations encoded by the region."""

        optional_components = tuple(
            component
            for component in self.free_components
            if component not in self.irreversible_components
        )
        initial_values = {
            component: cast(int, initial_state[component].value)
            for component in components
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
        initial_state: Hypercube,
        target_state: Hypercube,
        components: Tuple[str, ...],
    ) -> bool:
        """Test whether the region intersects a target hypercube."""

        for component, target_value in target_state.items():
            if not target_value.is_fixed:
                continue

            initial_value = cast(int, initial_state[component].value)
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
        initial_state: Hypercube,
        regions: Iterable[_MPTransitionRegion],
    ) -> None:

        self._components = components
        self._initial_state = initial_state
        self._regions = tuple(regions)

    def __iter__(self) -> Iterator[Dict[str, int]]:
        """Iterate over represented configurations."""

        for region in self._regions:
            yield from region.iter_configurations(
                self._initial_state,
                self._components,
            )

    def contains(
        self,
        target_state: HypercubeLike,
    ) -> bool:
        """Test whether a target hypercube intersects the represented space."""

        target_hypercube = _validate_hypercube(
            self._components,
            target_state,
            name="target_state",
        )

        return any(
            region.contains(
                self._initial_state,
                target_hypercube,
                self._components,
            )
            for region in self._regions
        )

    def enumerate(self) -> Tuple[Dict[str, int], ...]:
        """Return all configurations represented by the transition space."""

        return tuple(self)


def _most_permissive_transition_space(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    *,
    rule_bdds: Optional[Dict[str, ROBDD]] = None,
) -> _MostPermissiveTransitionSpace:
    """Return the compact one-step most-permissive transition space."""

    components = tuple(network.keys())
    initial_hypercube = _validate_complete_hypercube(
        components,
        initial_state,
        name="initial_state",
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
        initial_state=initial_hypercube,
        regions=regions,
    )


def _iter_most_permissive_transition_regions(
    network: "BooleanNetwork",
    initial_state: Hypercube,
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
            initial_state,
            relaxed_components=closure_components,
            rule_bdds=rule_bdds,
        )
        free_components = _free_components(components, closed_hypercube)
        irreversible_components = _irreversible_components(
            network,
            initial_state,
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

    for component in components:
        component_symbol = _asp_string(clingo, component)
        facts.append(f"node({component_symbol}).")

    facts.extend(_boolean_function_asp_facts(network, clingo=clingo))

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
    initial_state: Hypercube,
    closed_hypercube: Hypercube,
    *,
    free_components: Tuple[str, ...],
    rule_bdds: Dict[str, ROBDD],
) -> Tuple[str, ...]:
    """Return free components forced to leave their initial value."""

    irreversible_components = []
    for component in free_components:
        initial_value = cast(Literal[0, 1], initial_state[component].value)
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
    initial_state: HypercubeLike,
    *,
    relaxed_components: Iterable[str],
    depth: Optional[int] = None,
    rule_bdds: Optional[Dict[str, ROBDD]] = None,
) -> Hypercube:
    """Return the smallest closed hypercube obtained by relaxing components in K."""

    hypercube = Hypercube(initial_state)
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
    clingo: object,
) -> List[str]:
    """Encode unate rules as DNF and non-unate rules as exact ROBDDs."""

    facts = []

    for target, rule in network.items():
        target_symbol = _asp_string(clingo, target)
        if _expression_is_unate(network, rule):
            facts.append(f"unate({target_symbol}).")
            facts.extend(
                _positive_dnf_asp_facts(
                    network,
                    target=target,
                    clingo=clingo,
                )
            )
        else:
            facts.extend(
                _non_unate_bdd_asp_facts(
                    network,
                    target=target,
                    clingo=clingo,
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
    clingo: object,
) -> List[str]:
    """Encode positive DNF clauses for one unate function."""

    target_symbol = _asp_string(clingo, target)
    facts = []

    for clause_id, implicant in enumerate(
        dnf_implicants(network[target], value=1, ba=network.ba)
    ):
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

    return facts


def _non_unate_bdd_asp_facts(
    network: "BooleanNetwork",
    *,
    target: str,
    clingo: object,
) -> List[str]:
    """Encode one non-unate Boolean function as a reduced ordered BDD."""

    rule = network[target]
    robdd = ROBDD._from_rule(rule, ba=network.ba)
    target_symbol = _asp_string(clingo, target)
    facts = [
        f"non_unate({target_symbol}).",
        f"bdd_root({target_symbol}, {robdd._root_id}).",
    ]

    for node, source, low, high in robdd._iter_nodes():
        facts.append(
            "bdd_node("
            f"{target_symbol}, {node}, {_asp_string(clingo, source)}, "
            f"{low}, {high}"
            ")."
        )

    return facts
