#!/usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
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

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _reachability_with_most_permissive_hypercube_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
) -> bool:
    """Test existential MP reachability with hypercube regions."""

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

    initial_configurations = ConfigurationSet(components, [initial_hypercube])
    rule_bdds: Dict[str, ROBDD] = {}
    reachability = (
        _most_permissive_transition_space(
            network,
            configuration,
            rule_bdds=rule_bdds,
        ).contains(target_hypercube)
        for configuration in initial_configurations
    )
    return any(reachability)


def _reachability_with_most_permissive_asp_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
    *,
    quantifier: Literal["exists", "robust"],
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
    return all(
        _solve_most_permissive_reachability_asp(
            control,
            clingo=clingo,
            initial_state=Hypercube(configuration),
            target_state=target_hypercube,
        )
        for configuration in initial_configurations
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


def _reachable_configurations_with_most_permissive_hypercube_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
) -> Iterator[Dict[str, int]]:
    """Iterate over one-step most-permissive reachable configurations."""

    return iter(_most_permissive_transition_space(network, initial_state))


@dataclass(frozen=True)
class _MPTransitionRegion:
    """Compact encoding of one most-permissive transition region."""

    closure_components: Tuple[str, ...]
    free_components: Tuple[str, ...]
    closed_hypercube: Hypercube
    irreversible_components: Tuple[str, ...]

    def enumerate(
        self,
        initial_state: Hypercube,
        components: Tuple[str, ...],
    ) -> Tuple[Dict[str, int], ...]:
        """Return all configurations encoded by the region."""

        optional_components = tuple(
            component
            for component in self.free_components
            if component not in self.irreversible_components
        )
        configurations: List[Dict[str, int]] = []

        for length in range(len(optional_components) + 1):
            for selected_components in combinations(optional_components, length):
                flipped_components = set(self.irreversible_components)
                flipped_components.update(selected_components)

                configuration = {}
                for component in components:
                    initial_value = cast(int, initial_state[component].value)
                    if component in flipped_components:
                        configuration[component] = 1 - initial_value
                    else:
                        configuration[component] = initial_value

                configurations.append(configuration)

        return tuple(configurations)

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
            yield from region.enumerate(self._initial_state, self._components)

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
    if rule_bdds is None:
        rule_bdds = {}

    regions: List[_MPTransitionRegion] = []
    waiting: List[Tuple[str, ...]] = [components]
    scheduled: Set[Tuple[str, ...]] = {components}
    processed: Set[Tuple[str, ...]] = set()

    while waiting:
        closure_components = waiting.pop(0)
        scheduled.remove(closure_components)
        if closure_components in processed:
            continue

        processed.add(closure_components)
        closed_hypercube = _smallest_closed_hypercube(
            network,
            initial_hypercube,
            relaxed_components=closure_components,
            rule_bdds=rule_bdds,
        )
        free_components = _free_components(components, closed_hypercube)
        irreversible_components = _irreversible_components(
            network,
            initial_hypercube,
            closed_hypercube,
            free_components=free_components,
            rule_bdds=rule_bdds,
        )

        regions.append(
            _MPTransitionRegion(
                closure_components=closure_components,
                free_components=free_components,
                closed_hypercube=closed_hypercube,
                irreversible_components=irreversible_components,
            )
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

    return _MostPermissiveTransitionSpace(
        components=components,
        initial_state=initial_hypercube,
        regions=regions,
    )


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


def _validate_complete_hypercube(
    components: Tuple[str, ...],
    state: HypercubeLike,
    *,
    name: str,
) -> Hypercube:
    """Validate a complete Boolean configuration encoded as a hypercube."""

    hypercube = Hypercube(state)
    component_set = set(components)
    unknown_components = sorted(set(hypercube) - component_set)
    if unknown_components:
        raise ValueError(f"unknown components in {name}: {unknown_components}")

    invalid_components = [
        component
        for component in components
        if component not in hypercube or not hypercube[component].is_fixed
    ]
    if invalid_components:
        raise ValueError(f"{name} must define fixed values for all network components")

    return hypercube


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
        source: cast(int, source_value.value)
        for source, source_value in hypercube.items()
        if source_value.is_fixed
    }
    return robdd.can_take_value(fixed_values, value)


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
