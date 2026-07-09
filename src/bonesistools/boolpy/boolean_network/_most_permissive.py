#!/usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import (
    TYPE_CHECKING,
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
from ..boolean_algebra import Hypercube, prime_implicants
from ..boolean_algebra._structure import Implicants
from ..boolean_algebra._typing import HypercubeLike

if TYPE_CHECKING:
    from ._network import BooleanNetwork

ImplicantsByValue = Dict[Literal[0, 1], Implicants]


def _reachability_with_most_permissive_hypercube_backend(
    network: "BooleanNetwork",
    initial_state: HypercubeLike,
    target_state: HypercubeLike,
) -> bool:
    """Test one-step most-permissive reachability with hypercube regions."""

    return _most_permissive_transition_space(network, initial_state).contains(
        target_state
    )


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
        """Test whether a complete configuration belongs to the region."""

        for component in components:
            if component in self.free_components:
                continue

            if target_state[component].value != initial_state[component].value:
                return False

        for component in self.irreversible_components:
            initial_value = cast(int, initial_state[component].value)
            if target_state[component].value != 1 - initial_value:
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
        """Test whether a complete configuration is represented."""

        target_hypercube = _validate_complete_hypercube(
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
) -> _MostPermissiveTransitionSpace:
    """Return the compact one-step most-permissive transition space."""

    components = tuple(network.keys())
    initial_hypercube = _validate_complete_hypercube(
        components,
        initial_state,
        name="initial_state",
    )

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
        implicants: Dict[str, ImplicantsByValue] = {}
        closed_hypercube = _smallest_closed_hypercube(
            network,
            initial_hypercube,
            components_to_relax=closure_components,
            implicants=implicants,
        )
        free_components = _free_components(components, closed_hypercube)
        irreversible_components = _irreversible_components(
            network,
            initial_hypercube,
            closed_hypercube,
            free_components=free_components,
            implicants=implicants,
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
    implicants: Dict[str, ImplicantsByValue],
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
            implicants=implicants,
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
    components_to_relax: Iterable[str],
    depth: Optional[int] = None,
    implicants: Optional[Dict[str, ImplicantsByValue]] = None,
) -> Hypercube:
    """Return the smallest closed hypercube obtained by relaxing components in K."""

    hypercube = Hypercube(initial_state)
    fixed = {
        component
        for component in components_to_relax
        if component in hypercube and hypercube[component].is_fixed
    }
    if implicants is None:
        implicants = {}

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
