#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Set, Tuple

from ._bdd import _bdd_equivalence, _bdd_exclusive_or, _bdd_rule

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _bdd_asynchronous_transition_partitions(
    network: "BooleanNetwork",
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Tuple[Any, ...]:
    """Build disjunctive partitions of the asynchronous transition relation."""

    component_vars = dict(zip(components, current_vars))
    current_values = tuple(bdd.var(variable) for variable in current_vars)
    next_values = tuple(bdd.var(variable) for variable in next_vars)
    rules = tuple(
        _bdd_rule(network, network[component], bdd=bdd, variables=component_vars)
        for component in components
    )
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

    partitions = []
    for updated_index, (current_value, next_value, rule) in enumerate(
        zip(current_values, next_values, rules)
    ):
        transition = _bdd_exclusive_or(current_value, rule)
        transition &= _bdd_equivalence(next_value, rule)
        transition &= unchanged_prefixes[updated_index]
        transition &= unchanged_suffixes[updated_index + 1]
        partitions.append(transition)

    return tuple(partitions)


def _asynchronous_successor_configuration_bits(
    configuration: int,
    unstable_mask: int,
) -> Tuple[int, ...]:
    """Return successors obtained by updating one unstable component."""

    successors = []
    while unstable_mask:
        bit = unstable_mask & -unstable_mask
        successors.append(configuration ^ bit)
        unstable_mask ^= bit

    return tuple(successors)


def _tarjan_asynchronous_terminal_sccs(
    initial_configurations: Iterable[int],
    unstable_mask: Callable[[int], int],
) -> Tuple[Tuple[int, ...], ...]:
    """Find reachable terminal SCCs with an iterative on-the-fly Tarjan search."""

    visited: Set[int] = set()
    indices: Dict[int, int] = {}
    lowlinks: Dict[int, int] = {}
    component_stack: List[int] = []
    has_external_successor: Set[int] = set()
    terminal_components: List[Tuple[int, ...]] = []
    next_index = 0

    for initial_configuration in initial_configurations:
        if initial_configuration in visited:
            continue

        visited.add(initial_configuration)
        indices[initial_configuration] = next_index
        lowlinks[initial_configuration] = next_index
        next_index += 1
        component_stack.append(initial_configuration)
        dfs_configurations = [initial_configuration]
        remaining_masks = [unstable_mask(initial_configuration)]

        while dfs_configurations:
            configuration = dfs_configurations[-1]
            remaining_mask = remaining_masks[-1]
            if remaining_mask:
                bit = remaining_mask & -remaining_mask
                remaining_masks[-1] ^= bit
                successor = configuration ^ bit

                if successor not in visited:
                    visited.add(successor)
                    indices[successor] = next_index
                    lowlinks[successor] = next_index
                    next_index += 1
                    component_stack.append(successor)
                    dfs_configurations.append(successor)
                    remaining_masks.append(unstable_mask(successor))
                elif successor in indices:
                    lowlinks[configuration] = min(
                        lowlinks[configuration],
                        indices[successor],
                    )
                else:
                    has_external_successor.add(configuration)
                continue

            dfs_configurations.pop()
            remaining_masks.pop()

            if lowlinks[configuration] == indices[configuration]:
                component = []
                is_terminal = True
                while True:
                    member = component_stack.pop()
                    component.append(member)
                    is_terminal &= member not in has_external_successor
                    has_external_successor.discard(member)
                    del indices[member]
                    del lowlinks[member]
                    if member == configuration:
                        break

                if is_terminal:
                    terminal_components.append(tuple(component))

            if dfs_configurations:
                parent = dfs_configurations[-1]
                if configuration in indices:
                    lowlinks[parent] = min(
                        lowlinks[parent],
                        lowlinks[configuration],
                    )
                else:
                    has_external_successor.add(parent)

    return tuple(terminal_components)
