#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Tuple

from ._bdd import _bdd_equivalence, _bdd_exclusive_or, _bdd_rule

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _bdd_asynchronous_transition_relation(
    network: "BooleanNetwork",
    *,
    bdd: Any,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Any:
    """Build the asynchronous transition relation as one BDD."""

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

    relation = bdd.false
    for updated_index, (current_value, next_value, rule) in enumerate(
        zip(current_values, next_values, rules)
    ):
        transition = _bdd_exclusive_or(current_value, rule)
        transition &= _bdd_equivalence(next_value, rule)
        transition &= unchanged_prefixes[updated_index]
        transition &= unchanged_suffixes[updated_index + 1]
        relation |= transition

    return relation


def _asynchronous_successor_state_bits(
    state: int,
    unstable_mask: int,
) -> Tuple[int, ...]:
    """Return successors obtained by updating one unstable component."""

    successors = []
    while unstable_mask:
        bit = unstable_mask & -unstable_mask
        successors.append(state ^ bit)
        unstable_mask ^= bit

    return tuple(successors)
