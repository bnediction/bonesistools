#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Tuple

from ._bdd import _bdd_equivalence, _bdd_forward_quantification, _bdd_rule

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import _BDDManager, _BDDTransitionRelationNode


def _bdd_general_transition_partitions(
    network: "BooleanNetwork",
    *,
    bdd: _BDDManager,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Tuple[Tuple[_BDDTransitionRelationNode, ...], Tuple[Tuple[str, ...], ...]]:
    """Build conjunctive partitions for general reachability analyses."""

    component_vars = dict(zip(components, current_vars))
    rules = tuple(
        _bdd_rule(network, network[component], bdd=bdd, variables=component_vars)
        for component in components
    )
    clusters = []
    for current_var, next_var, rule in zip(current_vars, next_vars, rules):
        current_value = bdd.var(current_var)
        next_value = bdd.var(next_var)
        # Identity transitions do not alter reachability or terminal SCCs and
        # avoid constructing one large disjunction of component updates.
        unchanged = _bdd_equivalence(next_value, current_value)
        updated = _bdd_equivalence(next_value, rule)
        clusters.append(unchanged | updated)

    forward_quantification = _bdd_forward_quantification(
        bdd,
        tuple(clusters),
        current_vars,
    )

    return tuple(clusters), forward_quantification


def _general_successor_configuration_bits(
    configuration: int,
    unstable_mask: int,
) -> Tuple[int, ...]:
    """Return successors from non-empty subsets of unstable components."""

    successors = []
    updated_mask = unstable_mask
    while updated_mask:
        successors.append(configuration ^ updated_mask)
        updated_mask = (updated_mask - 1) & unstable_mask

    return tuple(successors)
