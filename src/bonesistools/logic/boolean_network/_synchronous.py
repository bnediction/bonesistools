#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Tuple

from ._bdd import _bdd_equivalence, _bdd_forward_quantification, _bdd_rule

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import _BDDManager, _BDDTransitionRelationNode


def _bdd_synchronous_transition_partitions(
    network: "BooleanNetwork",
    *,
    bdd: _BDDManager,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
) -> Tuple[Tuple[_BDDTransitionRelationNode, ...], Tuple[Tuple[str, ...], ...]]:
    """Build conjunctive partitions for the synchronous transition relation."""

    component_vars = dict(zip(components, current_vars))
    clusters = tuple(
        _bdd_equivalence(
            bdd.var(next_var),
            _bdd_rule(
                network,
                network[component],
                bdd=bdd,
                variables=component_vars,
            ),
        )
        for component, next_var in zip(components, next_vars)
    )
    forward_quantification = _bdd_forward_quantification(
        bdd,
        clusters,
        current_vars,
    )

    return clusters, forward_quantification


def _synchronous_successor_configuration_bits(
    configuration: int,
    next_configuration: int,
) -> Tuple[int, ...]:
    """Return the deterministic synchronous successor when it changes."""

    if next_configuration == configuration:
        return ()

    return (next_configuration,)
