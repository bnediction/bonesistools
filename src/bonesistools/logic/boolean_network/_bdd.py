#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Mapping, Tuple

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import _BDDManager, _BDDNode, _BDDTransitionRelationNode


def _bdd_equivalence(left: _BDDNode, right: _BDDNode) -> _BDDNode:
    """Return a BDD encoding logical equivalence."""

    return (left & right) | (~left & ~right)


def _bdd_exclusive_or(left: _BDDNode, right: _BDDNode) -> _BDDNode:
    """Return a BDD encoding exclusive disjunction."""

    return (left & ~right) | (~left & right)


def _bdd_forward_quantification(
    bdd: _BDDManager,
    clusters: Tuple[_BDDTransitionRelationNode, ...],
    current_vars: Tuple[str, ...],
) -> Tuple[Tuple[str, ...], ...]:
    """Assign each current variable to its last dependent partition."""

    current_variables = set(current_vars)
    last_cluster = {variable: 0 for variable in current_vars}
    for index, cluster in enumerate(clusters):
        for variable in bdd.support(cluster) & current_variables:
            last_cluster[variable] = index

    return tuple(
        tuple(variable for variable in current_vars if last_cluster[variable] == index)
        for index in range(len(clusters))
    )


def _bdd_rule(
    network: "BooleanNetwork",
    rule: Any,
    *,
    bdd: _BDDManager,
    variables: Mapping[str, str],
) -> _BDDNode:
    """Convert a Boolean rule into a BDD over current-state variables."""

    if rule is network.ba.TRUE:
        return bdd.true

    if rule is network.ba.FALSE:
        return bdd.false

    if isinstance(rule, network.ba.Symbol):
        return bdd.var(variables[str(rule)])

    if isinstance(rule, network.ba.NOT):
        return ~_bdd_rule(network, rule.args[0], bdd=bdd, variables=variables)

    if isinstance(rule, network.ba.AND):
        value = bdd.true
        for child in rule.args:
            value &= _bdd_rule(network, child, bdd=bdd, variables=variables)
        return value

    if isinstance(rule, network.ba.OR):
        value = bdd.false
        for child in rule.args:
            value |= _bdd_rule(network, child, bdd=bdd, variables=variables)
        return value

    raise TypeError(f"unsupported Boolean expression type: {type(rule)}")
