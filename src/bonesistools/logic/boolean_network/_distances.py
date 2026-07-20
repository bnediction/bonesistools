#!/usr/bin/env python

from __future__ import annotations

from typing import (
    Any,
    Collection,
    Dict,
    FrozenSet,
    List,
    Sequence,
    Tuple,
    cast,
)

import numpy as np
from boolean import BooleanAlgebra, Expression
from boolean.boolean import _FALSE, _TRUE

from ..boolean_algebra import ROBDD, equivalence
from ._network import BooleanNetwork
from ._typing import (
    BooleanNetworkLike,
    BooleanNetworkMetric,
    is_boolean_network_like,
)


def distance(
    bn1: BooleanNetworkLike,
    bn2: BooleanNetworkLike,
    metric: BooleanNetworkMetric = "hamming",
) -> float:
    """
    Compare two Boolean networks through their Boolean update functions.

    Networks must contain the same components. Logical representations are
    compared semantically, so syntactically different but equivalent update
    rules are treated as identical.

    Examples
    --------
    Only the update function of `A` differs between these networks. It differs
    on one of the four configurations of `A` and `B`:

    >>> from bonesistools.logic.bn import BooleanNetwork, distance
    >>> bn1 = BooleanNetwork({"A": "A & B", "B": "B"})
    >>> bn2 = BooleanNetwork({"A": "A", "B": "B"})
    >>> distance(bn1, bn2, metric="hamming")
    0.125
    >>> distance(bn1, bn2, metric="equivalence")
    0.5

    Parameters
    ----------
    bn1, bn2: BooleanNetworkLike
        Boolean networks or rule mappings to compare. Component order may
        differ, but component sets must be identical.
    metric: {"equivalence", "hamming"} (default: "hamming")
        Semantic comparison metric.

        - `equivalence` measures the proportion of components whose update
          functions are logically non-equivalent. Each component contributes
          either zero or one before averaging.
        - `hamming` measures the mean proportion of Boolean configurations on
          which corresponding update functions disagree. Each component
          receives equal weight.

    Returns
    -------
    float
        Semantic distance in the interval `[0, 1]`.

    Notes
    -----
    For example, `B & C` and `B` are non-equivalent, but differ only when
    `B = 1` and `C = 0`. Their local equivalence distance is therefore `1`,
    whereas their local Hamming distance is `1 / 4`.

    This function does not compare influence-graph structure, attractors,
    reachability, transition relations or trap spaces. Compare structure
    explicitly with:

    >>> import bonesistools as bt
    >>> bt.logic.ig.distance(
    ...     bn1.to_influence_graph(),
    ...     bn2.to_influence_graph(),
    ...     universe="union",
    ... )
    0.3333333333333333

    Raises
    ------
    TypeError
        If either argument is not a BooleanNetworkLike object.
    ValueError
        If the networks have different component sets or `metric` is
        unsupported.
    """

    network1, network2 = _coerce_networks(bn1, bn2)
    components = _common_components(network1, network2)
    agreement_count, disagreement_count = _network_comparison_counts(
        network1,
        network2,
        components,
        metric=metric,
    )
    comparison_count = agreement_count + disagreement_count

    if comparison_count == 0:
        return 0.0

    return disagreement_count / comparison_count


def similarity(
    bn1: BooleanNetworkLike,
    bn2: BooleanNetworkLike,
    metric: BooleanNetworkMetric = "hamming",
) -> float:
    """
    Compare two Boolean networks through their Boolean update functions.

    Networks must contain the same components. Logical representations are
    compared semantically, so syntactically different but equivalent update
    rules are treated as identical.

    Examples
    --------
    Only the update function of `A` differs between these networks. It agrees
    on three of the four configurations of `A` and `B`:

    >>> from bonesistools.logic.bn import BooleanNetwork, similarity
    >>> bn1 = BooleanNetwork({"A": "A & B", "B": "B"})
    >>> bn2 = BooleanNetwork({"A": "A", "B": "B"})
    >>> similarity(bn1, bn2, metric="hamming")
    0.875
    >>> similarity(bn1, bn2, metric="equivalence")
    0.5

    Parameters
    ----------
    bn1, bn2: BooleanNetworkLike
        Boolean networks or rule mappings to compare. Component order may
        differ, but component sets must be identical.
    metric: {"equivalence", "hamming"} (default: "hamming")
        Semantic comparison metric.

        - `equivalence` measures the proportion of components whose update
          functions are logically equivalent. Each component contributes
          either zero or one before averaging.
        - `hamming` measures the mean proportion of Boolean configurations on
          which corresponding update functions agree. Each component receives
          equal weight.

    Returns
    -------
    float
        Semantic similarity in the interval `[0, 1]`.

    Notes
    -----
    For example, `B & C` and `B` are non-equivalent, but agree on three of the
    four configurations of `B` and `C`. Their local equivalence similarity is
    therefore `0`, whereas their local Hamming similarity is `3 / 4`.

    This function does not compare influence-graph structure, attractors,
    reachability, transition relations or trap spaces. Compare structure
    explicitly with `bt.logic.ig.similarity()` on the corresponding influence
    graphs.

    Raises
    ------
    TypeError
        If either argument is not a BooleanNetworkLike object.
    ValueError
        If the networks have different component sets or `metric` is
        unsupported.
    """

    network1, network2 = _coerce_networks(bn1, bn2)
    components = _common_components(network1, network2)
    agreement_count, disagreement_count = _network_comparison_counts(
        network1,
        network2,
        components,
        metric=metric,
    )
    comparison_count = agreement_count + disagreement_count

    if comparison_count == 0:
        return 1.0

    return agreement_count / comparison_count


def _pairwise_distance_matrix(
    networks: Sequence[BooleanNetwork],
    *,
    metric: BooleanNetworkMetric,
) -> np.ndarray:
    """Return pairwise semantic distances between Boolean networks."""

    if metric not in ("equivalence", "hamming"):
        raise ValueError(f"unsupported Boolean-network comparison metric: {metric!r}")

    n_networks = len(networks)
    if n_networks == 0:
        return np.empty((0, 0), dtype=np.float64)

    components = tuple(sorted(networks[0]))
    component_set = set(components)
    if any(set(network) != component_set for network in networks[1:]):
        raise ValueError("Boolean networks must have the same component set")

    if n_networks == 1 or not components:
        return np.zeros((n_networks, n_networks), dtype=np.float64)

    rule_groups = _component_rule_groups(networks, components)
    if metric == "equivalence":
        return _equivalence_distance_matrix(
            networks,
            rule_groups,
        )

    return _hamming_distance_matrix(
        networks,
        rule_groups,
    )


def _component_rule_groups(
    networks: Sequence[BooleanNetwork],
    components: Sequence[str],
) -> List[Tuple[np.ndarray, Tuple[Expression, ...]]]:
    """Group structurally identical rules for every component."""

    groups = []
    for component in components:
        rule_indices: Dict[Expression, int] = {}
        rules = []
        indices = np.empty(len(networks), dtype=np.intp)

        for network_index, network in enumerate(networks):
            rule = network[component]
            if rule not in rule_indices:
                rule_indices[rule] = len(rules)
                rules.append(rule)
            indices[network_index] = rule_indices[rule]

        groups.append((indices, tuple(rules)))

    return groups


def _equivalence_distance_matrix(
    networks: Sequence[BooleanNetwork],
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
) -> np.ndarray:
    """Return pairwise proportions of non-equivalent update functions."""

    n_networks = len(networks)
    changed_counts = np.zeros((n_networks, n_networks), dtype=np.int64)
    comparison_cache: Dict[FrozenSet[Expression], bool] = {}
    ba = networks[0].ba

    for indices, rules in rule_groups:
        local_changes = np.zeros((len(rules), len(rules)), dtype=np.int64)

        for left_index, left in enumerate(rules):
            for right_index in range(left_index + 1, len(rules)):
                right = rules[right_index]
                key = frozenset((left, right))
                if key not in comparison_cache:
                    comparison_cache[key] = not equivalence(left, right, ba=ba)
                if comparison_cache[key]:
                    local_changes[left_index, right_index] = 1
                    local_changes[right_index, left_index] = 1

        changed_counts += local_changes[indices[:, None], indices]

    return changed_counts.astype(np.float64) / len(rule_groups)


def _hamming_distance_matrix(
    networks: Sequence[BooleanNetwork],
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
) -> np.ndarray:
    """Return mean pairwise functional disagreement proportions."""

    robdd_cache: Dict[Expression, ROBDD] = {}
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]] = {}
    ba = networks[0].ba
    max_variables = 0

    for _, rules in rule_groups:
        for left_index, left in enumerate(rules):
            for right in rules[left_index + 1 :]:
                _, n_variables = _robdd_disagreement_counts(
                    left,
                    right,
                    ba=ba,
                    robdd_cache=robdd_cache,
                    comparison_cache=comparison_cache,
                )
                max_variables = max(max_variables, n_variables)

    denominator = len(rule_groups) * (1 << max_variables)
    if denominator <= np.iinfo(np.int64).max:
        return _scaled_hamming_distance_matrix(
            len(networks),
            rule_groups,
            comparison_cache=comparison_cache,
            max_variables=max_variables,
            denominator=denominator,
        )

    return _unbounded_hamming_distance_matrix(
        len(networks),
        rule_groups,
        comparison_cache=comparison_cache,
    )


def _robdd_disagreement_counts(
    function1: Expression,
    function2: Expression,
    *,
    ba: BooleanAlgebra,
    robdd_cache: Dict[Expression, ROBDD],
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
) -> Tuple[int, int]:
    """Return disagreement and variable counts for two compiled functions."""

    key = frozenset((function1, function2))
    if key in comparison_cache:
        return comparison_cache[key]

    for function in (function1, function2):
        if function not in robdd_cache:
            robdd_cache[function] = ROBDD._from_expression(function, ba=ba)

    disagreement = robdd_cache[function1] ^ robdd_cache[function2]
    disagreement_count = disagreement.count(1)
    n_variables = len(disagreement.variables) if disagreement_count else 0
    comparison_cache[key] = disagreement_count, n_variables
    return comparison_cache[key]


def _scaled_hamming_distance_matrix(
    n_networks: int,
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
    max_variables: int,
    denominator: int,
) -> np.ndarray:
    """Accumulate Hamming numerators in a bounded integer matrix."""

    disagreement_counts = np.zeros((n_networks, n_networks), dtype=np.int64)

    for indices, rules in rule_groups:
        local_counts = np.zeros((len(rules), len(rules)), dtype=np.int64)
        for left_index, left in enumerate(rules):
            for right_index in range(left_index + 1, len(rules)):
                right = rules[right_index]
                count, n_variables = comparison_cache[frozenset((left, right))]
                scaled_count = count << (max_variables - n_variables)
                local_counts[left_index, right_index] = scaled_count
                local_counts[right_index, left_index] = scaled_count

        disagreement_counts += local_counts[indices[:, None], indices]

    return disagreement_counts.astype(np.float64) / denominator


def _unbounded_hamming_distance_matrix(
    n_networks: int,
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
) -> np.ndarray:
    """Accumulate exact Hamming counts without fixed-width integers."""

    distances = np.zeros((n_networks, n_networks), dtype=np.float64)

    for left_index in range(n_networks):
        for right_index in range(left_index + 1, n_networks):
            local_counts = []
            for indices, rules in rule_groups:
                left = rules[indices[left_index]]
                right = rules[indices[right_index]]
                if left == right:
                    local_counts.append((1, 0))
                    continue

                disagreement, n_variables = comparison_cache[frozenset((left, right))]
                assignment_count = 1 << n_variables
                local_counts.append((assignment_count - disagreement, disagreement))

            agreement, disagreement = _normalize_local_counts(local_counts)
            distance = disagreement / (agreement + disagreement)
            distances[left_index, right_index] = distance
            distances[right_index, left_index] = distance

    return distances


def _coerce_networks(
    bn1: BooleanNetworkLike,
    bn2: BooleanNetworkLike,
) -> Tuple[BooleanNetwork, BooleanNetwork]:
    """Return BooleanNetwork representations of two network-like objects."""

    for name, network in (("bn1", bn1), ("bn2", bn2)):
        if not is_boolean_network_like(network):
            raise TypeError(
                f"unsupported argument type for '{name}': expected "
                f"BooleanNetworkLike but received {type(network)}"
            )

    network1 = bn1 if isinstance(bn1, BooleanNetwork) else BooleanNetwork(bn1)
    network2 = (
        bn2 if isinstance(bn2, BooleanNetwork) else BooleanNetwork(bn2, ba=network1.ba)
    )

    return network1, network2


def _common_components(
    bn1: BooleanNetwork,
    bn2: BooleanNetwork,
) -> Tuple[str, ...]:
    """Validate and return a deterministic common component order."""

    components1 = set(bn1)
    components2 = set(bn2)

    if components1 != components2:
        raise ValueError("Boolean networks must have the same component set")

    return tuple(sorted(components1))


def _network_comparison_counts(
    bn1: BooleanNetwork,
    bn2: BooleanNetwork,
    components: Collection[str],
    *,
    metric: BooleanNetworkMetric,
) -> Tuple[int, int]:
    """Return exact aggregate agreement and disagreement counts."""

    if metric == "equivalence":
        agreement_count = sum(
            equivalence(bn1[component], bn2[component], ba=bn1.ba)
            for component in components
        )
        return agreement_count, len(components) - agreement_count

    if metric == "hamming":
        expression_cache: Dict[Expression, Expression] = {}
        local_counts = [
            _function_comparison_counts(
                bn1[component],
                bn2[component],
                ba=bn1.ba,
                other_ba=bn2.ba,
                expression_cache=expression_cache,
            )
            for component in components
        ]
        return _normalize_local_counts(local_counts)

    raise ValueError(f"unsupported Boolean-network comparison metric: {metric!r}")


def _function_comparison_counts(
    function1: Expression,
    function2: Expression,
    *,
    ba: BooleanAlgebra,
    other_ba: BooleanAlgebra,
    expression_cache: Dict[Expression, Expression],
) -> Tuple[int, int]:
    """Return exact agreement and disagreement counts for two functions."""

    if function1 == function2:
        return 1, 0

    if other_ba is not ba:
        function2 = _rebase_expression(
            function2,
            source_ba=other_ba,
            target_ba=ba,
            cache=expression_cache,
        )

    disagreement = cast(
        Expression,
        ba.OR(
            ba.AND(function1, ba.NOT(function2)),
            ba.AND(ba.NOT(function1), function2),
        ).simplify(),
    )

    if disagreement == ba.FALSE:
        return 1, 0
    if disagreement == ba.TRUE:
        return 0, 1

    variables = tuple(sorted(str(symbol) for symbol in disagreement.symbols))
    robdd = ROBDD._from_expression(
        disagreement,
        order=variables,
        ba=ba,
    )
    assignment_count = 1 << len(variables)
    disagreement_count = robdd.count(1)

    return assignment_count - disagreement_count, disagreement_count


def _rebase_expression(
    expression: Expression,
    *,
    source_ba: BooleanAlgebra,
    target_ba: BooleanAlgebra,
    cache: Dict[Expression, Expression],
) -> Expression:
    """Rebuild an expression with another compatible Boolean algebra."""

    if expression in cache:
        return cache[expression]

    expression_any: Any = expression

    if isinstance(expression_any, _FALSE):
        rebased = cast(Expression, target_ba.FALSE)
    elif isinstance(expression_any, _TRUE):
        rebased = cast(Expression, target_ba.TRUE)
    elif isinstance(expression_any, source_ba.Symbol):
        rebased = cast(Expression, target_ba.Symbol(str(expression_any.obj)))
    elif isinstance(expression_any, source_ba.NOT):
        rebased = cast(
            Expression,
            target_ba.NOT(
                _rebase_expression(
                    expression_any.args[0],
                    source_ba=source_ba,
                    target_ba=target_ba,
                    cache=cache,
                )
            ),
        )
    elif isinstance(expression_any, source_ba.AND):
        rebased = cast(
            Expression,
            target_ba.AND(
                *(
                    _rebase_expression(
                        child,
                        source_ba=source_ba,
                        target_ba=target_ba,
                        cache=cache,
                    )
                    for child in expression_any.args
                )
            ),
        )
    elif isinstance(expression_any, source_ba.OR):
        rebased = cast(
            Expression,
            target_ba.OR(
                *(
                    _rebase_expression(
                        child,
                        source_ba=source_ba,
                        target_ba=target_ba,
                        cache=cache,
                    )
                    for child in expression_any.args
                )
            ),
        )
    else:
        raise TypeError(f"unsupported Boolean expression type: {type(expression)}")

    cache[expression] = rebased
    return rebased


def _normalize_local_counts(
    local_counts: List[Tuple[int, int]],
) -> Tuple[int, int]:
    """Aggregate local proportions with equal component weights."""

    if not local_counts:
        return 0, 0

    common_count = max(
        agreement_count + disagreement_count
        for agreement_count, disagreement_count in local_counts
    )
    agreement_count = 0
    disagreement_count = 0

    for local_agreement, local_disagreement in local_counts:
        local_count = local_agreement + local_disagreement
        scale = common_count // local_count
        agreement_count += local_agreement * scale
        disagreement_count += local_disagreement * scale

    return agreement_count, disagreement_count
