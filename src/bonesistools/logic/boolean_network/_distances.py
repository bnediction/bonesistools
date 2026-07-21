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
from scipy.sparse import block_diag, csr_matrix

from ..boolean_algebra import ROBDD, equivalence
from ._network import BooleanNetwork
from ._typing import (
    BooleanNetworkLike,
    BooleanNetworkMetric,
    is_boolean_network_like,
)

_ROBDDNodeTable = Dict[int, Tuple[str, int, int]]
_ROBDDStructure = Tuple[int, _ROBDDNodeTable]


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
    varying_rule_groups = tuple(
        rule_group for rule_group in rule_groups if len(rule_group[1]) > 1
    )
    if not varying_rule_groups:
        return np.zeros((n_networks, n_networks), dtype=np.float64)

    if metric == "equivalence":
        return _equivalence_distance_matrix(
            networks,
            varying_rule_groups,
            n_components=len(components),
        )

    return _hamming_distance_matrix(
        networks,
        varying_rule_groups,
        n_components=len(components),
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
    *,
    n_components: int,
) -> np.ndarray:
    """Return pairwise proportions of non-equivalent update functions."""

    n_networks = len(networks)
    rule_classes = tuple(
        rule_class
        for rule_class in _equivalence_rule_classes(
            rule_groups,
            ba=networks[0].ba,
        )
        if rule_class[1] > 1
    )
    if not rule_classes:
        return np.zeros((n_networks, n_networks), dtype=np.float64)

    if _use_sparse_equivalence_matrix(n_networks, rule_classes):
        changed_counts = _sparse_equivalence_counts(n_networks, rule_classes)
    else:
        changed_counts = _dense_equivalence_counts(n_networks, rule_classes)

    return changed_counts.astype(np.float64) / n_components


def _equivalence_rule_classes(
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    ba: BooleanAlgebra,
) -> List[Tuple[np.ndarray, int]]:
    """Group semantically equivalent rules for every varying component."""

    comparison_cache: Dict[FrozenSet[Expression], bool] = {}
    class_groups = []

    for indices, rules in rule_groups:
        representatives = []
        rule_classes = np.empty(len(rules), dtype=np.intp)

        for rule_index, rule in enumerate(rules):
            for class_index, representative in enumerate(representatives):
                key = frozenset((rule, representative))
                if key not in comparison_cache:
                    comparison_cache[key] = equivalence(rule, representative, ba=ba)
                if comparison_cache[key]:
                    rule_classes[rule_index] = class_index
                    break
            else:
                rule_classes[rule_index] = len(representatives)
                representatives.append(rule)

        class_groups.append((rule_classes[indices], len(representatives)))

    return class_groups


def _use_sparse_equivalence_matrix(
    n_networks: int,
    rule_classes: Sequence[Tuple[np.ndarray, int]],
) -> bool:
    """Select sparse assembly from its estimated work relative to dense."""

    n_groups = len(rule_classes)
    agreement_entries = 0
    for indices, n_classes in rule_classes:
        frequencies = np.bincount(indices, minlength=n_classes)
        agreement_entries += sum(int(frequency) ** 2 for frequency in frequencies)

    dense_work = n_groups * n_networks * n_networks
    estimated_sparse_ratio = (
        1 / n_networks + 1 / n_groups + agreement_entries / dense_work
    )
    return estimated_sparse_ratio < 0.13


def _dense_equivalence_counts(
    n_networks: int,
    rule_classes: Sequence[Tuple[np.ndarray, int]],
) -> np.ndarray:
    """Accumulate non-equivalent rule counts in a dense matrix."""

    changed_counts = np.zeros((n_networks, n_networks), dtype=np.int64)
    for indices, _ in rule_classes:
        changed_counts += indices[:, None] != indices
    return changed_counts


def _sparse_equivalence_counts(
    n_networks: int,
    rule_classes: Sequence[Tuple[np.ndarray, int]],
) -> np.ndarray:
    """Accumulate non-equivalent rule counts through sparse agreements."""

    encoded = _rule_class_indicator_matrix(n_networks, rule_classes)
    agreement_counts = (encoded @ encoded.T).toarray()
    return len(rule_classes) - agreement_counts


def _hamming_distance_matrix(
    networks: Sequence[BooleanNetwork],
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    n_components: int,
) -> np.ndarray:
    """Return mean pairwise functional disagreement proportions."""

    robdd_cache: Dict[Expression, ROBDD] = {}
    robdd_structure_cache: Dict[ROBDD, _ROBDDStructure] = {}
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
                    robdd_structure_cache=robdd_structure_cache,
                    comparison_cache=comparison_cache,
                )
                max_variables = max(max_variables, n_variables)

    denominator = n_components * (1 << max_variables)
    if denominator <= np.iinfo(np.int64).max:
        local_counts, weight_entries, expanded_entries = _scaled_hamming_local_counts(
            rule_groups,
            comparison_cache=comparison_cache,
            max_variables=max_variables,
        )
        if _use_sparse_hamming_matrix(
            len(networks),
            len(rule_groups),
            weight_entries=weight_entries,
            expanded_entries=expanded_entries,
        ):
            return _sparse_scaled_hamming_distance_matrix(
                len(networks),
                rule_groups,
                local_counts,
                denominator=denominator,
            )

        return _dense_scaled_hamming_distance_matrix(
            len(networks),
            rule_groups,
            local_counts,
            denominator=denominator,
        )

    return _unbounded_hamming_distance_matrix(
        len(networks),
        rule_groups,
        comparison_cache=comparison_cache,
        max_variables=max_variables,
        denominator=denominator,
    )


def _robdd_disagreement_counts(
    function1: Expression,
    function2: Expression,
    *,
    ba: BooleanAlgebra,
    robdd_cache: Dict[Expression, ROBDD],
    robdd_structure_cache: Dict[ROBDD, _ROBDDStructure],
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
) -> Tuple[int, int]:
    """Return disagreement and variable counts for two compiled functions."""

    key = frozenset((function1, function2))
    if key in comparison_cache:
        return comparison_cache[key]

    for function in (function1, function2):
        if function not in robdd_cache:
            robdd_cache[function] = ROBDD._from_expression(function, ba=ba)

    robdd1 = robdd_cache[function1]
    robdd2 = robdd_cache[function2]
    disagreement_count = _count_robdd_disagreements(
        robdd1,
        robdd2,
        structure_cache=robdd_structure_cache,
    )
    n_variables = (
        len(set(robdd1.variables) | set(robdd2.variables)) if disagreement_count else 0
    )
    comparison_cache[key] = disagreement_count, n_variables
    return comparison_cache[key]


def _count_robdd_disagreements(
    robdd1: ROBDD,
    robdd2: ROBDD,
    *,
    structure_cache: Dict[ROBDD, _ROBDDStructure],
) -> int:
    """Count assignments on which two compiled functions disagree."""

    variables = tuple(sorted(set(robdd1.variables) | set(robdd2.variables)))
    variable_indices = {variable: index for index, variable in enumerate(variables)}
    root1, nodes1 = _robdd_structure(robdd1, cache=structure_cache)
    root2, nodes2 = _robdd_structure(robdd2, cache=structure_cache)
    n_variables = len(variables)
    count_cache: Dict[Tuple[int, int, int], int] = {}

    def count(node1: int, node2: int, next_variable: int) -> int:
        key = node1, node2, next_variable
        if key in count_cache:
            return count_cache[key]

        if node1 <= 1 and node2 <= 1:
            result = 1 << (n_variables - next_variable) if node1 != node2 else 0
            count_cache[key] = result
            return result

        variable1 = variable_indices[nodes1[node1][0]] if node1 > 1 else n_variables
        variable2 = variable_indices[nodes2[node2][0]] if node2 > 1 else n_variables
        variable = min(variable1, variable2)

        if variable1 == variable:
            _, low1, high1 = nodes1[node1]
        else:
            low1 = high1 = node1

        if variable2 == variable:
            _, low2, high2 = nodes2[node2]
        else:
            low2 = high2 = node2

        result = (1 << (variable - next_variable)) * (
            count(low1, low2, variable + 1) + count(high1, high2, variable + 1)
        )
        count_cache[key] = result
        return result

    return count(root1, root2, 0)


def _robdd_structure(
    robdd: ROBDD,
    *,
    cache: Dict[ROBDD, _ROBDDStructure],
) -> _ROBDDStructure:
    """Return cached decision nodes indexed by their identifiers."""

    if robdd not in cache:
        cache[robdd] = (
            robdd._root_id,
            {
                node: (variable, low, high)
                for node, variable, low, high in robdd._iter_nodes()
            },
        )
    return cache[robdd]


def _scaled_hamming_local_counts(
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
    max_variables: int,
) -> Tuple[List[np.ndarray], int, int]:
    """Build scaled local counts and sparse-work statistics."""

    local_matrices = []
    weight_entries = 0
    expanded_entries = 0

    for indices, rules in rule_groups:
        frequencies = np.bincount(indices, minlength=len(rules))
        local_counts = np.zeros((len(rules), len(rules)), dtype=np.int64)
        for left_index, left in enumerate(rules):
            for right_index in range(left_index + 1, len(rules)):
                right = rules[right_index]
                count, n_variables = comparison_cache[frozenset((left, right))]
                scaled_count = count << (max_variables - n_variables)
                local_counts[left_index, right_index] = scaled_count
                local_counts[right_index, left_index] = scaled_count
                if scaled_count:
                    weight_entries += 2
                    expanded_entries += int(
                        frequencies[left_index] + frequencies[right_index]
                    )

        local_matrices.append(local_counts)

    return local_matrices, weight_entries, expanded_entries


def _use_sparse_hamming_matrix(
    n_networks: int,
    n_groups: int,
    *,
    weight_entries: int,
    expanded_entries: int,
) -> bool:
    """Select sparse assembly from its estimated work relative to dense."""

    dense_work = n_groups * n_networks * n_networks
    estimated_sparse_ratio = (
        60 / n_networks
        + 1 / n_groups
        + weight_entries / dense_work
        + expanded_entries / dense_work
    )
    return estimated_sparse_ratio < 0.25


def _dense_scaled_hamming_distance_matrix(
    n_networks: int,
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    local_counts: Sequence[np.ndarray],
    *,
    denominator: int,
) -> np.ndarray:
    """Accumulate Hamming numerators in a bounded integer matrix."""

    disagreement_counts = np.zeros((n_networks, n_networks), dtype=np.int64)

    for (indices, _), component_counts in zip(rule_groups, local_counts):
        disagreement_counts += component_counts[indices[:, None], indices]

    return disagreement_counts.astype(np.float64) / denominator


def _sparse_scaled_hamming_distance_matrix(
    n_networks: int,
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    local_counts: Sequence[np.ndarray],
    *,
    denominator: int,
) -> np.ndarray:
    """Accumulate Hamming numerators through sparse matrix products."""

    rule_classes = tuple((indices, len(rules)) for indices, rules in rule_groups)
    encoded = _rule_class_indicator_matrix(n_networks, rule_classes)
    weights = block_diag(
        (csr_matrix(local_counts[0]), *local_counts[1:]),
        format="csr",
        dtype=np.int64,
    )
    disagreement_counts = (encoded @ weights @ encoded.T).toarray()
    return disagreement_counts.astype(np.float64) / denominator


def _rule_class_indicator_matrix(
    n_networks: int,
    rule_classes: Sequence[Tuple[np.ndarray, int]],
) -> csr_matrix:
    """Encode one selected rule class per component and network."""

    n_groups = len(rule_classes)
    class_counts = [n_classes for _, n_classes in rule_classes]
    offsets = np.cumsum(
        np.array([0, *class_counts[:-1]], dtype=np.intp),
    )
    columns = (
        np.column_stack([indices for indices, _ in rule_classes]) + offsets[None, :]
    ).ravel()
    rows = np.repeat(np.arange(n_networks), n_groups)
    return csr_matrix(
        (
            np.ones(rows.size, dtype=np.int64),
            (rows, columns),
        ),
        shape=(n_networks, sum(class_counts)),
        dtype=np.int64,
    )


def _unbounded_hamming_distance_matrix(
    n_networks: int,
    rule_groups: Sequence[Tuple[np.ndarray, Tuple[Expression, ...]]],
    *,
    comparison_cache: Dict[FrozenSet[Expression], Tuple[int, int]],
    max_variables: int,
    denominator: int,
) -> np.ndarray:
    """Accumulate exact Hamming counts without fixed-width integers."""

    distances = np.zeros((n_networks, n_networks), dtype=np.float64)

    for left_index in range(n_networks):
        for right_index in range(left_index + 1, n_networks):
            disagreement_count = 0
            for indices, rules in rule_groups:
                left = rules[indices[left_index]]
                right = rules[indices[right_index]]
                if left == right:
                    continue

                disagreement, n_variables = comparison_cache[frozenset((left, right))]
                disagreement_count += disagreement << (max_variables - n_variables)

            distance = disagreement_count / denominator
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
        robdd_cache: Dict[Expression, ROBDD] = {}
        robdd_structure_cache: Dict[ROBDD, _ROBDDStructure] = {}
        local_counts = [
            _function_comparison_counts(
                bn1[component],
                bn2[component],
                ba=bn1.ba,
                other_ba=bn2.ba,
                expression_cache=expression_cache,
                robdd_cache=robdd_cache,
                robdd_structure_cache=robdd_structure_cache,
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
    robdd_cache: Dict[Expression, ROBDD],
    robdd_structure_cache: Dict[ROBDD, _ROBDDStructure],
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

    for function in (function1, function2):
        if function not in robdd_cache:
            robdd_cache[function] = ROBDD._from_expression(function, ba=ba)

    robdd1 = robdd_cache[function1]
    robdd2 = robdd_cache[function2]
    disagreement_count = _count_robdd_disagreements(
        robdd1,
        robdd2,
        structure_cache=robdd_structure_cache,
    )
    if disagreement_count == 0:
        return 1, 0

    n_variables = len(set(robdd1.variables) | set(robdd2.variables))
    assignment_count = 1 << n_variables

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
