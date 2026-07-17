#!/usr/bin/env python

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    Iterable,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
)

from ._bdd import _bdd_equivalence, _bdd_exclusive_or, _bdd_rule

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import (
        _BDDConfigurationSetNode,
        _BDDManager,
        _BDDNode,
        _BDDTransitionRelationNode,
    )


class _AsynchronousTransition(NamedTuple):
    """Functional representation of one asynchronous update label."""

    current_variable: str
    current_value: _BDDNode
    rule: _BDDNode


def _bdd_asynchronous_transition_partitions(
    network: "BooleanNetwork",
    *,
    bdd: _BDDManager,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
) -> Tuple[_AsynchronousTransition, ...]:
    """Compile one local rule for each asynchronous transition label."""

    component_vars = dict(zip(components, current_vars))
    return tuple(
        _AsynchronousTransition(
            current_variable=variable,
            current_value=bdd.var(variable),
            rule=_bdd_rule(
                network,
                network[component],
                bdd=bdd,
                variables=component_vars,
            ),
        )
        for component, variable in zip(components, current_vars)
    )


def _bdd_asynchronous_relation_partitions(
    network: "BooleanNetwork",
    *,
    bdd: _BDDManager,
    components: Tuple[str, ...],
    current_vars: Tuple[str, ...],
    next_vars: Tuple[str, ...],
    functional_transitions: Optional[Tuple[_AsynchronousTransition, ...]] = None,
) -> Tuple[_BDDTransitionRelationNode, ...]:
    """Build full relation partitions for universal transition queries."""

    if functional_transitions is None:
        functional_transitions = _bdd_asynchronous_transition_partitions(
            network,
            bdd=bdd,
            components=components,
            current_vars=current_vars,
        )

    current_values = tuple(bdd.var(variable) for variable in current_vars)
    next_values = tuple(bdd.var(variable) for variable in next_vars)
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
    for updated_index, (transition, next_value) in enumerate(
        zip(functional_transitions, next_values)
    ):
        relation = _bdd_exclusive_or(
            transition.current_value,
            transition.rule,
        )
        relation &= _bdd_equivalence(next_value, transition.rule)
        relation &= unchanged_prefixes[updated_index]
        relation &= unchanged_suffixes[updated_index + 1]
        partitions.append(relation)

    return tuple(partitions)


def _bdd_asynchronous_transition_sources(
    states: _BDDConfigurationSetNode,
    transition: _AsynchronousTransition,
) -> _BDDConfigurationSetNode:
    """Return states enabling one functional asynchronous transition."""

    return states & _bdd_exclusive_or(
        transition.current_value,
        transition.rule,
    )


def _bdd_asynchronous_successors(
    bdd: _BDDManager,
    states: _BDDConfigurationSetNode,
    transition: _AsynchronousTransition,
) -> _BDDConfigurationSetNode:
    """Return successors through one functional asynchronous transition."""

    upward = states & ~transition.current_value & transition.rule
    downward = states & transition.current_value & ~transition.rule
    successors = bdd.false
    if upward != bdd.false:
        successors |= (
            bdd.exist(
                (transition.current_variable,),
                upward,
            )
            & transition.current_value
        )
    if downward != bdd.false:
        successors |= (
            bdd.exist(
                (transition.current_variable,),
                downward,
            )
            & ~transition.current_value
        )

    return successors


def _bdd_asynchronous_predecessors(
    bdd: _BDDManager,
    states: _BDDConfigurationSetNode,
    transition: _AsynchronousTransition,
) -> _BDDConfigurationSetNode:
    """Return predecessors through one functional asynchronous transition."""

    target_low = bdd.let({transition.current_variable: False}, states)
    target_high = bdd.let({transition.current_variable: True}, states)
    return (~transition.current_value & transition.rule & target_high) | (
        transition.current_value & ~transition.rule & target_low
    )


def _bdd_asynchronous_forward_chaining(
    bdd: _BDDManager,
    initial: _BDDConfigurationSetNode,
    transitions: Tuple[_AsynchronousTransition, ...],
    *,
    within: Optional[_BDDConfigurationSetNode] = None,
) -> _BDDConfigurationSetNode:
    """
    Return the asynchronous forward closure by transition-partition chaining.

    Each partition immediately sees configurations reached by earlier partitions
    in the same pass. For each partition, `processed[index]` records the
    configurations already propagated through it, so a newly reached
    configuration is processed only by partitions that have not seen it yet. If
    `within` is provided, successors outside that region are discarded.

    This implements transition-partition chaining as described in Section 2.5
    of [1], refined with one workset per partition. It is not the recursive MDD
    saturation algorithm that is the main subject of that article.

    References
    ----------
    [1] Ciardo, G., Marmorstein, R., & Siminiceanu, R. (2006). The saturation
    algorithm for symbolic state-space exploration. International Journal on
    Software Tools for Technology Transfer, 8(1), 4-25.
    """

    reachable = initial
    processed = [bdd.false for _ in transitions]

    while True:
        changed = False
        for index, transition in enumerate(transitions):
            fresh = reachable & ~processed[index]
            if fresh == bdd.false:
                continue

            processed[index] |= fresh
            successors = _bdd_asynchronous_successors(
                bdd,
                fresh,
                transition,
            )
            if within is not None:
                successors &= within
            new = successors & ~reachable
            if new != bdd.false:
                reachable |= new
                changed = True

        if not changed:
            return reachable


def _bdd_asynchronous_backward_chaining(
    bdd: _BDDManager,
    initial: _BDDConfigurationSetNode,
    transitions: Tuple[_AsynchronousTransition, ...],
    *,
    within: _BDDConfigurationSetNode,
) -> _BDDConfigurationSetNode:
    """
    Return the asynchronous backward closure by transition chaining.

    Each transition receives only target configurations that it has not
    processed before. Newly discovered predecessors are immediately available
    to subsequent transitions. The sweep direction alternates after each pass
    to reduce dependence on component order.

    This is the backward counterpart of the transition-partition chaining
    described in Section 2.5 of [1].

    References
    ----------
    [1] Ciardo, G., Marmorstein, R., & Siminiceanu, R. (2006). The saturation
    algorithm for symbolic state-space exploration. International Journal on
    Software Tools for Technology Transfer, 8(1), 4-25.
    """

    coreachable = initial
    processed = [bdd.false for _ in transitions]
    indexed_transitions = tuple(enumerate(transitions))
    reverse = False

    while True:
        changed = False
        ordered_transitions = (
            reversed(indexed_transitions) if reverse else iter(indexed_transitions)
        )
        for index, transition in ordered_transitions:
            fresh = coreachable & ~processed[index]
            if fresh == bdd.false:
                continue

            processed[index] |= fresh
            predecessors = _bdd_asynchronous_predecessors(
                bdd,
                fresh,
                transition,
            )
            new = predecessors & within & ~coreachable
            if new != bdd.false:
                coreachable |= new
                changed = True

        if not changed:
            return coreachable
        reverse = not reverse


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
