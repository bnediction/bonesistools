#!/usr/bin/env python

"""
Shared Booleanization helpers for multi-valued logical models.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, FrozenSet, Iterable, List, Mapping, Tuple

from ..boolean_network import BooleanNetwork


@dataclass(frozen=True)
class _LogicalModel:
    """Store desired threshold rules before Boolean regularization."""

    max_levels: Mapping[str, int]
    input_components: FrozenSet[str]
    threshold_rules: Mapping[str, Mapping[int, str]]
    component_names: Mapping[str, str]


def _booleanize_logical_model(model: _LogicalModel) -> BooleanNetwork:
    """Convert desired multi-valued threshold rules to a Boolean network."""

    rules = {}

    for source_component, max_level in model.max_levels.items():
        component = model.component_names[source_component]
        desired_rules = model.threshold_rules.get(source_component, {})

        if max_level < 1:
            raise ValueError(
                f"unsupported maximum level {max_level!r} for component "
                f"{source_component!r}: expected a positive integer"
            )

        for threshold in range(1, max_level + 1):
            boolean_component = _boolean_component(
                component,
                threshold=threshold,
                max_level=max_level,
            )

            if source_component in model.input_components:
                if max_level == 1 or threshold == 1:
                    rules[boolean_component] = boolean_component
                else:
                    rules[boolean_component] = _join_and(
                        _threshold_component(component, level)
                        for level in range(1, threshold + 1)
                    )
                continue

            try:
                desired_rule = desired_rules[threshold]
            except KeyError:
                raise ValueError(
                    f"missing desired threshold rule {threshold} for "
                    f"component {source_component!r}"
                )

            if max_level == 1:
                rules[boolean_component] = desired_rule
            else:
                rules[boolean_component] = _regularized_threshold_rule(
                    desired_rule,
                    component=component,
                    threshold=threshold,
                    max_level=max_level,
                )

    return BooleanNetwork(rules)


def _normalize_component_names(components: Iterable[str]) -> Dict[str, str]:
    """Normalize source component identifiers for Boolean expressions."""

    component_names = {
        source: re.sub(r"[^A-Za-z0-9_]", "_", source) for source in components
    }

    for source, component in component_names.items():
        if not component or component[0].isdigit():
            component_names[source] = f"_{component}"

    duplicates = [
        component
        for component in set(component_names.values())
        if list(component_names.values()).count(component) > 1
    ]

    if duplicates:
        raise ValueError(
            "Boolean component-name normalization creates duplicate names: "
            + ", ".join(sorted(duplicates))
        )

    return component_names


def _boolean_component(component: str, *, threshold: int, max_level: int) -> str:
    """Return the Boolean component representing one activity threshold."""

    if max_level == 1:
        return component

    return _threshold_component(component, threshold)


def _threshold_component(component: str, threshold: int) -> str:
    """Return the identifier of one threshold component."""

    return f"{component}_b{threshold}"


def _regularized_threshold_rule(
    desired_rule: str,
    *,
    component: str,
    threshold: int,
    max_level: int,
) -> str:
    """Regularize one desired threshold rule to preserve admissible updates."""

    terms = []

    if desired_rule != "0":
        desired_terms = [
            _threshold_component(component, level) for level in range(1, threshold)
        ]

        if desired_rule != "1":
            desired_terms.append(f"({desired_rule})")

        terms.append(_join_and(desired_terms))

    if threshold < max_level:
        terms.append(
            _join_and(
                _threshold_component(component, level)
                for level in range(1, threshold + 2)
            )
        )

    return _join_or(terms)


def _level_range_rule(
    component: str,
    *,
    minimum: int,
    maximum: int,
    max_level: int,
) -> str:
    """Encode a multi-valued level interval with Boolean threshold variables."""

    minimum = max(minimum, 0)
    maximum = min(maximum, max_level)

    if minimum > maximum:
        return "0"

    if minimum == 0 and maximum == max_level:
        return "1"

    if max_level == 1:
        if minimum == maximum == 0:
            return f"~{component}"
        if minimum == maximum == 1:
            return component

    level_conditions = []

    for level in range(minimum, maximum + 1):
        if level == 0:
            level_conditions.append(f"~{_threshold_component(component, 1)}")
            continue

        terms = [
            _threshold_component(component, threshold)
            for threshold in range(1, level + 1)
        ]

        if level < max_level:
            terms.append(f"~{_threshold_component(component, level + 1)}")

        level_conditions.append(_join_and(terms))

    rule = _join_or(level_conditions)
    return rule if len(level_conditions) == 1 else f"({rule})"


def _encode_state(
    state: Mapping[str, int],
    *,
    max_levels: Mapping[str, int],
    component_names: Mapping[str, str],
) -> Dict[str, int]:
    """Encode a multi-valued state as Boolean threshold assignments."""

    encoded = {}

    for source_component, value in state.items():
        if source_component not in max_levels:
            raise ValueError(f"unknown state component {source_component!r}")

        max_level = max_levels[source_component]

        if value < 0 or value > max_level:
            raise ValueError(
                f"unsupported state value {value!r} for {source_component!r}: "
                f"expected a value between 0 and {max_level}"
            )

        component = component_names[source_component]

        if max_level == 1:
            encoded[component] = value
            continue

        for threshold in range(1, max_level + 1):
            encoded[_threshold_component(component, threshold)] = int(
                value >= threshold
            )

    return encoded


def _syntactic_influences(
    boolean_network: BooleanNetwork,
) -> List[Tuple[str, str, int]]:
    """Extract signed rule literals without simplifying Boolean expressions."""

    influences = set()

    def collect(expression: Any, target: str, sign: int = 1) -> None:
        if isinstance(expression, boolean_network.ba.Symbol):
            influences.add((str(expression.obj), target, sign))
            return

        if isinstance(expression, boolean_network.ba.NOT):
            collect(expression.args[0], target, -sign)
            return

        for argument in expression.args:
            collect(argument, target, sign)

    for target, rule in boolean_network.items():
        collect(rule, target)

    return sorted(influences)


def _negate_rule(rule: str) -> str:
    """Negate a Boolean rule while simplifying constants and literals."""

    if rule == "0":
        return "1"

    if rule == "1":
        return "0"

    if rule.startswith("~") and " " not in rule:
        return rule[1:]

    return f"~({rule})"


def _join_or(expressions: Iterable[str]) -> str:
    """Join Boolean expressions with disjunction."""

    expressions = [expression for expression in expressions if expression]

    if not expressions:
        return "0"

    if len(expressions) == 1:
        return expressions[0]

    return " | ".join(f"({expression})" for expression in expressions)


def _join_and(expressions: Iterable[str]) -> str:
    """Join Boolean expressions with conjunction."""

    expressions = [expression for expression in expressions if expression != "1"]

    if "0" in expressions:
        return "0"

    if not expressions:
        return "1"

    return " & ".join(
        f"({expression})" if " | " in expression else expression
        for expression in expressions
    )
