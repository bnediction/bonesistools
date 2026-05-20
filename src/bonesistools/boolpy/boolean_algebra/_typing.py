#!/usr/bin/env python

from typing import Any, Mapping, Union

from boolean import Expression
from ._boolean import PartialBoolean

BooleanRule = Union[Expression, bool, int, str]

PartialBooleanLike = Union[PartialBoolean, bool, int, str]
HypercubeLike = Mapping[str, PartialBooleanLike]


def is_boolean_expression_available() -> bool:
    """
    Whether the `boolean.Expression` type is available.
    """
    return Expression is not type(NotImplemented)


def is_boolean_expression_like(obj: Any) -> bool:
    """
    Test whether an object behaves as a Boolean expression.
    """
    return is_boolean_expression_available() and isinstance(obj, Expression)


def is_boolean_rule_like(obj: Any) -> bool:
    """
    Test whether an object can be interpreted as a Boolean rule.
    """
    return (
        is_boolean_expression_like(obj)
        or isinstance(obj, bool)
        or obj in [0, 1]
        or isinstance(obj, str)
    )


def is_partial_boolean_like(value: Any) -> bool:
    """
    Test whether an object can be coerced into a PartialBoolean.
    """

    try:
        PartialBoolean(value)
        return True

    except (TypeError, ValueError):
        return False


def is_hypercube_like(obj: Any) -> bool:
    """
    Test whether an object behaves as a hypercube mapping.

    A hypercube-like object maps string components to PartialBoolean-like
    values.
    """

    if not isinstance(obj, Mapping):
        return False

    try:
        return all(
            isinstance(component, str) and is_partial_boolean_like(value)
            for component, value in obj.items()
        )

    except Exception:
        return False
