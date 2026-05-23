#!/usr/bin/env python

"""
Typing aliases and runtime validators for Boolean algebra objects.
"""

from typing import (
    Any,
    Mapping,
    Union,
)

from boolean import Expression
from ._boolean import PartialBoolean

BooleanRule = Union[Expression, bool, int, str]

PartialBooleanLike = Union[PartialBoolean, bool, int, str]
ConfigurationLike = Mapping[str, bool]
HypercubeLike = Mapping[str, PartialBooleanLike]


def is_boolean_expression_available() -> bool:
    """
    Test whether the `boolean.Expression` type is available.

    Examples
    --------
    >>> is_boolean_expression_available()
    True

    Returns
    -------
    bool
        True if the `boolean.Expression` type can be used for runtime checks.
    """
    return Expression is not type(NotImplemented)


def is_boolean_expression_like(obj: Any) -> bool:
    """
    Test whether an object behaves as a Boolean expression.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> is_boolean_expression_like(ba.parse("A & B"))
    True
    >>> is_boolean_expression_like("A & B")
    False

    Parameters
    ----------
    obj: Any
        Object to test.

    Returns
    -------
    bool
        True if `obj` is a `boolean.Expression` instance.
    """
    return is_boolean_expression_available() and isinstance(obj, Expression)


def is_boolean_rule_like(obj: Any) -> bool:
    """
    Test whether an object can be interpreted as a Boolean rule.

    Boolean-rule-like objects include `boolean.Expression` instances, Boolean
    constants, numeric constants 0 and 1, and strings.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> is_boolean_rule_like(ba.parse("A & B"))
    True
    >>> is_boolean_rule_like("A & B")
    True
    >>> is_boolean_rule_like(2)
    False

    Parameters
    ----------
    obj: Any
        Object to test.

    Returns
    -------
    bool
        True if `obj` can be interpreted as a Boolean rule.
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

    Examples
    --------
    >>> is_partial_boolean_like(0)
    True
    >>> is_partial_boolean_like("*")
    True
    >>> is_partial_boolean_like(2)
    False

    Parameters
    ----------
    value: Any
        Object to test.

    Returns
    -------
    bool
        True if `value` can be converted to a PartialBoolean.
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

    Examples
    --------
    >>> is_hypercube_like({"A": 0, "B": "*"})
    True
    >>> is_hypercube_like({"A": 2})
    False
    >>> is_hypercube_like({1: 0})
    False

    Parameters
    ----------
    obj: Any
        Object to test.

    Returns
    -------
    bool
        True if `obj` is a mapping from string components to
        PartialBoolean-like values.
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
