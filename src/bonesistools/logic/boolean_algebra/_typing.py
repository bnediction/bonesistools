#!/usr/bin/env python

"""
Typing aliases and runtime validators for Boolean algebra objects.
"""

from collections.abc import Mapping as _MappingABC
from typing import (
    Any as _Any,
)
from typing import (
    Mapping as _Mapping,
)
from typing import (
    Union as _Union,
)

from boolean import Expression
from typing_extensions import TypeGuard as _TypeGuard

from ._boolean import PartialBoolean
from ._kleene import KleeneValue, KleeneValueLike

BooleanRule = _Union[Expression, bool, int, str]

PartialBooleanLike = _Union[PartialBoolean, bool, int, float, str]
ConfigurationLike = _Mapping[str, _Union[bool, int]]
HypercubeLike = _Mapping[str, PartialBooleanLike]


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


def is_boolean_expression_like(obj: _Any) -> _TypeGuard[Expression]:
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


def is_boolean_rule_like(obj: _Any) -> _TypeGuard[BooleanRule]:
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


def is_partial_boolean_like(value: _Any) -> _TypeGuard[PartialBooleanLike]:
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


def is_kleene_value_like(value: _Any) -> _TypeGuard[KleeneValueLike]:
    """
    Test whether an object can be coerced into a KleeneValue.

    Examples
    --------
    >>> is_kleene_value_like(0)
    True
    >>> is_kleene_value_like("*")
    True
    >>> is_kleene_value_like(2)
    False

    Parameters
    ----------
    value: Any
        Object to test.

    Returns
    -------
    bool
        True if `value` can be converted to a KleeneValue.
    """

    try:
        KleeneValue(value)
        return True

    except (TypeError, ValueError):
        return False


def is_configuration_like(obj: _Any) -> _TypeGuard[ConfigurationLike]:
    """
    Test whether an object behaves as a concrete Boolean configuration.

    A configuration-like object maps string components to fixed Boolean values.
    Unlike hypercubes, configurations cannot contain the free value `"*"`.

    Examples
    --------
    >>> is_configuration_like({"A": 0, "B": True})
    True
    >>> is_configuration_like({"A": "*"})
    False
    >>> is_configuration_like(3)
    False

    Parameters
    ----------
    obj: Any
        Object to test.

    Returns
    -------
    bool
        True if `obj` is a mapping from string components to fixed Boolean
        values.
    """

    if not isinstance(obj, _MappingABC):
        return False

    try:
        return all(
            isinstance(component, str) and _is_configuration_value_like(value)
            for component, value in obj.items()
        )

    except Exception:
        return False


def is_hypercube_like(obj: _Any) -> _TypeGuard[HypercubeLike]:
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

    if not isinstance(obj, _MappingABC):
        return False

    try:
        return all(
            isinstance(component, str) and is_partial_boolean_like(value)
            for component, value in obj.items()
        )

    except Exception:
        return False


def _is_configuration_value_like(value: _Any) -> bool:
    try:
        return isinstance(value, (bool, int)) and value in [0, 1]

    except Exception:
        return False
