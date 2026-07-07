#!/usr/bin/env python

from __future__ import annotations

import math
from numbers import Number
from typing import Any, Union

KleeneValueLike = Union[bool, int, float, str, "KleeneValue"]


class KleeneValue:
    """
    Truth-lattice value in strong Kleene three-valued logic.

    A KleeneValue represents a truth value in K3 logic:

        0   false
        "*" unknown
        1   true

    The Kleene truth lattice is the chain:

        0 < * < 1

    Lattice operations are available as instance methods:

        x.meet(y) greatest lower bound
        x.join(y) least upper bound

    Because the Kleene truth lattice is a chain, meet and join coincide with
    min and max in the truth order.

    Kleene logical operations are also implemented through Python operators:

        ~x    Kleene negation
        x & y equivalent to x.meet(y)
        x | y equivalent to x.join(y)

    `KleeneValue("*")` denotes an unknown truth value, not the set {0, 1}.
    It is therefore distinct from `PartialBoolean("*")`, even though the same
    external symbol "*" is used.

    Parameters
    ----------
    value: bool, number, str or KleeneValue
        Kleene truth value. Supported values are `False`, `True`, `0`, `1`,
        `NaN` and `"*"`.

    Raises
    ------
    ValueError
        If the provided value is not supported.
    """

    __slots__ = ("_value",)

    def __init__(self, value: Any) -> None:

        self._value = self._coerce_value(value)

    def __repr__(self) -> str:

        return f"KleeneValue({self._value!r})"

    def __str__(self) -> str:

        return str(self._value)

    def __bool__(self) -> bool:

        raise TypeError(
            "KleeneValue cannot be converted to bool. Use &, | and ~ instead."
        )

    def __eq__(self, other: object) -> bool:

        if not isinstance(other, KleeneValue):
            try:
                other = KleeneValue(other)

            except (TypeError, ValueError):
                return NotImplemented

        return self._value == other._value

    def __ne__(self, other: object) -> bool:

        result = self.__eq__(other)

        if result is NotImplemented:
            return NotImplemented

        return not result

    def __hash__(self) -> int:

        return hash(self._value)

    def __lt__(self, other: object) -> bool:

        other = self._coerce_other(other)

        if not isinstance(other, KleeneValue):
            return NotImplemented

        return self.rank < other.rank

    def __le__(self, other: object) -> bool:

        other = self._coerce_other(other)

        if not isinstance(other, KleeneValue):
            return NotImplemented

        return self.rank <= other.rank

    def __gt__(self, other: object) -> bool:

        other = self._coerce_other(other)

        if not isinstance(other, KleeneValue):
            return NotImplemented

        return self.rank > other.rank

    def __ge__(self, other: object) -> bool:

        other = self._coerce_other(other)

        if not isinstance(other, KleeneValue):
            return NotImplemented

        return self.rank >= other.rank

    def __invert__(self) -> "KleeneValue":

        if self.value == 0:
            return KleeneValue(1)

        if self.value == 1:
            return KleeneValue(0)

        return KleeneValue("*")

    def __and__(self, other: object) -> "KleeneValue":

        return self.meet(other)

    def __or__(self, other: object) -> "KleeneValue":

        return self.join(other)

    def to_partial_boolean(self):
        """
        Convert to a PartialBoolean preserving the external symbol.

        This conversion changes interpretation: `KleeneValue("*")` means an
        unknown truth value, while `PartialBoolean("*")` means the admissible
        Boolean set {0, 1}.

        Returns
        -------
        PartialBoolean
            Partial Boolean value with the same external symbol.
        """

        from ._boolean import PartialBoolean

        return PartialBoolean(self._value)

    @property
    def value(self) -> Union[int, str]:
        """
        Underlying Kleene value.

        Returns
        -------
        int or str
            Canonical value, either 0, 1 or "*".
        """

        return self._value

    @property
    def rank(self) -> int:
        """
        Truth-order rank.

        Returns
        -------
        int
            Rank in the truth order 0 < * < 1.
        """

        return {0: 0, "*": 1, 1: 2}[self._value]

    @property
    def is_unknown(self) -> bool:
        """
        Whether the Kleene truth value is unknown.

        Returns
        -------
        bool
            True if the value is "*".
        """

        return self._value == "*"

    def meet(self, other: object) -> "KleeneValue":
        """
        Return the greatest lower bound with another Kleene truth value.

        In the Kleene truth chain 0 < * < 1, meet coincides with the minimum
        value in the truth order.

        Parameters
        ----------
        other: object
            Kleene-like value to combine with the current value.

        Returns
        -------
        KleeneValue
            Greatest lower bound of `self` and `other`.
        """

        other = KleeneValue(other)

        return self if self <= other else other

    def join(self, other: object) -> "KleeneValue":
        """
        Return the least upper bound with another Kleene truth value.

        In the Kleene truth chain 0 < * < 1, join coincides with the maximum
        value in the truth order.

        Parameters
        ----------
        other: object
            Kleene-like value to combine with the current value.

        Returns
        -------
        KleeneValue
            Least upper bound of `self` and `other`.
        """

        other = KleeneValue(other)

        return self if self >= other else other

    @staticmethod
    def _coerce_value(value: Any) -> Union[int, str]:

        if isinstance(value, KleeneValue):
            return value.value

        from ._boolean import PartialBoolean

        if isinstance(value, PartialBoolean):
            return value.value

        if isinstance(value, bool):
            return int(value)

        if isinstance(value, float) and math.isnan(value):
            return "*"

        if value in [0, 1]:
            return int(value)

        if value == "*":
            return "*"

        raise ValueError(
            "invalid argument value for 'value': "
            f"expected 0, 1, NaN, False, True or '*' but received {value!r}"
        )

    @staticmethod
    def _coerce_other(other: object) -> Union["KleeneValue", Any]:

        if isinstance(other, KleeneValue):
            return other

        try:
            return KleeneValue(other)

        except (TypeError, ValueError):
            return NotImplemented


def meet(x: object, y: object) -> KleeneValue:
    """
    Return the Kleene truth-order meet.

    This is a convenience wrapper around `KleeneValue(x).meet(y)`.
    Prefer the instance method when working directly with KleeneValue objects.

    Returns
    -------
    KleeneValue
        Greatest lower bound of `x` and `y`.
    """

    return KleeneValue(x).meet(y)


def join(x: object, y: object) -> KleeneValue:
    """
    Return the Kleene truth-order join.

    This is a convenience wrapper around `KleeneValue(x).join(y)`.
    Prefer the instance method when working directly with KleeneValue objects.

    Returns
    -------
    KleeneValue
        Least upper bound of `x` and `y`.
    """

    return KleeneValue(x).join(y)


def diff(v1: Any, v2: Any) -> int:
    """
    Return the discrete differential between two Kleene truth values.

    This differential is not defined on the set-theoretic order of
    `PartialBoolean`. It interprets input symbols as `KleeneValue` objects and
    uses the Kleene truth order:

        0 < * < 1

    The differential between two Kleene values is defined as:

        d(x1, x2) =  1  if x1 < x2
                     0  if x1 = x2
                    -1  if x1 > x2

    This is a directional heuristic. It should not be confused with the
    information/generalisation order of `PartialBoolean`, where 0 and 1 are
    incomparable and both are below "*".

    Examples
    --------
    >>> diff(0, 1)
    1

    >>> diff(1, 0)
    -1

    >>> diff(0, "*")
    1

    >>> diff("*", 1)
    1

    >>> diff("*", "*")
    0

    Parameters
    ----------
    v1, v2: bool, number, str, PartialBoolean or KleeneValue
        Values compared after conversion to KleeneValue.

    Returns
    -------
    int
        Differential from `v1` to `v2`:
            - 1 if `v1 < v2`,
            - 0 if `v1 == v2`,
            - -1 if `v1 > v2`.
    """

    kv1 = _as_kleene_value(v1)
    kv2 = _as_kleene_value(v2)

    if kv1 == kv2:
        return 0

    return 1 if kv1 < kv2 else -1


def _as_kleene_value(value: Any) -> KleeneValue:
    """
    Convert a value to a KleeneValue.

    Parameters
    ----------
    value: Any
        Value to convert. Supported values are `False`, `True`, `0`, `1`,
        `NaN`, `"*"`, PartialBoolean values and KleeneValue values.

    Returns
    -------
    KleeneValue
        Converted Kleene truth value.

    Raises
    ------
    TypeError
        If `value` has an unsupported type.
    ValueError
        If `value` is numeric but not 0, 1 or NaN.
    """

    try:
        return value if isinstance(value, KleeneValue) else KleeneValue(value)

    except ValueError:
        if isinstance(value, str):
            raise

        if isinstance(value, Number):
            raise ValueError(
                "invalid argument value for 'value': "
                f"expected 0, 1, False, True, NaN or '*' but received {value!r}"
            )

        raise TypeError(
            "unsupported argument type for 'value': "
            f"expected bool, number, str, PartialBoolean or KleeneValue "
            f"but received {type(value)}"
        )
