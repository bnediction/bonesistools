#!/usr/bin/env python

import math
from typing import Any, FrozenSet, Union

PartialBooleanValue = Union[bool, int, float, str]


class PartialBoolean:
    """
    Partial Boolean powerset abstraction over a Boolean value.

    A PartialBoolean is not a Kleene truth value and does not implement a
    three-valued logic. It is a powerset abstraction representing a non-empty
    subset of the Boolean domain B = {0, 1}:

        0   maps to {0}
        1   maps to {1}
        *   maps to {0, 1}

    Thus "*" denotes an unspecified, free, or unresolved Boolean value: both
    Boolean values remain admissible. This is the standard set-theoretic
    interpretation of Boolean hypercube coordinates.

    The induced information lattice is ordered by set inclusion:

        0 < *
        1 < *

    while 0 and 1 are incomparable.

    This is an information/generalisation lattice, not a truth lattice. In
    particular, it is distinct from the Kleene K3 truth lattice 0 < * < 1.

    Parameters
    ----------
    value: bool, number, str or PartialBoolean
        Partial Boolean value. Supported values are False, True, 0, 1, NaN and
        "*".

    Raises
    ------
    ValueError
        If the provided value is not supported.
    """

    __slots__ = ("_value",)

    def __init__(self, value: Any) -> None:

        self._value = self._coerce_value(value)

    def __repr__(self) -> str:

        return f"PartialBoolean({self._value!r})"

    def __str__(self) -> str:
        """
        Return the canonical string representation of the partial Boolean value.

        Returns
        -------
        str
            String representation of 0, 1 or "*".
        """

        return str(self._value)

    def __bool__(self) -> bool:

        if self.is_free:
            raise ValueError("cannot convert free PartialBoolean to bool")

        return bool(self._value)

    def __eq__(self, other: object) -> bool:
        """
        Test equality with another partial Boolean-like value.

        Parameters
        ----------
        other: object
            PartialBoolean-like value to compare against.

        Returns
        -------
        bool or NotImplemented
            True if both values are equal. Returns NotImplemented when `other`
            cannot be interpreted as a PartialBoolean value.
        """

        if not isinstance(other, PartialBoolean):

            try:
                other = PartialBoolean(other)

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
        """
        Test strict inclusion in another PartialBoolean.

        PartialBoolean values use the set-theoretic order:

            0 < *
            1 < *

        Examples
        --------
        >>> PartialBoolean(0) < PartialBoolean("*")
        True

        >>> PartialBoolean(1) < PartialBoolean("*")
        True

        >>> PartialBoolean(0) < PartialBoolean(1)
        False

        Parameters
        ----------
        other: object
            Partial Boolean value used for comparison.

        Returns
        -------
        bool or NotImplemented
            Whether the current admissible-value set is strictly included in
            `other`.
        """

        other = self._coerce_other(other)

        if not isinstance(other, PartialBoolean):
            return NotImplemented

        return self.as_set() < other.as_set()

    def __gt__(self, other: object) -> bool:
        """
        Test strict containment of another PartialBoolean.

        PartialBoolean values use the set-theoretic order:

            0 < *
            1 < *

        Examples
        --------
        >>> PartialBoolean("*") > PartialBoolean(0)
        True

        >>> PartialBoolean("*") > PartialBoolean(1)
        True

        >>> PartialBoolean(0) > PartialBoolean(1)
        False

        Parameters
        ----------
        other: object
            Partial Boolean value used for comparison.

        Returns
        -------
        bool or NotImplemented
            Whether the current admissible-value set strictly contains `other`.
        """

        other = self._coerce_other(other)

        if not isinstance(other, PartialBoolean):
            return NotImplemented

        return self.as_set() > other.as_set()

    def __le__(self, other: object) -> bool:
        """
        Test inclusion in another PartialBoolean.

        PartialBoolean values use the set-theoretic order:

            0 < *
            1 < *

        Examples
        --------
        >>> PartialBoolean(0) <= PartialBoolean("*")
        True

        >>> PartialBoolean("*") <= PartialBoolean("*")
        True

        >>> PartialBoolean(0) <= PartialBoolean(1)
        False

        Parameters
        ----------
        other: object
            Partial Boolean value used for comparison.

        Returns
        -------
        bool or NotImplemented
            Whether the current admissible-value set is included in `other`.
        """

        other = self._coerce_other(other)

        if not isinstance(other, PartialBoolean):
            return NotImplemented

        return self.as_set() <= other.as_set()

    def __ge__(self, other: object) -> bool:
        """
        Test containment of another PartialBoolean.

        PartialBoolean values use the set-theoretic order:

            0 < *
            1 < *

        Examples
        --------
        >>> PartialBoolean("*") >= PartialBoolean(0)
        True

        >>> PartialBoolean("*") >= PartialBoolean("*")
        True

        >>> PartialBoolean(0) >= PartialBoolean(1)
        False

        Parameters
        ----------
        other: object
            Partial Boolean value used for comparison.

        Returns
        -------
        bool or NotImplemented
            Whether the current admissible-value set contains `other`.
        """

        other = self._coerce_other(other)

        if not isinstance(other, PartialBoolean):
            return NotImplemented

        return self.as_set() >= other.as_set()

    def __contains__(self, other: object) -> bool:
        """
        Test whether another partial Boolean value is contained in this value.

        Returns
        -------
        bool
            Whether `other` is contained in the current PartialBoolean.
        """

        return self.contains(other)

    def to_kleene(self):
        """
        Convert to a KleeneValue preserving the external symbol.

        This conversion changes interpretation: `PartialBoolean("*")` means
        the admissible Boolean set {0, 1}, while `KleeneValue("*")` means an
        unknown truth value.

        Returns
        -------
        KleeneValue
            Kleene truth value with the same external symbol.
        """

        from ._kleene import KleeneValue

        return KleeneValue(self._value)

    @property
    def value(self) -> Union[int, str]:
        """
        Underlying PartialBoolean value.

        Returns
        -------
        int or str
            Canonical value, either 0, 1 or "*".
        """

        return self._value

    @property
    def is_fixed(self) -> bool:
        """
        Whether the PartialBoolean is fixed to a Boolean value.

        Returns
        -------
        bool
            True if the value is 0 or 1.
        """

        return self._value in [0, 1]

    @property
    def is_free(self) -> bool:
        """
        Whether the PartialBoolean is free.

        Returns
        -------
        bool
            True if the value is "*".
        """

        return self._value == "*"

    def as_set(self) -> FrozenSet[int]:
        """
        Return the admissible Boolean values represented by this object.

        Returns
        -------
        frozenset
            `frozenset({0})`, `frozenset({1})` or `frozenset({0, 1})`.
        """

        if self._value == "*":
            return frozenset({0, 1})

        if self._value == 0:
            return frozenset({0})

        return frozenset({1})

    def possibilities(self) -> FrozenSet[int]:
        """
        Return the admissible Boolean values represented by this object.

        Returns
        -------
        frozenset
            `frozenset({0})`, `frozenset({1})` or `frozenset({0, 1})`.
        """

        return self.as_set()

    def contains(self, other: object) -> bool:
        """
        Test whether the PartialBoolean contains another partial Boolean value.

        Containment follows set inclusion of admissible Boolean values:
        "*" contains 0, 1 and "*", while fixed values only contain themselves.

        Parameters
        ----------
        other: object
            PartialBoolean-like value to test.

        Returns
        -------
        bool
            Whether `other` is contained in the current PartialBoolean.

        Raises
        ------
        TypeError
            If `other` has an unsupported type.
        ValueError
            If `other` has an unsupported value.
        """

        other = other if isinstance(other, PartialBoolean) else PartialBoolean(other)

        return other.as_set() <= self.as_set()

    @staticmethod
    def _coerce_value(value: Any) -> Union[int, str]:
        """
        Convert supported values into canonical PartialBoolean values.

        Parameters
        ----------
        value: Any
            Value to convert.

        Returns
        -------
        int or str
            Canonical PartialBoolean value.

        Raises
        ------
        ValueError
            If the provided value is not supported.
        """

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
    def _coerce_other(other: object) -> Union["PartialBoolean", Any]:

        if isinstance(other, PartialBoolean):
            return other

        try:
            return PartialBoolean(other)

        except (TypeError, ValueError):
            return NotImplemented
