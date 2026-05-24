#!/usr/bin/env python

from typing import Any, Union

import math

PartialBooleanValue = Union[bool, int, float, str]


class PartialBoolean:
    """
    Partial Boolean value in {0, 1, *}.

    A PartialBoolean represents a Boolean abstraction where:
        - 0 denotes an inactive Boolean state,
        - 1 denotes an active Boolean state,
        - "*" denotes a free, unspecified or unresolved Boolean state.

    Two complementary semantics are supported:

        - an ensemble/set-theoretic interpretation through `contains`,
          where "*" represents the Boolean ensemble:

                * = {0, 1}

          while fixed values only contain themselves;

        - a biological ordering interpretation through comparison operators and
          differential calculus, following the convention:

                0 < * < 1

    These semantics are intentionally distinct. In particular, although "*"
    contains both 0 and 1 from an ensemble perspective, it is neither equal to
    0 nor to 1.

    Parameters
    ----------
    value: bool or int or str
        Partial Boolean value. Supported values are:
            - False, 0
            - True, 1
            - "*"

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

    def __lt__(self, other: "PartialBoolean") -> bool:
        """
        Return whether the current partial Boolean value is biologically lower
        than another one.

        Partial Boolean values follow the biological ordering convention:

            0 < * < 1

        Examples
        --------
        >>> PartialBoolean(0) < PartialBoolean("*")
        True

        >>> PartialBoolean("*") < PartialBoolean(1)
        True

        >>> PartialBoolean(1) < PartialBoolean(0)
        False

        Parameters
        ----------
        other: PartialBoolean
            Partial Boolean value used for comparison.

        Returns
        -------
        bool
            Whether the current value is biologically lower than `other`.
        """

        order = {
            PartialBoolean(0): 0,
            PartialBoolean(float("nan")): 1,
            PartialBoolean(1): 2,
        }

        return order[self] < order[other]

    def __gt__(self, other: "PartialBoolean") -> bool:
        """
        Return whether the current partial Boolean value is biologically greater
        than another one.

        Partial Boolean values follow the biological ordering convention:

            0 < * < 1

        Examples
        --------
        >>> PartialBoolean(1) > PartialBoolean("*")
        True

        >>> PartialBoolean("*") > PartialBoolean(0)
        True

        >>> PartialBoolean(0) > PartialBoolean(1)
        False

        Parameters
        ----------
        other: PartialBoolean
            Partial Boolean value used for comparison.

        Returns
        -------
        bool
            Whether the current value is biologically greater than `other`.
        """

        return other < self

    def __le__(self, other: "PartialBoolean") -> bool:
        """
        Return whether the current partial Boolean value is biologically lower
        than or equal to another one.

        Partial Boolean values follow the biological ordering convention:

            0 < * < 1

        Examples
        --------
        >>> PartialBoolean(0) <= PartialBoolean("*")
        True

        >>> PartialBoolean("*") <= PartialBoolean("*")
        True

        >>> PartialBoolean(1) <= PartialBoolean(0)
        False

        Parameters
        ----------
        other: PartialBoolean
            Partial Boolean value used for comparison.

        Returns
        -------
        bool
            Whether the current value is biologically lower than or equal to
            `other`.
        """

        return self == other or self < other

    def __ge__(self, other: "PartialBoolean") -> bool:
        """
        Return whether the current partial Boolean value is biologically greater
        than or equal to another one.

        Partial Boolean values follow the biological ordering convention:

            0 < * < 1

        Examples
        --------
        >>> PartialBoolean(1) >= PartialBoolean("*")
        True

        >>> PartialBoolean("*") >= PartialBoolean("*")
        True

        >>> PartialBoolean(0) >= PartialBoolean(1)
        False

        Parameters
        ----------
        other: PartialBoolean
            Partial Boolean value used for comparison.

        Returns
        -------
        bool
            Whether the current value is biologically greater than or equal to
            `other`.
        """

        return self == other or self > other

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

    def contains(self, other: object) -> bool:
        """
        Test whether the PartialBoolean contains another partial Boolean value.

        The containment relation follows the combinatorial interpretation:
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
        ValueError
            If `other` cannot be interpreted as a PartialBoolean value.
        """

        other = other if isinstance(other, PartialBoolean) else PartialBoolean(other)

        return self.is_free or self == other

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
