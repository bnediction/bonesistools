#!/usr/bin/env python

from typing import Any, Union

PartialBooleanValue = Union[bool, int, str]


class PartialBoolean:
    """
    Partial Boolean value in {0, 1, *}.

    A PartialBoolean represents a Boolean abstraction where:
        - 0 denotes an inactive Boolean state,
        - 1 denotes an active Boolean state,
        - "*" denotes a free, unspecified or unresolved Boolean state.

    No ordering is defined between PartialBoolean values. In particular,
    expressions such as `PartialBoolean(0) < PartialBoolean(1)` are not
    supported, because ordering depends on the intended semantics.

    The method `contains` provides the combinatorial interpretation used for
    hypercubes: "*" contains both Boolean values 0 and 1, while fixed values
    only contain themselves.

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

    def __init__(self, value: PartialBooleanValue) -> None:

        self._value = self._coerce_value(value)

    @property
    def value(self) -> Union[int, str]:
        """
        Underlying PartialBoolean value.
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

    def __repr__(self) -> str:

        return f"PartialBoolean({self._value!r})"

    def __str__(self) -> str:

        return str(self._value)

    def __bool__(self) -> bool:

        if self.is_free:
            raise ValueError("cannot convert free PartialBoolean to bool")

        return bool(self._value)

    def __eq__(self, other: object) -> bool:

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

        if value in [0, 1, "*"]:
            return value

        raise ValueError(
            "unsupported PartialBoolean value: "
            f"expected 0, 1, False, True or '*', "
            f"but received {value!r}"
        )
