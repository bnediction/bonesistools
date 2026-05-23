#!/usr/bin/env python

from numbers import Number
from typing import Optional, Union

import math

from ._boolean import PartialBoolean


class BooleanDifferentialCalculus:
    """
    Boolean differential calculus over partial Boolean values.

    This class provides utilities for comparing partial Boolean values under the
    biological ordering convention:

        0 < * < 1

    This ordering is used to compute Boolean differentials between two partial
    Boolean states and to infer pairwise predecessor relationships from signed
    regulatory effects.

    Values are converted to `PartialBoolean` before comparison. Numeric NaN
    values are interpreted as the free partial value `"*"`.

    Raises
    ------
    TypeError
        If a value cannot be converted to a PartialBoolean because its type is
        unsupported.
    ValueError
        If a numeric value is not a valid PartialBoolean value, or if `sign`
        is not -1 or 1.
    AssertionError
        If pairwise predecessor inference reaches an inconsistent internal
        state.
    """

    @staticmethod
    def _coerce_value(value: Union[bool, Number, PartialBoolean]) -> PartialBoolean:
        """
        Convert a value to a PartialBoolean.

        Parameters
        ----------
        value: bool, number or PartialBoolean
            Value to convert. Supported values are `False`, `True`, `0`, `1`,
            `NaN`, `"*"` and PartialBoolean values.

        Returns
        -------
        PartialBoolean
            Converted partial Boolean value.

        Raises
        ------
        TypeError
            If `value` has an unsupported type.
        ValueError
            If `value` is numeric but not 0, 1 or NaN.
        """

        if isinstance(value, PartialBoolean):
            return value

        if isinstance(value, bool):
            return PartialBoolean(value)

        if isinstance(value, Number):
            if value in [0, 1] or (
                isinstance(value, float)
                and math.isnan(value)
            ):
                return PartialBoolean(value)

            raise ValueError(
                "invalid argument value for 'value': "
                f"expected 0, 1, False, True or NaN but received {value!r}"
            )

        if value == "*":
            return PartialBoolean(value)

        raise TypeError(
            "unsupported argument type for 'value': "
            f"expected bool, number, '*' or PartialBoolean but received "
            f"{type(value)}"
        )

    @staticmethod
    def differential(v1, v2) -> int:
        """
        Return the Boolean differential between two partial Boolean values.

        Values are compared according to the biological ordering convention:

            0 < * < 1

        The returned differential is:
            - 0 if both values are equal,
            - 1 if `v1 < v2`,
            - -1 if `v1 > v2`.

        Examples
        --------
        >>> BooleanDifferentialCalculus.differential(0, 1)
        1

        >>> BooleanDifferentialCalculus.differential(1, 0)
        -1

        >>> BooleanDifferentialCalculus.differential(0, "*")
        1

        >>> BooleanDifferentialCalculus.differential("*", 1)
        1

        >>> BooleanDifferentialCalculus.differential("*", "*")
        0

        Parameters
        ----------
        v1, v2: bool, number, str or PartialBoolean
            Values to compare after conversion to PartialBoolean.

        Returns
        -------
        int
            Boolean differential from `v1` to `v2`.
        """

        v1 = BooleanDifferentialCalculus._coerce_value(v1)
        v2 = BooleanDifferentialCalculus._coerce_value(v2)

        if v1 == v2:
            return 0

        return 1 if v1 < v2 else -1

    @staticmethod
    def pairwise_predecessor_test(
        source_v1: Union[bool, Number, str, PartialBoolean],
        source_v2: Union[bool, Number, str, PartialBoolean],
        target_v1: Union[bool, Number, str, PartialBoolean],
        target_v2: Union[bool, Number, str, PartialBoolean],
        sign: int,
    ) -> Optional[bool]:
        """
        Infer pairwise predecessor order from signed regulatory consistency.

        The test compares two observations, denoted condition 1 and condition 2,
        using:
            - the partial Boolean state of a source node,
            - the partial Boolean state of a target node,
            - the signed influence of the source on the target.

        A predecessor relationship can be inferred only when the source state is
        fixed and identical in both observations. If the source differs between
        observations, or if the source is unspecified, no conclusion is made.

        Parameters
        ----------
        source_v1, source_v2: bool, number, str or PartialBoolean
            Source-node partial Boolean values in observations 1 and 2.
        target_v1, target_v2: bool, number, str or PartialBoolean
            Target-node partial Boolean values in observations 1 and 2.
        sign: {-1, 1}
            Sign of the source influence on the target.

        Returns
        -------
        bool or None
            `True` if observation 1 is inferred to precede observation 2,
            `False` if observation 2 is inferred to precede observation 1,
            and `None` if no conclusion can be made.

        Raises
        ------
        ValueError
            If `sign` is not -1 or 1.
        """

        if sign not in [-1, 1]:
            raise ValueError(
                "invalid argument value for 'sign': "
                f"expected -1 or 1 but received {sign!r}"
            )

        source_v1 = BooleanDifferentialCalculus._coerce_value(source_v1)
        source_v2 = BooleanDifferentialCalculus._coerce_value(source_v2)
        target_v1 = BooleanDifferentialCalculus._coerce_value(target_v1)
        target_v2 = BooleanDifferentialCalculus._coerce_value(target_v2)

        source_differential = BooleanDifferentialCalculus.differential(
            source_v1,
            source_v2,
        )
        target_differential = BooleanDifferentialCalculus.differential(
            target_v1,
            target_v2,
        )

        if source_differential != 0 or target_differential == 0:
            return None

        if source_v1 == source_v2 == PartialBoolean("*"):
            return None

        if source_v1 == source_v2 == PartialBoolean(1):
            return sign == target_differential

        if source_v1 == source_v2 == PartialBoolean(0):
            return sign != target_differential

        raise AssertionError("found inconsistent predecessor inference state")