#!/usr/bin/env python

from __future__ import annotations

from collections import Counter
from numbers import Number
from typing import Any, Iterable, Mapping, Optional, Tuple

from ._boolean import PartialBoolean
from ._typing import HypercubeLike, PartialBooleanLike


class PartialBooleanDifferential:
    """
    Discrete differential operator over partial Boolean values.

    This utility class defines a modest discrete differential calculus over
    partial Boolean values by using the biological ordering convention:

        0 < * < 1

    The differential between two partial Boolean values is defined as:

        d(x1, x2) =  1  if x1 < x2
                     0  if x1 = x2
                    -1  if x1 > x2

    This should be understood as a discrete directional comparison, not as a
    continuous differential calculus. The ordering is intentionally distinct
    from the ensemble/set-theoretic interpretation of `PartialBoolean.contains`,
    where "*" represents the Boolean ensemble {0, 1}.

    Examples
    --------
    >>> PartialBooleanDifferential.differential(0, 1)
    1

    >>> PartialBooleanDifferential.differential(1, 0)
    -1

    >>> PartialBooleanDifferential.differential(0, "*")
    1

    >>> PartialBooleanDifferential.differential("*", 1)
    1

    >>> PartialBooleanDifferential.differential("*", "*")
    0
    """

    def __init__(self) -> None:
        raise TypeError(
            "PartialBooleanDifferential is a static utility class "
            "and should not be instantiated."
        )

    @staticmethod
    def differential(v1: PartialBooleanLike, v2: PartialBooleanLike) -> int:
        """
        Return the discrete differential between two partial Boolean values.

        Values are compared according to the biological ordering convention:

            0 < * < 1

        Parameters
        ----------
        v1, v2: PartialBooleanLike
            Values compared after conversion to PartialBoolean.

        Returns
        -------
        int
            Differential from `v1` to `v2`:
                - 1 if `v1 < v2`,
                - 0 if `v1 == v2`,
                - -1 if `v1 > v2`.
        """

        v1 = PartialBooleanDifferential._coerce_value(v1)
        v2 = PartialBooleanDifferential._coerce_value(v2)

        if v1 == v2:
            return 0

        return 1 if v1 < v2 else -1

    @staticmethod
    def _coerce_value(value: PartialBooleanLike) -> PartialBoolean:
        """
        Convert a value to a PartialBoolean.

        Parameters
        ----------
        value: PartialBooleanLike
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

        try:
            return value if isinstance(value, PartialBoolean) else PartialBoolean(value)

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
                f"expected bool, number, str or PartialBoolean but received "
                f"{type(value)}"
            )


class BooleanPredecessorInference:
    """
    Infer predecessor relationships between partial Boolean configurations.

    This class implements the scBoolDiff predecessor logic. It combines:

        1. a discrete partial Boolean differential, based on:

               0 < * < 1

        2. signed regulatory interactions of the form:

               source --sign--> target

    to infer whether one partial Boolean configuration may precede another one.

    The method is intended as a regulatory-consistency heuristic. It should not
    be interpreted as a full probabilistic trajectory inference method. Its goal
    is to count whether observed Boolean changes are consistent with signed
    regulatory effects.

    For one signed interaction `(source, target, sign)`, a conclusion can be
    made only when:
        - the source has the same value in both configurations,
        - the source is fixed, i.e. not `"*"`,
        - the target changes between the two configurations.

    Examples
    --------
    >>> cell1 = {"A": 1, "B": 0, "C": 1, "D": 1}
    >>> cell2 = {"A": 1, "B": 1, "C": 1, "D": 0}

    >>> interactions = [
    ...     ("A", "B", {"sign": 1}),
    ...     ("C", "D", {"sign": -1}),
    ... ]

    The first interaction says that active `A` positively regulates `B`.
    Since `A` is active in both cells and `B` increases from 0 to 1, this
    supports `cell1` preceding `cell2`.

    The second interaction says that active `C` negatively regulates `D`.
    Since `C` is active in both cells and `D` decreases from 1 to 0, this
    also supports `cell1` preceding `cell2`.

    >>> BooleanPredecessorInference.predecessor(cell1, cell2, interactions)
    'first_precedes_second'
    """

    FIRST_PRECEDES_SECOND = "first_precedes_second"
    SECOND_PRECEDES_FIRST = "second_precedes_first"
    INCONCLUSIVE = "inconclusive"

    def __init__(self) -> None:
        raise TypeError(
            "BooleanPredecessorInference is a static utility class "
            "and should not be instantiated."
        )

    @staticmethod
    def pairwise_predecessor_test(
        source_v1: PartialBooleanLike,
        source_v2: PartialBooleanLike,
        target_v1: PartialBooleanLike,
        target_v2: PartialBooleanLike,
        sign: int,
    ) -> Optional[bool]:
        """
        Infer pairwise predecessor order from signed regulatory consistency.

        Parameters
        ----------
        source_v1, source_v2: bool, number, str or PartialBoolean
            Source-node partial Boolean values in configurations 1 and 2.
        target_v1, target_v2: bool, number, str or PartialBoolean
            Target-node partial Boolean values in configurations 1 and 2.
        sign: {-1, 1}
            Sign of the source influence on the target.

        Returns
        -------
        bool or None
            `True` if configuration 1 is inferred to precede configuration 2,
            `False` if configuration 2 is inferred to precede configuration 1,
            and `None` if no conclusion can be made.

        Raises
        ------
        ValueError
            If `sign` is not -1 or 1.
        """

        if sign not in {-1, 1}:
            raise ValueError(
                "invalid argument value for 'sign': "
                f"expected -1 or 1 but received {sign!r}"
            )

        source_v1 = BooleanPredecessorInference._coerce_value(source_v1)
        source_v2 = BooleanPredecessorInference._coerce_value(source_v2)
        target_v1 = BooleanPredecessorInference._coerce_value(target_v1)
        target_v2 = BooleanPredecessorInference._coerce_value(target_v2)

        source_differential = PartialBooleanDifferential.differential(
            source_v1,
            source_v2,
        )
        target_differential = PartialBooleanDifferential.differential(
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

    @staticmethod
    def predecessor_votes(
        cell1: HypercubeLike,
        cell2: HypercubeLike,
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]],
    ) -> Counter[str]:
        """
        Return predecessor votes between two partial Boolean configurations.

        Each signed interaction contributes one vote:
            - `first_precedes_second` if it supports `cell1` before `cell2`,
            - `second_precedes_first` if it supports the reverse order,
            - `inconclusive` if the interaction cannot support either order.

        Interactions involving genes missing from either configuration are
        counted as inconclusive.

        Parameters
        ----------
        cell1, cell2: HypercubeLike
            Partial Boolean configurations indexed by component name.
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]]
            Signed source-target interactions.

        Returns
        -------
        Counter[str]
            Vote counts for all predecessor categories.
        """

        votes = Counter(
            {
                BooleanPredecessorInference.FIRST_PRECEDES_SECOND: 0,
                BooleanPredecessorInference.SECOND_PRECEDES_FIRST: 0,
                BooleanPredecessorInference.INCONCLUSIVE: 0,
            }
        )

        for source, target, data in interactions:
            if source not in cell1 or source not in cell2:
                votes[BooleanPredecessorInference.INCONCLUSIVE] += 1
                continue

            if target not in cell1 or target not in cell2:
                votes[BooleanPredecessorInference.INCONCLUSIVE] += 1
                continue

            result = BooleanPredecessorInference.pairwise_predecessor_test(
                source_v1=cell1[source],
                source_v2=cell2[source],
                target_v1=cell1[target],
                target_v2=cell2[target],
                sign=BooleanPredecessorInference._edge_sign(data),
            )

            if result is True:
                votes[BooleanPredecessorInference.FIRST_PRECEDES_SECOND] += 1

            elif result is False:
                votes[BooleanPredecessorInference.SECOND_PRECEDES_FIRST] += 1

            else:
                votes[BooleanPredecessorInference.INCONCLUSIVE] += 1

        return votes

    @staticmethod
    def predecessor(
        cell1: HypercubeLike,
        cell2: HypercubeLike,
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]],
    ) -> Optional[str]:
        """
        Infer the most supported predecessor relationship between two configurations.

        Parameters
        ----------
        cell1, cell2: HypercubeLike
            Partial Boolean configurations indexed by component name.
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]]
            Signed source-target interactions.

        Returns
        -------
        str or None
            Most supported predecessor label, or None in case of a tie or no
            conclusive interaction.
        """

        votes = BooleanPredecessorInference.predecessor_votes(
            cell1,
            cell2,
            interactions,
        )

        forward = votes[BooleanPredecessorInference.FIRST_PRECEDES_SECOND]
        backward = votes[BooleanPredecessorInference.SECOND_PRECEDES_FIRST]

        if forward > backward:
            return BooleanPredecessorInference.FIRST_PRECEDES_SECOND

        if backward > forward:
            return BooleanPredecessorInference.SECOND_PRECEDES_FIRST

        return None

    @staticmethod
    def predecessor_score(
        cell1: HypercubeLike,
        cell2: HypercubeLike,
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]],
    ) -> float:
        """
        Return a normalized predecessor score between two configurations.

        The score is:

            (forward - backward) / conclusive

        where `conclusive = forward + backward`.

        Parameters
        ----------
        cell1, cell2: HypercubeLike
            Partial Boolean configurations indexed by component name.
        interactions: Iterable[Tuple[str, str, Mapping[str, Any]]]
            Signed source-target interactions.

        Returns
        -------
        float
            Normalized score in [-1, 1]. Returns 0.0 if no conclusive
            interaction is available.
        """

        votes = BooleanPredecessorInference.predecessor_votes(
            cell1,
            cell2,
            interactions,
        )

        forward = votes[BooleanPredecessorInference.FIRST_PRECEDES_SECOND]
        backward = votes[BooleanPredecessorInference.SECOND_PRECEDES_FIRST]
        conclusive = forward + backward

        if conclusive == 0:
            return 0.0

        return (forward - backward) / conclusive

    @staticmethod
    def _coerce_value(value: PartialBooleanLike) -> PartialBoolean:

        return PartialBooleanDifferential._coerce_value(value)

    @staticmethod
    def _edge_sign(data: Mapping[str, Any]) -> int:
        """
        Return and validate an interaction sign.

        Parameters
        ----------
        data: Mapping[str, Any]
            Edge attributes containing a `sign` value.

        Returns
        -------
        int
            Interaction sign, either -1 or 1.

        Raises
        ------
        ValueError
            If the sign is missing or is not -1 or 1.
        """

        sign = data.get("sign")

        if sign not in {-1, 1}:
            raise ValueError(
                "invalid interaction sign: expected -1 or 1 " f"but received {sign!r}"
            )

        return sign
