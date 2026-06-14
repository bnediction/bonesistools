#!/usr/bin/env python

import numbers
from typing import Callable, Iterable, Optional, TypeVar, Union, overload

import numpy as np

from ._compat import Literal
from ._typing import (
    DataFrameAxis,
    FileOrientation,
    RandomStateSeed,
)

T = TypeVar("T", bound=str)
F = TypeVar("F", bound=Callable[..., object])


def _format_choices(choices: Iterable[str], *, include_none: bool = False) -> str:
    choices = tuple(choices)
    formatted = tuple(repr(choice) for choice in choices)
    if include_none:
        formatted = formatted + ("None",)
    if len(formatted) == 1:
        return formatted[0]
    if len(formatted) == 2:
        return f"{formatted[0]} or {formatted[1]}"
    return f"{', '.join(formatted[:-1])} or {formatted[-1]}"


@overload
def _as_literal(
    value: str,
    *,
    choices: Iterable[T],
    name: str,
    allow_none: Literal[False] = False,
) -> T: ...


@overload
def _as_literal(
    value: Optional[str],
    *,
    choices: Iterable[T],
    name: str,
    allow_none: Literal[True],
) -> Optional[T]: ...


def _as_literal(
    value: Optional[str],
    *,
    choices: Iterable[T],
    name: str,
    allow_none: bool = False,
) -> Optional[T]:
    if value is None and allow_none:
        return None

    if not isinstance(value, str):
        expected = f"{str} or None" if allow_none else str
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {expected} but received {type(value)}"
        )

    resolved_choices = tuple(choices)
    for choice in resolved_choices:
        if value == choice:
            return choice

    expected = _format_choices(resolved_choices, include_none=allow_none)
    raise ValueError(
        f"invalid argument value for '{name}': "
        f"expected one of {expected} but received {value!r}"
    )


@overload
def _as_string(
    value: str,
    name: str,
    *,
    allow_none: Literal[False] = False,
) -> str: ...


@overload
def _as_string(
    value: Optional[str],
    name: str,
    *,
    allow_none: Literal[True],
) -> Optional[str]: ...


def _as_string(
    value: Optional[str],
    name: str,
    *,
    allow_none: bool = False,
) -> Optional[str]:
    if value is None and allow_none:
        return None

    if not isinstance(value, str):
        expected = f"{str} or None" if allow_none else str
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {expected} but received {type(value)}"
        )

    return value


@overload
def _as_boolean(
    value: bool,
    name: str,
    *,
    allow_none: Literal[False] = False,
) -> bool: ...


@overload
def _as_boolean(
    value: Optional[bool],
    name: str,
    *,
    allow_none: Literal[True],
) -> Optional[bool]: ...


def _as_boolean(
    value: Optional[bool],
    name: str,
    *,
    allow_none: bool = False,
) -> Optional[bool]:
    if value is None and allow_none:
        return None

    if not isinstance(value, bool):
        expected = f"{bool} or None" if allow_none else bool
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {expected} but received {type(value)}"
        )

    return value


@overload
def _as_callable(
    value: F,
    name: str,
    *,
    allow_none: Literal[False] = False,
) -> F: ...


@overload
def _as_callable(
    value: Optional[F],
    name: str,
    *,
    allow_none: Literal[True],
) -> Optional[F]: ...


def _as_callable(
    value: Optional[F],
    name: str,
    *,
    allow_none: bool = False,
) -> Optional[F]:
    if value is None and allow_none:
        return None

    if not callable(value):
        expected = "callable object or None" if allow_none else "callable object"
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {expected} but received {type(value)}"
        )

    return value


def _as_positive_number(value: float, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {float} but received {type(value)}"
        )

    value = float(value)
    if value <= 0:
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected positive value but received {value!r}"
        )

    return value


def _as_non_negative_number(value: float, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {float} but received {type(value)}"
        )

    value = float(value)
    if value < 0:
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected non-negative value but received {value!r}"
        )

    return value


def _as_positive_integer(value: Union[int, float], name: str) -> int:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {int} but received {type(value)}"
        )

    if value <= 0:
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected non-null positive value but received {value!r}"
        )

    if isinstance(value, float) and not value.is_integer():
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected integer but received {value!r}"
        )

    return int(value)


def _as_non_negative_integer(value: Union[int, float], name: str) -> int:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {int} but received {type(value)}"
        )

    if value < 0:
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected non-negative value but received {value!r}"
        )

    if isinstance(value, float) and not value.is_integer():
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected integer but received {value!r}"
        )

    return int(value)


def _as_probability(value: float, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(
            f"unsupported argument type for '{name}': "
            f"expected {float} but received {type(value)}"
        )

    value = float(value)
    if not 0 <= value <= 1:
        raise ValueError(
            f"invalid argument value for '{name}': "
            f"expected value between 0 and 1 but received {value!r}"
        )

    return value


def _as_seed(seed: RandomStateSeed) -> np.random.RandomState:
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, numbers.Integral):
        return np.random.RandomState(int(seed))
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        f"invalid argument value for 'seed': "
        f"expected None, np.random, integer seed or np.random.RandomState "
        f"but received {seed!r}"
    )


def _as_dataframe_axis(axis: Union[int, str]) -> DataFrameAxis:
    if isinstance(axis, str) and axis == "index":
        return "index"
    if isinstance(axis, str) and axis == "columns":
        return "columns"
    if isinstance(axis, int) and not isinstance(axis, bool) and axis == 0:
        return "index"
    if isinstance(axis, int) and not isinstance(axis, bool) and axis == 1:
        return "columns"

    raise ValueError(
        f"invalid argument value for 'axis': "
        f"expected 0, 1, 'index' or 'columns' but received {axis!r}"
    )


def _as_orientation(orientation: str) -> FileOrientation:
    if orientation == "rows" or orientation == "columns":
        return orientation

    raise ValueError(
        f"invalid argument value for 'orientation': "
        f"expected 'rows' or 'columns' but received {orientation!r}"
    )
