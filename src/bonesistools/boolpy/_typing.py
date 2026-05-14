#!/usr/bin/env python

import importlib

from typing import (
    Union,
)

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

try:
    from mpbn import MPBooleanNetwork
    _mpbn_is_available = True
except ImportError:
    _mpbn_is_available = False
    MPBooleanNetwork = type(NotImplemented)

AxisInt = int
Axis = Union[AxisInt, Literal["obs", "var"]]
