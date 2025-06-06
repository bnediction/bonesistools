#!/usr/bin/env python

from typing import (
    Union,
)
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

AxisInt = int
Axis = Union[AxisInt, Literal["obs", "var"]]
