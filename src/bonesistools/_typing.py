#!/usr/bin/env python

from types import ModuleType
from typing import Optional, Union

import numpy as np

from ._compat import Literal

RandomStateSeed = Optional[Union[int, np.random.RandomState, ModuleType]]
AutoInteger = Union[int, Literal["auto"]]
DataFrameAxis = Literal["index", "columns"]
FileOrientation = Literal["rows", "columns"]
