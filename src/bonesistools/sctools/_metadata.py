#!/usr/bin/env python

import numbers
from typing import Optional, Union

import numpy as np

from .._typing import RandomStateSeed


def _format_random_state(seed: RandomStateSeed) -> Optional[Union[int, str]]:
    if seed is None:
        return None
    if seed is np.random:
        return "np.random"
    if isinstance(seed, numbers.Integral):
        return int(seed)
    if isinstance(seed, np.random.RandomState):
        return "RandomState"
    return str(seed)
