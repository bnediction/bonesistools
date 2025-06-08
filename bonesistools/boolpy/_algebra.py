#!/usr/bin/env python

from typing import Union
from numbers import Number
from ._boolean import PartialBoolean

import math

class BooleanDifferentialCalculus(object):

    def __init__(self) -> None:
        pass

    def __conversion__(self, value):
        if isinstance(value, bool):
            return PartialBoolean(value)
        elif isinstance(value, Number):
            if value in [0, 1] or math.isnan(value):
                return PartialBoolean(value)
            else:
                raise TypeError(f"unsupported type conversion for '{value}'")
        elif isinstance(value, PartialBoolean):
            return value
        else:
            raise TypeError(f"unsupported type conversion for '{value}'")

    def differential(self, v1, v2):
        _v1 = self.__conversion__(v1)
        _v2 = self.__conversion__(v2)
        if _v1 == _v2:
            return 0
        elif _v1 < _v2:
            return 1
        elif _v1 > _v2:
            return -1
        
    def pairwise_predecessor_test(
        self,
        source_v1: Union[bool, PartialBoolean],
        source_v2: Union[bool, PartialBoolean],
        target_v1: Union[bool, PartialBoolean],
        target_v2: Union[bool, PartialBoolean],
        sign: int
    ) -> Union[bool, None]:
        """
        By assuming there is two conditions 1 and 2, estimate which one preceeds the other one with respect to two nodes.

        Parameters
        ----------
        source_v1, source_v2, target_v1, target_v2: bool | PartialBoolean
            Partial Boolean values.
        sign: -1 | 1
            Specify the sign effect of source upon target.

        Returns
        -------
        Return True if first condition preceeds the second one, False in the contrary case and None if no conclusion is possible.
        """

        source_v1 = self.__conversion__(source_v1)
        source_v2 = self.__conversion__(source_v2)
        target_v1 = self.__conversion__(target_v1)
        target_v2 = self.__conversion__(target_v2)
        source_differential = self.differential(source_v1, source_v2)
        target_differential = self.differential(target_v1, target_v2)
        if sign not in [-1, 1]:
            raise ValueError(f"invalid argument value for 'sign': {sign}")
        if source_differential != 0 or target_differential == 0:
            return None
        elif source_differential == 0:
            if source_v1 == source_v2 == PartialBoolean(1):
                return True if sign == target_differential else False
            if source_v1 == source_v2 == PartialBoolean(0):
                return False if sign == target_differential else True
            if source_v1 == source_v2 == PartialBoolean(float("nan")):
                return None
        else:
            raise AssertionError("found incoherence")
