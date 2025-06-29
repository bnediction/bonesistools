#!/usr/bin/env python

from typing import Union

import math

class PartialBoolean:

    def __init__(self, value: Union[bool, float, int]):
        if isinstance(value, bool):
            self.__value = 1 if value is True else 0
        elif value in [0, 1] or math.isnan(value):
            self.__value = value
        else:
            raise ValueError(f"invalid argument value for 'value': {value}")
    
    def __repr__(self):
        return f"Boolean({self.__value})"
    
    def __str__(self):
        return f"{self.__value}"

    @property
    def get(self):
        return self.__value

    @property
    def set(self, value):
        if isinstance(value, bool):
            self.__value = 1 if value is True else 0
        elif value in [0, 1] or math.isnan(value):
            self.__value = value
        else:
            raise ValueError(f"invalid argument value for 'value': {value}")
    
    def __eq__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'==' not supported between instances of {PartialBoolean} and {type(other)}")
        elif math.isnan(self.__value):
            return True if math.isnan(other.__value) else False
        else:
            return self.__value == other.__value

    def __ne__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'!=' not supported between instances of {PartialBoolean} and {type(other)}")
        else:
            return not self == other
    
    def __lt__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'<' not supported between instances of {PartialBoolean} and {type(other)}")
        elif math.isnan(self.__value):
            if math.isnan(other.__value):
                return False
            else:
                return False if other.__value == 0 else True
        elif math.isnan(other.__value):
            return True if self.__value == 0 else False
        else:
            return self.__value < other.__value
    
    def __le__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'<=' not supported between instances of {PartialBoolean} and {type(other)}")
        elif math.isnan(self.__value):
            if math.isnan(other.__value):
                return True
            else:
                return False if other.__value == 0 else True
        elif math.isnan(other.__value):
            return True if self.__value == 0 else False
        else:
            return self.__value <= other.__value
    
    def __gt__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'>' not supported between instances of {PartialBoolean} and {type(other)}")
        elif math.isnan(self.__value):
            if math.isnan(other.__value):
                return False
            else:
                return True if other.__value == 0 else False
        elif math.isnan(other.__value):
            return False if self.__value == 0 else True
        else:
            return self.__value > other.__value
    
    def __ge__(self, other):
        if not isinstance(other, PartialBoolean):
            raise TypeError(f"'>=' not supported between instances of {PartialBoolean} and {type(other)}")
        elif math.isnan(self.__value):
            if math.isnan(other.__value):
                return True
            else:
                return True if other.__value == 0 else False
        elif math.isnan(other.__value):
            return False if self.__value == 0 else True
        else:
            return self.__value >= other.__value
