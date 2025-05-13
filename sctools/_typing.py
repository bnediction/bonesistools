#!/usr/bin/env python

import importlib
import types
from typing import (
    Optional,
    Union,
    List,
    Tuple
)
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from anndata import AnnData
from pandas import DataFrame

from functools import wraps

try:
    _mudata_is_available = importlib.util.find_spec("mudata") is not None
except:
    _mudata_is_available = importlib.find_loader("mudata") is not None

class UnionType(object):
    
    def __init__(self, *args):
        self.types = args

    def __str__(self):
        return ",".join(self.types)

    __repr__ = __str__

DataFrameList = List[DataFrame]
Suffixes = Tuple[Optional[str], Optional[str]]

AxisInt = int
Axis = Union[AxisInt, Literal["obs", "var"]]
Keys = Union[List[str],str]

AnnDataList = List[AnnData]

if _mudata_is_available:
    from mudata import MuData
    MuDataList = List[MuData]
    ScData = Union[AnnData,MuData]
else:
    MuData = type(NotImplemented)
    MuDataList = List[type(NotImplemented)]
    ScData = AnnData

def type_checker(
    function: types.FunctionType=None,
    **options
):
    """
    Check if the argument types passed to function are correct.
    
    Parameters
    ----------
    function
        called function
    **options
        key-value data-structure where 'key' is the argument name
        and 'value' is the expected argument type
    
    Returns
    -------
    Raise an error if at least one argument type is not correct.
    """

    if function is not None:

        @wraps(function)
        def wrapper(*args, **kwargs):
            if len(options) == 0:
                raise Exception("Expected verification arguments")
            function_code = function.__code__
            arg_names = function_code.co_varnames
            for key, value in options.items():
                idx = arg_names.index(key)
                if (len(args) > idx):
                    arg = args[idx]
                else:
                    if key in kwargs:
                        arg = kwargs.get(key)
                if isinstance(value, UnionType):
                    types_match = False
                    for dtype in value.types:
                        if isinstance(arg, dtype):
                            types_match = True
                    if types_match == False:
                        raise TypeError(f"unsupported argument type for '{key}': expected '{value}' but received '{type(key)}'")
                elif not isinstance(arg, value):
                    raise TypeError(f"unsupported argument type for '{key}': expected '{value.__name__}' but received '{type(key)}'")
            output = function(*args, **kwargs)
            return output
        return wrapper
    
    else:

        @wraps(function)
        def partial_wrapper(function):
            return type_checker(function, **options)

        return partial_wrapper

def anndata_checker(
    function: types.FunctionType=None,
    n: int=1
):
    """
    Check if the first 'n' arguments are AnnData instances.
    
    Parameters
    ----------
    function
        called function
    n
        number of arguments to test
    
    Returns
    -------
    Raise an error if at least one of the first 'n' arguments is not an AnnData instance.
    """
    
    if function is not None:

        @wraps(function)
        def wrapper(*args, **kwargs):
            iterator = iter(list(args) + list(kwargs.values()))
            for i in range(n):
                value = next(iterator)
                if not isinstance(value, AnnData):
                    raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected '{AnnData}' but received '{type(value)}'")
            return function(*args, **kwargs)
        return wrapper

    else:

        @wraps(function)
        def partial_wrapper(function):
            return anndata_checker(function, n)
        return partial_wrapper

if _mudata_is_available:

    def mudata_checker(
        function: types.FunctionType=None,
        n: int=1
    ):
        """
        Check if the first 'n' arguments are MuData instances.

        Parameters
        ----------
        function
            called function
        n
            number of arguments to test

        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not a MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not isinstance(value, MuData):
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected '{MuData}' but received '{type(value)}'")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper
    
    def anndata_or_mudata_checker(
        function: types.FunctionType=None,
        n: int=1
    ):
        """
        Check if the first 'n' arguments are AnnData or MuData instances.

        Parameters
        ----------
        function
            called function
        n
            number of arguments to test

        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not an AnnData or a MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not (isinstance(value, AnnData) or isinstance(value, MuData)):
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected '{AnnData}' or '{MuData}' but received '{type(value)}'")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper

else:

    def mudata_checker(
        function: types.FunctionType=None,
        n: int=1
    ):
        """
        Check if the first 'n' arguments are MuData instances.
        
        Parameters
        ----------
        function
            called function
        n
            number of arguments to test
        
        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not a MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                raise ModuleNotFoundError("no module named 'mudata'")
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper

    def anndata_or_mudata_checker(
        function: types.FunctionType=None,
        n: int=1
    ):
        """
        Check if the first 'n' arguments are AnnData or MuData instances.

        Parameters
        ----------
        function
            called function
        n
            number of arguments to test

        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not an AnnData or a MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not isinstance(value, AnnData):
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected '{AnnData}' but received '{type(value)}'")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper
