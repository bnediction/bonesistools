#!/usr/bin/env python

import importlib
from typing import (
    Optional,
    Union,
    Callable,
    Sequence,
    List,
    Tuple
)
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal # type: ignore

from anndata import AnnData
from pandas import DataFrame
import numpy as np

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
Keys = Union[str, Sequence[str]]

AnnDataList = List[AnnData]

if _mudata_is_available:
    from mudata import MuData # type: ignore
    MuDataList = List[MuData]
    ScData = Union[AnnData, MuData]
else:
    MuData = type(NotImplemented)
    MuDataList = List[type(NotImplemented)]
    ScData = AnnData

Metric = Literal[
    "cityblock",
    "cosine",
    "euclidean",
    "l1",
    "l2",
    "manhattan",
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule"
]
Metric_Function = Callable[[np.ndarray, np.ndarray], float]

Shortest_Path_Method = Literal["dijkstra", "bellman-ford"]

def type_checker(
    function: Optional[Callable] = None,
    **options
):
    """
    Check if the argument types passed to 'function' are correct.
    
    Parameters
    ----------
    function: Callable
        function for which argument types are checked.
    **options
        key-value data-structure where 'key' is the argument name
        and 'value' is the expected argument type.
    
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
                        raise TypeError(f"unsupported argument type for key '{key}': expected '{value}' but received {type(key)}")
                elif not isinstance(arg, value):
                    raise TypeError(f"unsupported argument type for key '{key}': expected {value.__name__} but received {type(key)}")
            output = function(*args, **kwargs)
            return output
        return wrapper
    
    else:

        @wraps(function)
        def partial_wrapper(function):
            return type_checker(function, **options)

        return partial_wrapper

def anndata_checker(
    function: Optional[Callable] = None,
    n: int = 1
):
    """
    Check if the first arguments are AnnData instances.
    
    Parameters
    ----------
    function: Callable
        Function for which the first argument types are expected being AnnData instances.
    n: int (default: 1)
        Number of arguments to test.
    
    Returns
    -------
    Raise an error if at least one of the first 'n' arguments is not a AnnData instance.
    """
    
    if function is not None:

        @wraps(function)
        def wrapper(*args, **kwargs):
            iterator = iter(list(args) + list(kwargs.values()))
            for i in range(n):
                value = next(iterator)
                if not isinstance(value, AnnData):
                    raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected {AnnData} but received {type(value)}")
            return function(*args, **kwargs)
        return wrapper

    else:

        @wraps(function)
        def partial_wrapper(function):
            return anndata_checker(function, n)
        return partial_wrapper

if _mudata_is_available:

    def mudata_checker(
        function: Optional[Callable] = None,
        n: int = 1
    ):
        """
        Check if the first arguments are MuData instances.
    
        Parameters
        ----------
        function: Callable
            Function for which the first argument types are expected being MuData instances.
        n: int (default: 1)
            Number of arguments to test.
    
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
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected {MuData} but received {type(value)}")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper
    
    def anndata_or_mudata_checker(
        function: Optional[Callable] = None,
        n: int = 1
    ):
        """
        Check if the first arguments are AnnData or MuData instances.
    
        Parameters
        ----------
        function: Callable
            Function for which the first argument types are expected being AnnData or MuData instances.
        n: int (default: 1)
            Number of arguments to test.
    
        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not a AnnData or MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not (isinstance(value, AnnData) or isinstance(value, MuData)):
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected {AnnData} or {MuData} but received {type(value)}")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper

else:

    def mudata_checker(
        function: Optional[Callable] = None,
        n: int = 1
    ):
        """
        Check if the first arguments are MuData instances.
    
        Parameters
        ----------
        function: Callable
            Function for which the first argument types are expected being MuData instances.
        n: int (default: 1)
            Number of arguments to test.
    
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
        function: Optional[Callable] = None,
        n: int = 1
    ):
        """
        Check if the first arguments are AnnData or MuData instances.
    
        Parameters
        ----------
        function: Callable
            Function for which the first argument types are expected being AnnData or MuData instances.
        n: int (default: 1)
            Number of arguments to test.
    
        Returns
        -------
        Raise an error if at least one of the first 'n' arguments is not a AnnData or MuData instance.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not isinstance(value, AnnData):
                        raise TypeError(f"unsupported argument type for '{i+1}'-th argument: expected {AnnData} but received {type(value)}")
                return function(*args, **kwargs)
            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)
            return partial_wrapper
