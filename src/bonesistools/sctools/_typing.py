#!/usr/bin/env python

import importlib
from typing import (
    Callable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from anndata import AnnData
from pandas import DataFrame
import numpy as np

from functools import wraps

try:
    _mudata_is_available = importlib.util.find_spec("mudata") is not None
except AttributeError:
    _mudata_is_available = importlib.find_loader("mudata") is not None


class UnionType(object):
    """
    Container used by `type_checker` to accept several possible types.
    """

    def __init__(self, *args):
        self.types = args

    def __str__(self):
        return " or ".join(
            getattr(dtype, "__name__", str(dtype)) for dtype in self.types
        )

    __repr__ = __str__


DataFrameList = List[DataFrame]
Suffixes = Tuple[Optional[str], Optional[str]]

AxisInt = int
Axis = Union[AxisInt, Literal["obs", "var"]]
Keys = Union[str, Sequence[str]]

AnnDataList = List[AnnData]

if _mudata_is_available:
    from mudata import MuData  # type: ignore

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
    "yule",
]
Metric_Function = Callable[[np.ndarray, np.ndarray], float]

Shortest_Path_Method = Literal["dijkstra", "bellman-ford"]


def type_checker(function: Optional[Callable] = None, **options):
    """
    Decorate a function with runtime argument type checks.

    Parameters
    ----------
    function: Callable, optional
        Function to decorate. If None, return a decorator.
    **options: type
        Mapping from argument names to expected runtime types.

    Returns
    -------
    Callable
        Decorated function, or decorator if `function` is None.

    Raises
    ------
    Exception
        If no verification options are provided.
    TypeError
        If one of the checked arguments has an unsupported type.
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
                if len(args) > idx:
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
                        raise TypeError(
                            f"unsupported argument type for '{key}': "
                            f"expected {value} but received {type(arg)}"
                        )
                elif not isinstance(arg, value):
                    raise TypeError(
                        f"unsupported argument type for '{key}': "
                        f"expected {value.__name__} but received {type(arg)}"
                    )
            output = function(*args, **kwargs)
            return output

        return wrapper

    else:

        @wraps(function)
        def partial_wrapper(function):
            return type_checker(function, **options)

        return partial_wrapper


def anndata_checker(function: Optional[Callable] = None, n: int = 1):
    """
    Decorate a function by checking that the first arguments are AnnData objects.

    Parameters
    ----------
    function: Callable, optional
        Function to decorate. If None, return a decorator.
    n: int (default: 1)
        Number of arguments to test.

    Returns
    -------
    Callable
        Decorated function, or decorator if `function` is None.

    Raises
    ------
    TypeError
        If one of the first `n` arguments is not an AnnData object.
    """

    if function is not None:

        @wraps(function)
        def wrapper(*args, **kwargs):
            iterator = iter(list(args) + list(kwargs.values()))
            for i in range(n):
                value = next(iterator)
                if not isinstance(value, AnnData):
                    raise TypeError(
                        f"unsupported argument type for '{i+1}'-th argument: "
                        f"expected {AnnData} but received {type(value)}"
                    )
            return function(*args, **kwargs)

        return wrapper

    else:

        @wraps(function)
        def partial_wrapper(function):
            return anndata_checker(function, n)

        return partial_wrapper


if _mudata_is_available:

    def mudata_checker(function: Optional[Callable] = None, n: int = 1):
        """
        Decorate a function by checking that the first arguments are MuData objects.

        Parameters
        ----------
        function: Callable, optional
            Function to decorate. If None, return a decorator.
        n: int (default: 1)
            Number of arguments to test.

        Returns
        -------
        Callable
            Decorated function, or decorator if `function` is None.

        Raises
        ------
        TypeError
            If one of the first `n` arguments is not a MuData object.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not isinstance(value, MuData):
                        raise TypeError(
                            f"unsupported argument type for '{i+1}'-th argument: "
                            f"expected {MuData} but received {type(value)}"
                        )
                return function(*args, **kwargs)

            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)

            return partial_wrapper

    def anndata_or_mudata_checker(function: Optional[Callable] = None, n: int = 1):
        """
        Decorate a function by checking that the first arguments are AnnData or MuData objects.

        Parameters
        ----------
        function: Callable, optional
            Function to decorate. If None, return a decorator.
        n: int (default: 1)
            Number of arguments to test.

        Returns
        -------
        Callable
            Decorated function, or decorator if `function` is None.

        Raises
        ------
        TypeError
            If one of the first `n` arguments is neither an AnnData nor a
            MuData object.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not (isinstance(value, AnnData) or isinstance(value, MuData)):
                        raise TypeError(
                            f"unsupported argument type for '{i+1}'-th argument: "
                            f"expected {AnnData} or {MuData} but received "
                            f"{type(value)}"
                        )
                return function(*args, **kwargs)

            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)

            return partial_wrapper

else:

    def mudata_checker(function: Optional[Callable] = None, n: int = 1):
        """
        Decorate a function by checking that the first arguments are MuData objects.

        Parameters
        ----------
        function: Callable, optional
            Function to decorate. If None, return a decorator.
        n: int (default: 1)
            Number of arguments to test.

        Returns
        -------
        Callable
            Decorated function, or decorator if `function` is None.

        Raises
        ------
        ModuleNotFoundError
            If the decorated function is called while mudata is not installed.
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

    def anndata_or_mudata_checker(function: Optional[Callable] = None, n: int = 1):
        """
        Decorate a function by checking that the first arguments are AnnData or MuData objects.

        Parameters
        ----------
        function: Callable, optional
            Function to decorate. If None, return a decorator.
        n: int (default: 1)
            Number of arguments to test.

        Returns
        -------
        Callable
            Decorated function, or decorator if `function` is None.

        Raises
        ------
        TypeError
            If one of the first `n` arguments is not an AnnData object.
        """

        if function is not None:

            @wraps(function)
            def wrapper(*args, **kwargs):
                iterator = iter(list(args) + list(kwargs.values()))
                for i in range(n):
                    value = next(iterator)
                    if not isinstance(value, AnnData):
                        raise TypeError(
                            f"unsupported argument type for '{i+1}'-th argument: "
                            f"expected {AnnData} but received {type(value)}"
                        )
                return function(*args, **kwargs)

            return wrapper

        else:

            @wraps(function)
            def partial_wrapper(function):
                return mudata_checker(function, n)

            return partial_wrapper
