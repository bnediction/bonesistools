#!/usr/bin/env python

import importlib as _importlib
from functools import wraps
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Collection,
    List,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    overload,
)

import numpy as np
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import spmatrix
from typing_extensions import ParamSpec, TypeAlias

from bonesistools._compat import Literal

P = ParamSpec("P")
R = TypeVar("R")

_importlib_util = getattr(_importlib, "util", None)

if _importlib_util is not None:
    _mudata_is_available = _importlib_util.find_spec("mudata") is not None
else:
    _find_loader = getattr(_importlib, "find_loader", None)
    _mudata_is_available = (
        False if _find_loader is None else _find_loader("mudata") is not None
    )


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


DataFrameList: TypeAlias = List[DataFrame]
Suffixes: TypeAlias = Tuple[Optional[str], Optional[str]]

AnnDataAxis: TypeAlias = Literal["obs", "var"]
AnnDataAxisWithBoth: TypeAlias = Literal["obs", "var", "both"]
AnnDataAxisWithInteger: TypeAlias = Union[AnnDataAxis, Literal[0, 1]]
Axis: TypeAlias = AnnDataAxis
Keys: TypeAlias = Union[str, Sequence[str]]

AnnDataList: TypeAlias = List[AnnData]
Matrix: TypeAlias = Union[np.ndarray, spmatrix]
VarSubset: TypeAlias = Optional[Union[str, Collection[str]]]

if TYPE_CHECKING:
    MuData = object
    MuDataList: TypeAlias = List[Any]
    ScData: TypeAlias = Any
else:
    if _mudata_is_available:
        from mudata import MuData  # type: ignore

        MuDataList = List[MuData]
        ScData = Union[AnnData, MuData]
    else:
        MuData = type(NotImplemented)
        MuDataList = List[type(NotImplemented)]
        ScData = AnnData

Metric: TypeAlias = Literal[
    "braycurtis",
    "canberra",
    "chebyshev",
    "cityblock",
    "correlation",
    "cosine",
    "dice",
    "euclidean",
    "hamming",
    "haversine",
    "jaccard",
    "l1",
    "l2",
    "mahalanobis",
    "manhattan",
    "matching",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalsneath",
    "sqeuclidean",
    "wminkowski",
    "yule",
]
UMAPMetric: TypeAlias = Literal[
    "braycurtis",
    "canberra",
    "chebyshev",
    "cityblock",
    "correlation",
    "cosine",
    "dice",
    "euclidean",
    "hamming",
    "haversine",
    "hellinger",
    "jaccard",
    "kulsinski",
    "l1",
    "l2",
    "ll_dirichlet",
    "mahalanobis",
    "manhattan",
    "matching",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "wminkowski",
    "yule",
]
Metric_Function: TypeAlias = Callable[[np.ndarray, np.ndarray], float]

Shortest_Path_Method: TypeAlias = Literal["dijkstra", "bellman-ford"]


@overload
def type_checker(
    function: Callable[P, R],
    **options: Any,
) -> Callable[P, R]: ...


@overload
def type_checker(
    function: None = None,
    **options: Any,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def type_checker(
    function: Optional[Callable[P, R]] = None,
    **options: Any,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
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
    """

    if function is not None:

        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            if len(options) == 0:
                raise Exception("Expected verification arguments")
            function_code = function.__code__
            arg_names = function_code.co_varnames
            for key, value in options.items():
                idx = arg_names.index(key)
                arg = None
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
                    if not types_match:
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

        def partial_wrapper(function: Callable[P, R]) -> Callable[P, R]:
            return type_checker(function, **options)

        return partial_wrapper


@overload
def anndata_checker(function: Callable[P, R], n: int = 1) -> Callable[P, R]: ...


@overload
def anndata_checker(
    function: None = None,
    n: int = 1,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def anndata_checker(
    function: Optional[Callable[P, R]] = None,
    n: int = 1,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
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

    """

    if function is not None:

        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
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

        def partial_wrapper(function: Callable[P, R]) -> Callable[P, R]:
            return anndata_checker(function, n)

        return partial_wrapper


@overload
def mudata_checker(function: Callable[P, R], n: int = 1) -> Callable[P, R]: ...


@overload
def mudata_checker(
    function: None = None,
    n: int = 1,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def mudata_checker(
    function: Optional[Callable[P, R]] = None,
    n: int = 1,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
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
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            if not _mudata_is_available:
                raise ModuleNotFoundError("no module named 'mudata'")

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

    def partial_wrapper(function: Callable[P, R]) -> Callable[P, R]:
        return mudata_checker(function, n)

    return partial_wrapper


@overload
def anndata_or_mudata_checker(
    function: Callable[P, R],
    n: int = 1,
) -> Callable[P, R]: ...


@overload
def anndata_or_mudata_checker(
    function: None = None,
    n: int = 1,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def anndata_or_mudata_checker(
    function: Optional[Callable[P, R]] = None,
    n: int = 1,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
    """
    Decorate a function by checking that the first arguments are
    AnnData or MuData objects.

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

    """

    if function is not None:

        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            iterator = iter(list(args) + list(kwargs.values()))
            for i in range(n):
                value = next(iterator)
                if _mudata_is_available:
                    valid = isinstance(value, (AnnData, MuData))
                    expected = f"{AnnData} or {MuData}"
                else:
                    valid = isinstance(value, AnnData)
                    expected = str(AnnData)

                if not valid:
                    raise TypeError(
                        f"unsupported argument type for '{i+1}'-th argument: "
                        f"expected {expected} but received {type(value)}"
                    )
            return function(*args, **kwargs)

        return wrapper

    def partial_wrapper(function: Callable[P, R]) -> Callable[P, R]:
        return anndata_or_mudata_checker(function, n)

    return partial_wrapper
