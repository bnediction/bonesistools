#!/usr/bin/env python

import importlib
from functools import wraps
from typing import Callable, Optional, TypeVar, Union, overload

from typing_extensions import ParamSpec

P = ParamSpec("P")
R = TypeVar("R")


@overload
def require_dependency(
    function: Callable[P, R],
    *,
    module: str,
    package: str,
    extra: str,
) -> Callable[P, R]: ...


@overload
def require_dependency(
    function: None = None,
    *,
    module: str,
    package: str,
    extra: str,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def require_dependency(
    function: Optional[Callable[P, R]] = None,
    *,
    module: str,
    package: str,
    extra: str,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
    """
    Decorate a function with an optional dependency check.

    The wrapped function imports `module` lazily when called. If the import
    fails, an ImportError explains which package and optional extra should be
    installed.

    Parameters
    ----------
    function: Callable, optional
        Function to decorate. If None, return a decorator.
    module: str
        Python module imported to check dependency availability.
    package: str
        Package name displayed in the ImportError message.
    extra: str
        bonesistools optional extra displayed in the ImportError message.

    Returns
    -------
    Callable
        Decorated function if `function` is provided, otherwise a decorator.
    """

    def decorator(function: Callable[P, R]) -> Callable[P, R]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            try:
                importlib.import_module(module)
            except ImportError as error:
                raise ImportError(
                    f"{package} is required for this function. "
                    f"Install bonesistools with the {extra} extra or install {package}."
                ) from error
            return function(*args, **kwargs)

        return wrapper

    return decorator(function) if function is not None else decorator


@overload
def require_sklearn(function: Callable[P, R]) -> Callable[P, R]: ...


@overload
def require_sklearn(
    function: None = None,
) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def require_sklearn(
    function: Optional[Callable[P, R]] = None,
) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
    """
    Decorate a function requiring scikit-learn.

    Parameters
    ----------
    function: Callable, optional
        Function to decorate. If None, return a decorator.

    Returns
    -------
    Callable
        Decorated function if `function` is provided, otherwise a decorator.
    """

    return require_dependency(
        function,
        module="sklearn",
        package="scikit-learn",
        extra="sctools",
    )
