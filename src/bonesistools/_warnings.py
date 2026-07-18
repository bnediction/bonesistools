#!/usr/bin/env python

import warnings
from functools import wraps
from typing import Any, Callable, Optional, TypeVar, cast

_F = TypeVar("_F", bound=Callable[..., Any])


def _warn_deprecated(
    name: str,
    *,
    replacement: Optional[str] = None,
    stacklevel: int = 2,
) -> None:

    message = f"{name} is deprecated and will be removed in 2.0.0"
    if replacement is not None:
        message = f"{message}; use {replacement} instead"
    warnings.warn(f"{message}.", FutureWarning, stacklevel=stacklevel)


def _warn_deprecated_argument(
    old_name: str,
    new_name: str,
    *,
    stacklevel: int = 2,
) -> None:

    _warn_deprecated(
        f"`{old_name}`",
        replacement=f"`{new_name}`",
        stacklevel=stacklevel,
    )


def _deprecated(*, replacement: Optional[str] = None) -> Callable[[_F], _F]:
    """Mark a callable as deprecated while preserving its public signature."""

    def decorator(function: _F) -> _F:
        @wraps(function)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            _warn_deprecated(
                f"`{function.__qualname__}()`",
                replacement=replacement,
                stacklevel=3,
            )
            return function(*args, **kwargs)

        return cast(_F, wrapped)

    return decorator


def _rename_deprecated_arguments(**replacements: str) -> Callable[[_F], _F]:
    """Accept deprecated keyword names while exposing the current signature."""

    def decorator(function: _F) -> _F:
        @wraps(function)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            for old_name, new_name in replacements.items():
                if old_name not in kwargs:
                    continue

                _warn_deprecated_argument(old_name, new_name, stacklevel=4)
                if new_name in kwargs:
                    raise TypeError(
                        "invalid argument combination: use either "
                        f"'{old_name}' or '{new_name}', not both"
                    )
                kwargs[new_name] = kwargs.pop(old_name)

            return function(*args, **kwargs)

        return cast(_F, wrapped)

    return decorator
