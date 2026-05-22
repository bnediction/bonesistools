#!/usr/bin/env python

import importlib
from functools import wraps


def require_dependency(function=None, *, module: str, package: str, extra: str):
    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
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


def require_sklearn(function=None):
    return require_dependency(
        function,
        module="sklearn",
        package="scikit-learn",
        extra="sctools",
    )
