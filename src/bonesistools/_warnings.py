#!/usr/bin/env python

import warnings
from typing import Optional


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
