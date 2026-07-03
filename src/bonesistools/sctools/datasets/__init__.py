#!/usr/bin/env python

"""
Registered example datasets for single-cell workflows.

The `datasets` namespace exposes AnnData-oriented dataset loaders and
lightweight objects used in examples, tests and tutorials. Prefer
`datasets.load(...)` for registered datasets.
"""

from __future__ import annotations

import warnings
from typing import List as _List

from anndata import AnnData

from ._geo import from_geo
from ._registry import available, clear, info, load

del annotations

__all__ = [
    "available",
    "clear",
    "from_geo",
    "info",
    "load",
    "nestorowa",
]


def nestorowa(
    quiet: bool = False,
) -> AnnData:
    """
    Load the Nestorowa dataset.

    Deprecated. Use `datasets.load("nestorowa")` instead.
    """

    warnings.warn(
        "`datasets.nestorowa()` is deprecated and will be removed in a future "
        'release. Use `datasets.load("nestorowa")` instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    return load("nestorowa", quiet=quiet)


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
