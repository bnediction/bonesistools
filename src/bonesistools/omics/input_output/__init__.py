#!/usr/bin/env python

"""
Input/output helpers for single-cell workflows.

The `io` namespace centralizes AnnData dataset loading, GEO import and matrix
export helpers.
"""

from typing import List as _List

from ._geo import from_geo
from ._registry import available, clear, info, load
from ._write import to_csv, to_mtx, to_npz

__all__ = [
    "available",
    "clear",
    "from_geo",
    "info",
    "load",
    "to_csv",
    "to_mtx",
    "to_npz",
]


def nestorowa(quiet: bool = False):
    """
    Load the Nestorowa dataset.

    Deprecated. Use `bt.omics.io.load("nestorowa")` instead.
    """

    from ..._warnings import _warn_deprecated

    _warn_deprecated(
        "`bt.omics.io.nestorowa()`",
        replacement='`bt.omics.io.load("nestorowa")`',
        stacklevel=2,
    )

    return load("nestorowa", quiet=quiet)


def __dir__() -> _List[str]:
    hidden = {"nestorowa"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
