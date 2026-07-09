#!/usr/bin/env python

"""
Input/output helpers for single-cell workflows.

The `io` namespace centralizes AnnData dataset loading, GEO import and matrix
export helpers.
"""

from typing import List as _List

from ..datasets._geo import from_geo
from ..datasets._registry import available, clear, info, load
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

    Deprecated. Use `bt.sct.io.load("nestorowa")` instead.
    """

    from ..._warnings import _warn_deprecated

    _warn_deprecated(
        "`bt.sct.io.nestorowa()`",
        replacement='`bt.sct.io.load("nestorowa")`',
        stacklevel=2,
    )

    return load("nestorowa", quiet=quiet)


def __dir__() -> _List[str]:
    hidden = {"nestorowa"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
