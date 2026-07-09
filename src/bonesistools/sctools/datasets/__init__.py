#!/usr/bin/env python

"""
Deprecated dataset namespace for single-cell workflows.

Use `bt.sct.io` for dataset loading, GEO import and cache management.
"""

from __future__ import annotations

from typing import List as _List

from anndata import AnnData as _AnnData

from ..._warnings import _warn_deprecated
from ._geo import from_geo as _from_geo
from ._registry import available as _available
from ._registry import clear as _clear
from ._registry import info as _info
from ._registry import load as _load

del annotations

__all__ = [
    "available",
    "clear",
    "from_geo",
    "info",
    "load",
    "nestorowa",
]


def load(
    name: str,
    quiet: bool = False,
) -> _AnnData:
    """
    Load a built-in single-cell dataset.

    Deprecated. Use `bt.sct.io.load` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.load()`",
        replacement="`bt.sct.io.load()`",
        stacklevel=2,
    )

    return _load(name, quiet=quiet)


def info(name: str):
    """
    Return metadata for a registered dataset.

    Deprecated. Use `bt.sct.io.info` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.info()`",
        replacement="`bt.sct.io.info()`",
        stacklevel=2,
    )

    return _info(name)


def available():
    """
    List registered datasets.

    Deprecated. Use `bt.sct.io.available` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.available()`",
        replacement="`bt.sct.io.available()`",
        stacklevel=2,
    )

    return _available()


def clear(*names: str) -> None:
    """
    Remove cached dataset files.

    Deprecated. Use `bt.sct.io.clear` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.clear()`",
        replacement="`bt.sct.io.clear()`",
        stacklevel=2,
    )

    _clear(*names)


def from_geo(
    accession: str,
    cache_dir=None,
    quiet: bool = False,
) -> _AnnData:
    """
    Download a GEO dataset and return it as an AnnData object.

    Deprecated. Use `bt.sct.io.from_geo` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.from_geo()`",
        replacement="`bt.sct.io.from_geo()`",
        stacklevel=2,
    )

    return _from_geo(accession, cache_dir=cache_dir, quiet=quiet)


def nestorowa(
    quiet: bool = False,
) -> _AnnData:
    """
    Load the Nestorowa dataset.

    Deprecated. Use `bt.sct.io.load("nestorowa")` instead.
    """

    _warn_deprecated(
        "`bt.sct.datasets.nestorowa()`",
        replacement='`bt.sct.io.load("nestorowa")`',
        stacklevel=2,
    )

    return _load("nestorowa", quiet=quiet)


def __dir__() -> _List[str]:
    return sorted(name for name in set(globals()) - set(__all__) if name[0] != "_")
