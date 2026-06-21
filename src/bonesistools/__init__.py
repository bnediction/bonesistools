#!/usr/bin/env python

"""
BoNesisTools provides bioinformatics utilities for upstream and downstream
analyses of the BoNesis framework.

Packages
--------
sct
    Single-cell and multimodal annotated data tools.
bpy
    Boolean modelling utilities.
dbs
    Biological database interfaces.
grn
    Deprecated alias for `bpy.ig`.

Credits: BNeDiction; PEPR Santé Numérique 2030.
"""

from __future__ import annotations

import importlib as _importlib
import sys as _sys
from types import ModuleType as _ModuleType
from typing import TYPE_CHECKING as _TYPE_CHECKING
from typing import List as _List

from . import boolpy as bpy
from . import databases as dbs
from . import sctools as sct
from ._metadata import package_version as _package_version
from ._warnings import _warn_deprecated

if _TYPE_CHECKING:
    from . import grntools as grn

del annotations

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"
__version__ = _package_version()

__all__ = [
    "__version__",
    "sct",
    "bpy",
    "dbs",
    "grn",
]

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["sct", "bpy", "dbs"]}
)

_sys.modules.update(
    {
        f"{__name__}.sct.{alias}": getattr(sct, alias)
        for alias in ["pp", "tl", "pl", "datasets"]
    }
)

_sys.modules.update(
    {
        f"{__name__}.bpy.{alias}": getattr(bpy, alias)
        for alias in ["ba", "bn", "ig", "pl"]
    }
)


def __getattr__(name: str) -> _ModuleType:
    if name == "grn":
        _warn_deprecated(
            "`bt.grn`",
            replacement="`bt.bpy.ig`",
            stacklevel=2,
        )

        grn = _importlib.import_module(".grntools", __name__)
        _sys.modules[f"{__name__}.grn"] = grn

        return grn

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
