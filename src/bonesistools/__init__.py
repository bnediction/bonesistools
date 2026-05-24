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
import warnings as _warnings
from types import ModuleType
from typing import TYPE_CHECKING

from . import sctools as sct
from . import boolpy as bpy
from . import databases as dbs

if TYPE_CHECKING:
    from . import grntools as grn

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"

__all__ = [
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
    {f"{__name__}.bpy.{alias}": getattr(bpy, alias) for alias in ["ba", "bn", "ig"]}
)


def __getattr__(name: str) -> ModuleType:
    if name == "grn":
        message = (
            "`bt.grn` is deprecated and will be removed in a future release; "
            "use `bt.bpy.ig` instead."
        )
        _warnings.warn(
            message,
            FutureWarning,
            stacklevel=2,
        )

        grn = _importlib.import_module(".grntools", __name__)
        _sys.modules[f"{__name__}.grn"] = grn

        return grn

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
