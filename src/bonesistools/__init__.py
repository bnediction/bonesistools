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

Credits: BNeDiction; PEPR Santé Numérique 2030.
"""

from __future__ import annotations

import sys as _sys
from importlib import import_module as _import_module
from typing import List as _List

from ._metadata import package_version as _package_version

bpy = _import_module(f"{__name__}.boolpy")
dbs = _import_module(f"{__name__}.databases")
sct = _import_module(f"{__name__}.sctools")

boolpy = bpy
databases = dbs
sctools = sct

del annotations

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"
__version__ = _package_version()

__all__ = [
    "__version__",
    "sct",
    "bpy",
    "dbs",
]

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["sct", "bpy", "dbs"]}
)

_sys.modules.update(
    {
        f"{__name__}.sct.{alias}": getattr(sct, alias)
        for alias in ["pp", "tl", "io", "pl"]
    }
)
_sys.modules[f"{__name__}.sct.datasets"] = sct._datasets

_sys.modules.update(
    {
        f"{__name__}.bpy.{alias}": getattr(bpy, alias)
        for alias in ["ba", "bn", "ig", "io"]
    }
)


def __dir__() -> _List[str]:
    hidden = {"boolpy", "databases", "sctools"}
    return sorted((set(globals()) | set(__all__)) - hidden)
