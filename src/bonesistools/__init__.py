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
from typing import List as _List

from . import boolpy as bpy
from . import databases as dbs
from . import sctools as sct
from ._metadata import package_version as _package_version

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
        for alias in ["pp", "tl", "pl", "datasets"]
    }
)

_sys.modules.update(
    {
        f"{__name__}.bpy.{alias}": getattr(bpy, alias)
        for alias in ["ba", "bn", "ig", "pl"]
    }
)


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
