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

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"

import importlib as _importlib
import sys as _sys
import warnings as _warnings

from . import sctools as sct
from . import boolpy as bpy
from . import databases as dbs

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["sct", "bpy", "dbs"]}
)

def __getattr__(name):
    if name == "grn":
        _warnings.warn(
            "`bt.grn` is deprecated and will be removed in a future release; "
            "use `bt.bpy.ig` instead.",
            FutureWarning,
            stacklevel=2,
        )
        grn = _importlib.import_module(".grntools", __name__)
        _sys.modules[f"{__name__}.grn"] = grn
        return grn

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "pp",
    "tl",
    "pl",
    "datasets",
    "typing",
]