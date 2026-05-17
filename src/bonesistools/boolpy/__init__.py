#!/usr/bin/env python

"""
Utilities for Boolean modelling.

The `bpy` package provides tools for Boolean algebra, Boolean networks
and interaction/influence graphs.

Sub-packages
------------
ba
    Boolean algebra and partial Boolean abstractions.
bn
    Boolean networks and logical modelling utilities.
ig
    Interaction and influence graph utilities.
"""

import warnings as _warnings
import sys as _sys

from . import boolean_algebra as ba
from . import boolean_network as bn
from . import interaction_graph as ig

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["ba", "bn", "ig"]}
)

_DEPRECATED = {
    "BooleanDifferentialCalculus": ("ba", "BooleanDifferentialCalculus"),
    "PartialBoolean": ("ba", "PartialBoolean"),
    "BooleanNetworkEnsemble": ("bn", "BooleanNetworkEnsemble"),
    "Hypercube": ("bn", "Hypercube"),
    "HypercubeCollection": ("bn", "HypercubeCollection"),
    "bn_to_pydot": ("bn", "bn_to_pydot"),
}

def __getattr__(name):
    if name in _DEPRECATED:
        module_alias, attr = _DEPRECATED[name]
        _warnings.warn(
            f"`bt.bpy.{name}` is deprecated; use `bt.bpy.{module_alias}.{attr}` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(globals()[module_alias], attr)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "ba",
    "bn",
    "ig",
]