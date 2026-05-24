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

from __future__ import annotations

import sys as _sys
import warnings as _warnings
from types import ModuleType
from typing import cast

from . import boolean_algebra as ba
from . import boolean_network as bn
from . import interaction_graph as ig

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["ba", "bn", "ig"]}
)

_MODULES: dict[str, ModuleType] = {
    "ba": ba,
    "bn": bn,
    "ig": ig,
}

_DEPRECATED: dict[str, tuple[str, str]] = {
    "PartialBoolean": ("ba", "PartialBoolean"),
    "BooleanNetworkEnsemble": ("bn", "BooleanNetworkEnsemble"),
    "Hypercube": ("ba", "Hypercube"),
    "HypercubeCollection": ("ba", "HypercubeCollection"),
    "bn_to_pydot": ("bn", "bn_to_pydot"),
}

__all__ = [
    "ba",
    "bn",
    "ig",
]


def __getattr__(name: str) -> object:
    if name in _DEPRECATED:
        module_alias, attr = _DEPRECATED[name]
        _warnings.warn(
            f"`bt.bpy.{name}` is deprecated; use `bt.bpy.{module_alias}.{attr}` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return cast(object, getattr(_MODULES[module_alias], attr))

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
