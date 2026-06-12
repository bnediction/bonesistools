#!/usr/bin/env python

"""
Utilities for Boolean modelling.

The `bpy` package provides tools for Boolean algebra, Boolean networks,
influence graphs and plotting.

Sub-packages
------------
ba
    Boolean algebra and partial Boolean abstractions.
bn
    Boolean networks and logical modelling utilities.
ig
    Influence graph utilities.
pl
    Plotting utilities.
"""

from __future__ import annotations

import sys as _sys
import warnings as _warnings
from types import ModuleType as _ModuleType
from typing import Dict as _Dict
from typing import List as _List
from typing import Tuple as _Tuple
from typing import cast as _cast

from . import boolean_algebra as ba
from . import boolean_network as bn
from . import influence_graph as ig
from . import plotting as pl

del annotations

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["ba", "bn", "ig", "pl"]}
)

_MODULES: _Dict[str, _ModuleType] = {
    "ba": ba,
    "bn": bn,
    "ig": ig,
    "pl": pl,
}

_DEPRECATED: _Dict[str, _Tuple[str, str]] = {
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
    "pl",
]


def __getattr__(name: str) -> object:
    if name in _DEPRECATED:
        module_alias, attr = _DEPRECATED[name]
        _warnings.warn(
            f"`bt.bpy.{name}` is deprecated and will be removed in 2.0.0; use "
            f"`bt.bpy.{module_alias}.{attr}` instead.",
            FutureWarning,
            stacklevel=2,
        )
        return _cast(object, getattr(_MODULES[module_alias], attr))

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
