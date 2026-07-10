#!/usr/bin/env python

"""
Utilities for Boolean modelling.

The `logic` package provides tools for Boolean algebra, Boolean networks,
and influence graphs.

Sub-packages
------------
ba
    Boolean algebra and partial Boolean abstractions.
bn
    Boolean networks and logical modelling utilities.
ig
    Influence graph utilities.
io
    Input/output helpers.
"""

from __future__ import annotations

import sys as _sys
from types import ModuleType as _ModuleType
from typing import Dict as _Dict
from typing import List as _List
from typing import Tuple as _Tuple
from typing import cast as _cast

from .._warnings import _warn_deprecated
from . import boolean_algebra as ba
from . import boolean_network as bn
from . import influence_graph as ig
from . import input_output as io

for _name in [
    "boolean_algebra",
    "boolean_network",
    "influence_graph",
    "input_output",
]:
    globals().pop(_name, None)

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["ba", "bn", "ig", "io"]}
)

_MODULES: _Dict[str, _ModuleType] = {
    "ba": ba,
    "bn": bn,
    "ig": ig,
    "io": io,
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
    "io",
]


def __getattr__(name: str) -> object:
    if name in _DEPRECATED:
        module_alias, attr = _DEPRECATED[name]
        _warn_deprecated(
            f"`bt.logic.{name}`",
            replacement=f"`bt.logic.{module_alias}.{attr}`",
            stacklevel=2,
        )
        return _cast(object, getattr(_MODULES[module_alias], attr))

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> _List[str]:
    return sorted(name for name in set(globals()) | set(__all__) if name[0] != "_")
