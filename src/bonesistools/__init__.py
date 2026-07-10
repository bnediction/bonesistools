#!/usr/bin/env python

"""
BoNesisTools provides bioinformatics utilities for upstream and downstream
analyses of the BoNesis framework.

Packages
--------
omics
    Single-cell and multimodal annotated data tools.
logic
    Boolean modelling utilities.
resources
    Biological database and prior-knowledge resources.

Credits: BNeDiction; PEPR Santé Numérique 2030.
"""

from __future__ import annotations

import sys as _sys
from types import ModuleType as _ModuleType
from typing import Dict as _Dict
from typing import List as _List
from typing import Tuple as _Tuple

from . import logic, omics, resources
from ._metadata import package_version as _package_version
from ._warnings import _warn_deprecated

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"
__version__ = _package_version()

__all__ = [
    "__version__",
    "omics",
    "logic",
    "resources",
]

_DEPRECATED_ALIASES: _Dict[str, _Tuple[str, _ModuleType]] = {
    "sct": ("omics", omics),
    "bpy": ("logic", logic),
    "dbs": ("resources", resources),
}

_sys.modules.update(
    {
        f"{__name__}.{alias}": module
        for alias, (_, module) in _DEPRECATED_ALIASES.items()
    }
)
_sys.modules[f"{__name__}.sctools"] = omics
_sys.modules[f"{__name__}.boolpy"] = logic
_sys.modules[f"{__name__}.databases"] = resources

_sys.modules.update(
    {
        f"{__name__}.omics.{alias}": getattr(omics, alias)
        for alias in ["pp", "tl", "io", "pl"]
    }
)
_sys.modules[f"{__name__}.omics.datasets"] = getattr(omics, "_datasets")
_sys.modules.update(
    {
        f"{__name__}.sctools.{alias}": getattr(omics, alias)
        for alias in ["pp", "tl", "io", "pl"]
    }
)
_sys.modules[f"{__name__}.sctools.datasets"] = getattr(omics, "_datasets")
_sys.modules.update(
    {
        f"{__name__}.sct.{alias}": getattr(omics, alias)
        for alias in ["pp", "tl", "io", "pl"]
    }
)
_sys.modules[f"{__name__}.sct.datasets"] = getattr(omics, "_datasets")

_sys.modules.update(
    {
        f"{__name__}.logic.{alias}": getattr(logic, alias)
        for alias in ["ba", "bn", "ig", "io"]
    }
)
_sys.modules.update(
    {
        f"{__name__}.boolpy.{alias}": getattr(logic, alias)
        for alias in ["ba", "bn", "ig", "io"]
    }
)

_sys.modules.update(
    {
        f"{__name__}.resources.{alias}": getattr(resources, alias)
        for alias in ["hcop", "ncbi", "omnipath"]
    }
)
_sys.modules.update(
    {
        f"{__name__}.dbs.{alias}": getattr(resources, alias)
        for alias in ["hcop", "ncbi", "omnipath"]
    }
)
_sys.modules.update(
    {
        f"{__name__}.databases.{alias}": getattr(resources, alias)
        for alias in ["hcop", "ncbi", "omnipath"]
    }
)
_sys.modules.update(
    {
        f"{__name__}.bpy.{alias}": getattr(logic, alias)
        for alias in ["ba", "bn", "ig", "io"]
    }
)


def __getattr__(name: str) -> object:
    if name in _DEPRECATED_ALIASES:
        replacement, module = _DEPRECATED_ALIASES[name]
        _warn_deprecated(
            f"`bt.{name}`",
            replacement=f"`bt.{replacement}`",
            stacklevel=2,
        )
        return module

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> _List[str]:
    return sorted(name for name in (set(globals()) | set(__all__)) if name[0] != "_")
