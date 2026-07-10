#!/usr/bin/env python

"""
Utilities for single-cell and multimodal annotated data.

The `omics` package provides preprocessing, analysis and plotting tools
for AnnData and MuData objects, following a Scanpy-like API while
offering additional and complementary features.

Sub-packages
------------
pp
    Preprocessing utilities.
tl
    Analysis and inference tools.
io
    Input/output helpers.
pl
    Plotting utilities.
"""

import sys as _sys
from importlib import import_module as _import_module
from typing import List as _List

pp = _import_module(f"{__name__}.preprocessing")
tl = _import_module(f"{__name__}.tools")
io = _import_module(f"{__name__}.input_output")
pl = _import_module(f"{__name__}.plotting")
_datasets = _import_module(f"{__name__}.datasets")

for _name in ["preprocessing", "tools", "input_output", "plotting", "datasets"]:
    globals().pop(_name, None)

datasets = _datasets

__all__ = [
    "pp",
    "tl",
    "io",
    "pl",
]

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["pp", "tl", "io", "pl"]}
)
_sys.modules[f"{__name__}.datasets"] = _datasets


def __dir__() -> _List[str]:
    hidden = {"datasets", "input_output", "plotting", "preprocessing", "tools"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
