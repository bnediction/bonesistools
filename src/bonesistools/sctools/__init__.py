#!/usr/bin/env python

"""
Utilities for single-cell and multimodal annotated data.

The `sct` package provides preprocessing, analysis and plotting tools
for AnnData and MuData objects, following a Scanpy-like API while
offering additional and complementary features.

Sub-packages
------------
pp
    Preprocessing utilities.
tl
    Analysis and inference tools.
pl
    Plotting utilities.
datasets
    Bundled example datasets.
"""

import sys as _sys
from typing import List

from . import _typing as typing
from . import datasets
from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl

__all__ = [
    "pp",
    "tl",
    "pl",
    "datasets",
    "typing",
]

_sys.modules.update(
    {
        f"{__name__}.{alias}": globals()[alias]
        for alias in ["pp", "tl", "pl", "datasets"]
    }
)


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
