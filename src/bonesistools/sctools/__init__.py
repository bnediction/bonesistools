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

from . import _typing as typing

from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets

import sys as _sys

_sys.modules.update(
    {f"{__name__}.{alias}": globals()[alias] for alias in ["tl", "pp", "pl"]}
)
