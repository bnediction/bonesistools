#!/usr/bin/env python

"""
sctools handles and performs operations on unimodal and multimodal annotated data matrix (anndata/mudata). \
It offers a range of efficient features not available in Scanpy package, proposes some corrections \
on Scanpy-based functions when inconsistencies have been highlighted, and matplotlib-based visualizations.
"""

from ._typing import (
    UnionType,
    type_checker,
    anndata_checker,
    mudata_checker,
    anndata_or_mudata_checker,
    AnnDataList,
    MuDataList,
    DataFrameList,
    AxisInt,
    Axis,
    Keys,
    Suffixes,
    ScData,
    Metric
)

from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl

import sys

sys.modules.update({f"{__name__}.{alias}": globals()[alias] for alias in ["tl", "pp", "pl"]})
