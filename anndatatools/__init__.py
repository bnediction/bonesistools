#!/usr/bin/env python

"""
anndatatools is a python package for handling and performing operations on annotated data (anndata).
It offers a range of efficient features not available in Scanpy package, proposes some corrections
on Scanpy-based functions when inconsistencies have been highlighted, and matplotlib-based visualizations.

credits: "BNeDiction; PEPR Santé Numérique 2030"
"""

__version__ = "1.0.0"
__credits__ = "BNeDiction; PEPR Santé Numérique 2030"

from ._check_anndata import _adata_arg_checking

from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl

import sys

sys.modules.update({f"{__name__}.{alias}": globals()[alias] for alias in ["tl", "pp", "pl"]})
