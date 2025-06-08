#!/usr/bin/env python

"""
Bonesistools is a python package proposing functionnalities and analysis upstream and downstream of BoNesis framework.

credits: "BNeDiction; PEPR Santé Numérique 2030"
"""

__credits__ = "BNeDiction; PEPR Santé Numérique 2030"

from . import sctools as sct
from . import boolpy as bpy
from . import databases as dbs
from . import grntools as grn

import sys as _sys

_sys.modules.update({f"{__name__}.{alias}": globals()[alias] for alias in ["sct", "bpy", "dbs", "grn"]})
