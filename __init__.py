#!/usr/bin/env python

from . import anndatatools as adt
from . import boolpy as bpy
from . import databases as dbs
from . import grntools as grn
from . import utils

import sys as _sys

_sys.modules.update({f"{__name__}.{alias}": globals()[alias] for alias in ["adt", "bpy", "dbs", "grn"]})
