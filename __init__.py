#!/usr/bin/env python

from . import anndatatools as at
from . import boolpy as bp
from . import databases as db
from . import grntools as gt
from . import utils

import sys as _sys

_sys.modules.update({f"{__name__}.{alias}": globals()[alias] for alias in ["at", "bp", "db", "gt"]})
