#!/usr/bin/env python

"""
Interfaces to OmniPath-derived prior knowledge resources.

The `omnipath` sub-package provides utilities for retrieving and
constructing biological interaction priors and regulatory domains
from OmniPath-related resources such as CollecTRI and DoRothEA.
"""

from ._collectri import load_collectri_grn
from ._dorothea import load_dorothea_grn
