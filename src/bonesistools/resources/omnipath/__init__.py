#!/usr/bin/env python

"""
Interfaces to OmniPath-derived prior knowledge resources.

The `omnipath` sub-package provides utilities for retrieving and
constructing biological interaction priors and regulatory domains
from OmniPath-related resources such as CollecTRI and DoRothEA.
"""

from typing import List as _List

from ._collectri import collectri
from ._collectri import load_collectri_grn as load_collectri_grn
from ._dorothea import dorothea
from ._dorothea import load_dorothea_grn as load_dorothea_grn

__all__ = [
    "collectri",
    "dorothea",
]


def __dir__() -> _List[str]:
    hidden = {"load_collectri_grn", "load_dorothea_grn"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
