#!/usr/bin/env python

"""
Interfaces to biological databases and prior knowledge resources.

The `resources` package provides utilities for retrieving, parsing and
standardising biological information from external resources.

Current sub-packages
--------------------
ncbi
    Gene nomenclature and synonym handling utilities.
hcop
    HCOP orthology translation utilities.
omnipath
    Access to regulatory and interaction prior knowledge databases.
"""

from typing import List as _List

from . import hcop, ncbi, omnipath

__all__ = [
    "hcop",
    "ncbi",
    "omnipath",
]


def __dir__() -> _List[str]:
    return sorted(name for name in set(globals()) | set(__all__) if name[0] != "_")
