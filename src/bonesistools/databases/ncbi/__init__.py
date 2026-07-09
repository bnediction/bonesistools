#!/usr/bin/env python

"""
Utilities for NCBI-based gene nomenclature and synonym handling.

The `ncbi` sub-package provides tools for resolving gene aliases,
standardising gene identifiers and handling ambiguous nomenclature
across heterogeneous biological resources.
"""

from typing import List as _List

from ._genesyn import GeneSynonyms as GeneSynonyms
from ._genesyn import genesyn

__all__ = [
    "genesyn",
]


def __dir__() -> _List[str]:
    hidden = {"GeneSynonyms"}
    return sorted((set(globals()) | set(__all__)) - hidden)
