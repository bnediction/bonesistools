#!/usr/bin/env python

"""
Utilities for NCBI-based gene nomenclature and synonym handling.

The `ncbi` sub-package provides tools for resolving gene aliases,
standardising gene identifiers and handling ambiguous nomenclature
across heterogeneous biological resources.
"""

from typing import List as _List

from ._genesyn import GeneSynonyms

__all__ = [
    "GeneSynonyms",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
