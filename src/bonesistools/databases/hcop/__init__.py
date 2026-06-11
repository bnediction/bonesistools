#!/usr/bin/env python

"""
Interfaces to HCOP orthology mappings.

The `hcop` sub-package provides utilities for translating human gene symbols
to supported target organisms with the HGNC Comparison of Orthology
Predictions resource.
"""

from typing import List as _List

from ._orthologs import organisms, orthologs

__all__ = [
    "organisms",
    "orthologs",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
