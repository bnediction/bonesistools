#!/usr/bin/env python

"""
Interfaces to HCOP orthology mappings.

The `hcop` sub-package provides utilities for translating gene symbols between
organisms supported by the HGNC Comparison of Orthology Predictions resource.
"""

from typing import List as _List

from ._orthologs import organisms, orthologs, versions

__all__ = [
    "organisms",
    "orthologs",
    "versions",
]


def __dir__() -> _List[str]:
    return sorted(name for name in set(globals()) | set(__all__) if name[0] != "_")
