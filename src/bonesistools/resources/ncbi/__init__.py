#!/usr/bin/env python

"""
Interfaces to NCBI Gene information resources.

The `ncbi` sub-package provides organism-specific gene metadata and tools for
resolving and converting gene identifiers.
"""

from typing import List as _List

from ._genesyn import GeneSynonyms as GeneSynonyms
from ._genesyn import genesyn as genesyn
from ._identifiers import GeneIdentifiers as GeneIdentifiers
from ._identifiers import identifiers

__all__ = [
    "GeneIdentifiers",
    "identifiers",
]


def __dir__() -> _List[str]:
    hidden = {"GeneSynonyms", "genesyn"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
