#!/usr/bin/env python

"""
Utilities for NCBI-based gene nomenclature and synonym handling.

The `ncbi` sub-package provides tools for resolving gene aliases,
standardising gene identifiers and handling ambiguous nomenclature
across heterogeneous biological resources.
"""

from ._genesyn import (
    InputIdentifierType,
    OutputIdentifierType,
    GeneSynonyms,
    support_legacy_gene_synonyms_args,
)
