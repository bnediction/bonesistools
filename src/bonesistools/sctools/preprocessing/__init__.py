#!/usr/bin/env python

"""
Preprocessing utilities for AnnData objects.

The `pp` namespace contains helpers for filtering observations and variables,
regressing out covariates, merging annotations and transferring data between
related AnnData objects.
"""

from ._simple import (
    filter_obs,
    filter_var,
    regress_out,
    merge,
    transfer_layer,
    transfer_obs_sti,
    transfer_obs_its,
)

from ._genename import (
    convert_gene_identifiers,
    standardize_gene_identifiers,
    var_names_merge_duplicates,
)

__all__ = [
    "filter_obs",
    "filter_var",
    "regress_out",
    "merge",
    "transfer_layer",
    "transfer_obs_sti",
    "transfer_obs_its",
    "convert_gene_identifiers",
    "standardize_gene_identifiers",
    "var_names_merge_duplicates",
]


def __dir__():
    return sorted(set(globals()) | set(__all__))
