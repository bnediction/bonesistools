#!/usr/bin/env python

"""
Preprocessing utilities for AnnData objects.

The `pp` namespace contains helpers for filtering observations and variables,
regressing out covariates, merging annotations and transferring data between
related AnnData objects.
"""

from typing import Any as _Any
from typing import List as _List

from ..._warnings import _warn_deprecated
from ..tools import regress_out as _regress_out
from ._duplicates import (
    merge_duplicate_vars,
    var_names_merge_duplicates,
)
from ._genename import (
    convert_gene_identifiers,
    standardize_gene_identifiers,
)
from ._simple import (
    filter_obs,
    filter_var,
    sort_anndata,
)
from ._transfer import (
    merge,
    transfer_layer,
    transfer_obs_its,
    transfer_obs_sti,
    transfer_obs_to_integrated,
    transfer_obs_to_specific,
)

__all__ = [
    "filter_obs",
    "filter_var",
    "merge",
    "sort_anndata",
    "transfer_layer",
    "transfer_obs_to_integrated",
    "transfer_obs_to_specific",
    "transfer_obs_sti",
    "transfer_obs_its",
    "convert_gene_identifiers",
    "standardize_gene_identifiers",
    "merge_duplicate_vars",
    "var_names_merge_duplicates",
    "regress_out",
]


def regress_out(*args: _Any, **kwargs: _Any) -> _Any:
    """
    Deprecated. Regress out unwanted sources of variation.

    Use `bonesistools.sctools.tools.regress_out` or `bt.sct.tl.regress_out`
    instead.
    """

    _warn_deprecated(
        "`bt.sct.pp.regress_out`",
        replacement="`bt.sct.tl.regress_out`",
        stacklevel=2,
    )
    return _regress_out(*args, **kwargs)


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
