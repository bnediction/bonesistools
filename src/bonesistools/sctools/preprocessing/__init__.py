#!/usr/bin/env python

"""
Preprocessing utilities for AnnData objects.

The `pp` namespace contains helpers for filtering observations and variables,
regressing out covariates, merging annotations and transferring data between
related AnnData objects.
"""

import warnings as _warnings
from typing import Any as _Any
from typing import List as _List

from ..tools import regress_out as _regress_out
from ._genename import (
    convert_gene_identifiers,
    standardize_gene_identifiers,
    var_names_merge_duplicates,
)
from ._simple import (
    filter_obs,
    filter_var,
    merge,
)
from ._transfer import (
    transfer_layer,
    transfer_obs_its,
    transfer_obs_sti,
)

__all__ = [
    "filter_obs",
    "filter_var",
    "merge",
    "transfer_layer",
    "transfer_obs_sti",
    "transfer_obs_its",
    "convert_gene_identifiers",
    "standardize_gene_identifiers",
    "var_names_merge_duplicates",
    "regress_out",
]


def regress_out(*args: _Any, **kwargs: _Any) -> _Any:
    """
    Deprecated. Regress out unwanted sources of variation.

    Use `bonesistools.sctools.tools.regress_out` or `bt.sct.tl.regress_out`
    instead.
    """

    _warnings.warn(
        "`bt.sct.pp.regress_out` is deprecated; use "
        "`bt.sct.tl.regress_out` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _regress_out(*args, **kwargs)


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
