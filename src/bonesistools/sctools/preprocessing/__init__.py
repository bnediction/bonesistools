#!/usr/bin/env python

"""
Preprocessing utilities for AnnData objects.

The `pp` namespace contains helpers for filtering observations and variables,
annotating features, regressing out covariates, merging annotations and
transferring data between related AnnData objects.
"""

from typing import Any as _Any
from typing import List as _List

from ..._warnings import _warn_deprecated
from ..tools import regress_out as _regress_out
from ._classification import (
    mitochondrial_genes,
    ribosomal_genes,
)
from ._duplicates import (
    merge_duplicate_vars,
)
from ._duplicates import var_names_merge_duplicates as var_names_merge_duplicates
from ._filter import (
    filter_obs,
    filter_var,
)
from ._genename import (
    convert_gene_identifiers,
)
from ._genename import standardize_gene_identifiers as standardize_gene_identifiers
from ._hvg import hvg
from ._qc import qc
from ._simple import (
    sort,
)
from ._simple import sort_anndata as sort_anndata
from ._transfer import (
    merge,
    transfer_layer,
    transfer_obs_to_integrated,
    transfer_obs_to_specific,
)
from ._transfer import transfer_obs_its as transfer_obs_its
from ._transfer import transfer_obs_sti as transfer_obs_sti
from ._transform import (
    log1p,
    normalize,
    scale,
)

__all__ = [
    "filter_obs",
    "filter_var",
    "normalize",
    "log1p",
    "scale",
    "qc",
    "mitochondrial_genes",
    "ribosomal_genes",
    "hvg",
    "merge",
    "sort",
    "transfer_layer",
    "transfer_obs_to_integrated",
    "transfer_obs_to_specific",
    "convert_gene_identifiers",
    "merge_duplicate_vars",
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
    hidden = {
        "regress_out",
        "sort_anndata",
        "standardize_gene_identifiers",
        "transfer_obs_its",
        "transfer_obs_sti",
        "var_names_merge_duplicates",
    }
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
