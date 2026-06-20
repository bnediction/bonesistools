#!/usr/bin/env python

from __future__ import annotations

from typing import Union

import numpy as np
from anndata import AnnData

from .._typing import (
    AnnDataAxisWithBoth,
    anndata_checker,
)
from .._validation import _as_anndata_axis
from ._transfer import merge as merge


@anndata_checker
def sort_anndata(
    adata: AnnData,
    on: AnnDataAxisWithBoth = "both",
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Sort observations and/or variables by their AnnData index names.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    on: {"obs", "var", "both"} (default: "both")
        Axis to sort. If `"obs"`, sort observations by `adata.obs_names`. If
        `"var"`, sort variables by `adata.var_names`. If `"both"`, sort both
        axes.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with axes reordered.
        Otherwise, reorders `adata` in place and returns None.

        Reordering is applied to:

        - observations: if `on="obs"` or `on="both"`;
        - variables: if `on="var"` or `on="both"`.

    """

    adata = adata.copy() if copy else adata

    on = _as_anndata_axis(on, allow_both=True)

    if on in {"obs", "both"}:
        obs_order = np.argsort(adata.obs_names.to_numpy(), kind="stable")
        adata._inplace_subset_obs(obs_order)

    if on in {"var", "both"}:
        var_order = np.argsort(adata.var_names.to_numpy(), kind="stable")
        adata._inplace_subset_var(var_order)

    return adata if copy else None
