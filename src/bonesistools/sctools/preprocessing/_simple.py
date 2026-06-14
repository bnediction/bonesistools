#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Callable, Union, cast

import numpy as np
from anndata import AnnData

from ..._validation import _as_callable, _as_string
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


@anndata_checker
def filter_obs(
    adata: AnnData,
    obs: str,
    function: Callable[[Any], Any],
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Filter observations based on a column in `adata.obs`.

    Observations are kept when `function(adata.obs[obs].values)` evaluates to
    True for the corresponding row.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    obs: str
        Column name in `adata.obs` used for filtering.
    function: Callable
        Function applied to `adata.obs[obs].values`. It must return a Boolean
        mask compatible with observations.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with observations filtered.
        Otherwise, subsets `adata` in place and returns None.

        Filtering is applied to:

        - observations: rows matching the Boolean mask.

    """

    adata = adata.copy() if copy else adata

    obs = _as_string(obs, "obs")
    function = _as_callable(function, "function")

    obs_subset = function(adata.obs[obs].values)
    adata._inplace_subset_obs(cast(Any, obs_subset))

    return adata if copy else None


@anndata_checker
def filter_var(
    adata: AnnData,
    var: str,
    function: Callable[[Any], Any],
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Filter variables based on a column in `adata.var`.

    Variables are kept when `function(adata.var[var].values)` evaluates to
    True for the corresponding column.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    var: str
        Column name in `adata.var` used for filtering.
    function: Callable
        Function applied to `adata.var[var].values`. It must return a Boolean
        mask compatible with variables.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with variables filtered.
        Otherwise, subsets `adata` in place and returns None.

        Filtering is applied to:

        - variables: columns matching the Boolean mask.

    """

    adata = adata.copy() if copy else adata

    var = _as_string(var, "var")
    function = _as_callable(function, "function")

    var_subset = function(adata.var[var].values)
    adata._inplace_subset_var(cast(Any, var_subset))

    return adata if copy else None
