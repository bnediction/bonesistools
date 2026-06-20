#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Callable, Union, cast

from anndata import AnnData

from ..._validation import _as_callable, _as_string
from .._typing import anndata_checker


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
