#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Callable, Union, cast

from anndata import AnnData
from pandas import DataFrame

from .._typing import (
    Axis,
    Suffixes,
    anndata_checker,
)


@anndata_checker
def filter_obs(
    adata: AnnData,
    obs: str,
    function: Callable[[object], object],
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
        Filtered AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    KeyError
        If `obs` is not found in `adata.obs`.
    TypeError
        If `obs` is not a string or `function` is not callable.
    """

    adata = adata.copy() if copy else adata

    if isinstance(obs, str):
        if obs not in adata.obs:
            raise KeyError(f"key '{obs}' not found in adata.obs")
    else:
        raise TypeError(
            f"unsupported argument type for 'obs': "
            f"expected {str} but received {type(obs)}"
        )

    if not callable(function):
        raise TypeError(
            f"unsupported argument type for 'function': "
            f"expected callable object but received {type(function)}"
        )

    obs_subset = function(adata.obs[obs].values)
    adata._inplace_subset_obs(cast(Any, obs_subset))

    return adata if copy else None


@anndata_checker
def filter_var(
    adata: AnnData,
    var: str,
    function: Callable[[object], object],
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
        Filtered AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    KeyError
        If `var` is not found in `adata.var`.
    TypeError
        If `var` is not a string or `function` is not callable.
    """

    adata = adata.copy() if copy else adata

    if isinstance(var, str):
        if var not in adata.var:
            raise KeyError(f"key '{var}' not found in adata.var")
    else:
        raise TypeError(
            f"unsupported argument type for 'var': "
            f"expected {str} but received {type(var)}"
        )

    if not callable(function):
        raise TypeError(
            f"unsupported argument type for 'function': "
            f"expected callable object but received {type(function)}"
        )

    var_subset = function(adata.var[var].values)
    adata._inplace_subset_var(cast(Any, var_subset))

    return adata if copy else None


@anndata_checker(n=2)
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = 0,
    suffixes: Suffixes = ("_x", "_y"),
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Merge annotation tables from two AnnData objects with an index-based join.

    The left AnnData object receives columns from the right AnnData object.
    Observation annotations are merged when `axis=0` or `"obs"`; variable
    annotations are merged when `axis=1` or `"var"`.

    Parameters
    ----------
    left_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object receiving new information.
    right_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object sending information.
    axis: {0, 1, "obs", "var"} (default: 0)
        If 0 or `"obs"`, merge `.obs`. If 1 or `"var"`, merge `.var`.
    suffixes: Tuple[str, str] (default: ('_x','_y'))
        Length-2 sequence where each element is a string indicating the suffix
        to add to overlapping column names in `left_ad` and `right_ad`,
        respectively.
    copy: bool (default: False)
        Return a copy instead of modifying `left_ad`.

    Returns
    -------
    AnnData or None
        Merged AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    ValueError
        If `axis` is not 0, 1, `"obs"` or `"var"`.
    """

    left_ad = left_ad.copy() if copy else left_ad

    if axis in [0, "obs"]:
        left_obs = cast(DataFrame, left_ad.obs)
        right_obs = cast(DataFrame, right_ad.obs)
        left_df = left_obs.copy()
        right_df = right_obs.copy()
    elif axis in [1, "var"]:
        left_var = cast(DataFrame, left_ad.var)
        right_var = cast(DataFrame, right_ad.var)
        left_df = left_var.copy()
        right_df = right_var.copy()
    else:
        raise ValueError(
            f"invalid argument value for 'axis': "
            f"expected 0, 1, 'obs' or 'var' but received {axis!r}"
        )

    df = left_df.merge(
        right=right_df, how="left", left_index=True, right_index=True, suffixes=suffixes
    )

    if axis in [0, "obs"]:
        left_ad.obs = df
    elif axis in [1, "var"]:
        left_ad.var = df

    return left_ad if copy else None
