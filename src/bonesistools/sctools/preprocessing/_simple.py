#!/usr/bin/env python

from typing import Callable, Union

from anndata import AnnData

from .._typing import (
    Axis,
    Suffixes,
    anndata_checker,
)


@anndata_checker
def filter_obs(
    adata: AnnData, obs: str, function: Callable, copy: bool = False
) -> Union[AnnData, None]:
    """
    Filter observations based on a column in 'adata.obs'.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    obs: str
        Column name in 'adata.obs' used for filtering.
    function: Callable
        Function to apply to the observation used for filtering.
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.

    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.

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
    adata._inplace_subset_obs(obs_subset)

    return adata if copy else None


@anndata_checker
def filter_var(
    adata: AnnData, var: str, function: Callable, copy: bool = False
) -> Union[AnnData, None]:
    """
    Filter variables based on a column in 'adata.var'.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    var: str
        Column name in 'adata.var' used for filtering.
    function: Callable
        Function to apply to the variable used for filtering.
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.

    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.

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
    adata._inplace_subset_var(var_subset)

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
    Merge dataframes from 'adata.obs' or 'adata.var' with an index-based join.

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
        to add to overlapping column names in 'left_ad' and 'right_ad' respectively.
    copy: bool (default: False)
        Return a copy instead of updating 'left_ad' object.

    Returns
    -------
    Depending on 'copy', update 'left_ad' or return AnnData object.

    Raises
    ------
    ValueError
        If `axis` is not 0, 1, `"obs"` or `"var"`.
    """

    left_ad = left_ad.copy() if copy else left_ad

    if axis in [0, "obs"]:
        left_df = left_ad.obs.copy()
        right_df = right_ad.obs.copy()
    elif axis in [1, "var"]:
        left_df = left_ad.var.copy()
        right_df = right_ad.var.copy()
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
