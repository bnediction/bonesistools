#!/usr/bin/env python

from __future__ import annotations

from typing import Any, List, Union, cast, overload

import pandas as pd
from anndata import AnnData
from numpy import nan
from pandas import DataFrame

from ..._compat import Literal
from ..._warnings import _warn_deprecated
from .._typing import (
    AnnDataList,
    Axis,
    DataFrameList,
    Keys,
    Suffixes,
    UnionType,
    anndata_checker,
    type_checker,
)
from .._validation import _as_anndata_axis
from ..tools import to_dataframe


@overload
def transfer_layer(
    left_ad: AnnData,
    right_ad: AnnData,
    layers: Keys,
    *,
    copy: Literal[False],
) -> None: ...


@overload
def transfer_layer(
    left_ad: AnnData,
    right_ad: AnnData,
    layers: Keys,
    copy: Literal[True] = True,
) -> AnnData: ...


@overload
def transfer_layer(
    left_ad: AnnData,
    right_ad: AnnData,
    layers: Keys,
    copy: bool = False,
) -> Union[AnnData, None]: ...


@anndata_checker(n=2)
def transfer_layer(
    left_ad: AnnData,
    right_ad: AnnData,
    layers: Keys,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Transfer layers from 'right_ad.layers' to 'left_ad.layers' by preserving
    the order of observations and variables.

    Missing observations or variables in `right_ad` are filled with NaN in the
    transferred layer. Extra observations or variables from `right_ad` are
    discarded so that the resulting layer follows `left_ad.obs_names` and
    `left_ad.var_names`.

    Parameters
    ----------
    left_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object receiving layer-based information.
    right_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object sending layer-based information.
    layers: Keys
        Layer name or sequence of layer names to transfer into `left_ad.layers`.
    copy: bool (default: False)
        Return a copy instead of modifying `left_ad`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `left_ad` with layers transferred.
        Otherwise, updates `left_ad` in place and returns None.

        Transferred layer data are stored in:

        - `left_ad.layers[layer]`: transferred matrix for each requested layer.

    """

    left_ad = left_ad.copy() if copy else left_ad

    if isinstance(layers, str):
        layers = [layers]
    elif not isinstance(layers, List):
        raise TypeError(
            f"unsupported argument type for 'layers': "
            f"expected {str} or {list} but received {type(layers)}"
        )

    for layer in layers:
        df = to_dataframe(adata=right_ad, obs=None, layer=layer)
        left_cols = set(left_ad.var.index)
        left_idx = set(left_ad.obs.index)
        right_cols = set(right_ad.var.index)
        right_idx = set(right_ad.obs.index)

        cols_to_add = list(left_cols.difference(right_cols))
        idx_to_add = list(left_idx.difference(right_idx))
        if cols_to_add:
            df.loc[:, cols_to_add] = nan
        df = pd.concat(
            [
                df,
                pd.DataFrame(
                    data=nan,
                    columns=df.columns,
                    index=cast(Any, idx_to_add),
                ),
            ]
        )

        cols_to_remove = list(right_cols.difference(left_cols))
        index_to_remove = list(right_idx.difference(left_idx))
        df.drop(columns=cols_to_remove, index=index_to_remove, inplace=True)

        left_ad.layers[layer] = cast(
            Any,
            df[left_ad.var.index.tolist()].reindex(left_ad.obs.index.tolist()),
        )

    return left_ad if copy else None


@overload
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = "obs",
    suffixes: Suffixes = ("_x", "_y"),
    *,
    copy: Literal[False],
) -> None: ...


@overload
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = "obs",
    suffixes: Suffixes = ("_x", "_y"),
    copy: Literal[True] = True,
) -> AnnData: ...


@overload
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = "obs",
    suffixes: Suffixes = ("_x", "_y"),
    copy: bool = False,
) -> Union[AnnData, None]: ...


@anndata_checker(n=2)
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = "obs",
    suffixes: Suffixes = ("_x", "_y"),
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Merge annotation tables from two AnnData objects with an index-based join.

    The left AnnData object receives columns from the right AnnData object.
    Observation annotations are merged when `axis="obs"`; variable annotations
    are merged when `axis="var"`.

    Parameters
    ----------
    left_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object receiving new information.
    right_ad: AnnData
        Unimodal annotated data matrix.
        It corresponds to the object sending information.
    axis: {"obs", "var"} (default: "obs")
        If `"obs"`, merge `.obs`. If `"var"`, merge `.var`. Deprecated values
        0 and 1 are still accepted as aliases for `"obs"` and `"var"`.
    suffixes: Tuple[str, str] (default: ('_x','_y'))
        Length-2 sequence where each element is a string indicating the suffix
        to add to overlapping column names in `left_ad` and `right_ad`,
        respectively.
    copy: bool (default: False)
        Return a copy instead of modifying `left_ad`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `left_ad` with annotations merged.
        Otherwise, updates `left_ad` in place and returns None.

        Merged annotations are stored in:

        - `left_ad.obs`: merged observation annotations if `axis="obs"`;
        - `left_ad.var`: merged variable annotations if `axis="var"`.

    """

    left_ad = left_ad.copy() if copy else left_ad

    axis = _as_anndata_axis(axis, allow_integer=True)
    if axis == "obs":
        left_obs = cast(DataFrame, left_ad.obs)
        right_obs = cast(DataFrame, right_ad.obs)
        left_df = left_obs.copy()
        right_df = right_obs.copy()
    elif axis == "var":
        left_var = cast(DataFrame, left_ad.var)
        right_var = cast(DataFrame, right_ad.var)
        left_df = left_var.copy()
        right_df = right_var.copy()

    df = left_df.merge(
        right=right_df, how="left", left_index=True, right_index=True, suffixes=suffixes
    )

    if axis == "obs":
        left_ad.obs = df
    elif axis == "var":
        left_ad.var = df

    return left_ad if copy else None


@overload
def transfer_obs_to_integrated(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    *,
    copy: Literal[False],
) -> None: ...


@overload
def transfer_obs_to_integrated(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: Literal[True] = True,
) -> AnnData: ...


@overload
def transfer_obs_to_integrated(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnData, None]: ...


@anndata_checker
def transfer_obs_to_integrated(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Transfer observation columns from specific AnnData objects to an
    integrated AnnData object.

    This function handles cases where identical observation indices are shared
    by several AnnData objects.

    The merge is performed on observation index and condition, allowing
    condition-specific AnnData objects to contain overlapping observation names.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
        It corresponds to the integrated dataset receiving information.
    adatas: AnnDataList
        List of unimodal annotated data matrices.
        It corresponds to the specific datasets sending information.
    obs: Keys
        Observation column name or names from specific `adatas.obs` to
        transfer.
    conditions: Keys
        Conditions associated with AnnData objects, ordered like `adatas`.
    condition_colname: str (default: 'condition')
        Column name in integrated `adata.obs` storing conditions.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with observation
        annotations transferred. Otherwise, updates `adata` in place and
        returns None.

        Transferred annotations are stored in:

        - `adata.obs`: integrated observation annotations.
    """

    adata = adata.copy() if copy else adata
    obs = [obs] if isinstance(obs, str) else obs

    all_samples_dfs: List[DataFrame] = []
    for _condition, _adata in zip(conditions, adatas):
        _adata_obs = cast(DataFrame, _adata.obs)
        _df = _adata_obs.loc[:, obs].copy()
        _df[condition_colname] = _condition
        all_samples_dfs.append(_df)
        del _df

    all_samples_df = pd.concat(all_samples_dfs)
    adata_obs = cast(DataFrame, adata.obs)
    one_sample_df = adata_obs.copy()

    index_name = __generate_unique_index_name([all_samples_df, one_sample_df])

    all_samples_df[index_name] = all_samples_df.index
    all_samples_df = all_samples_df.set_index([index_name, condition_colname])

    one_sample_df[index_name] = one_sample_df.index
    one_sample_df = one_sample_df.set_index([index_name, condition_colname])

    merge_df = one_sample_df.merge(
        right=all_samples_df, how="left", left_index=True, right_index=True
    )
    merge_df.reset_index(level=(condition_colname,), inplace=True)
    merge_df.index.name = None

    adata.obs = merge_df

    return adata if copy else None


@type_checker(dfs=UnionType(DataFrame, List))
def __generate_unique_index_name(dfs: Union[DataFrame, DataFrameList]) -> str:

    dfs = [dfs] if isinstance(dfs, DataFrame) else dfs
    column_names = set()
    for df in dfs:
        column_names.update(set(df.columns))
    index_name = "index"
    _i = 0
    while index_name in column_names:
        index_name = f"index_{_i}"
        _i += 1
    return index_name


@overload
def transfer_obs_to_specific(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    *,
    copy: Literal[False],
) -> None: ...


@overload
def transfer_obs_to_specific(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: Literal[True] = True,
) -> AnnDataList: ...


@overload
def transfer_obs_to_specific(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnDataList, None]: ...


@anndata_checker
def transfer_obs_to_specific(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnDataList, None]:
    """
    Transfer observation columns from an integrated AnnData object to
    specific AnnData objects.

    Rows are matched by observation index after splitting the integrated
    AnnData object by `condition_colname`.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
        It corresponds to the integrated dataset sending information.
    adatas: AnnDataList
        List of unimodal annotated data matrices.
        It corresponds to the specific datasets receiving information.
    obs: Keys
        Observation column name or names from integrated `adata.obs` to
        transfer.
    conditions: Keys
        Conditions associated with AnnData objects, ordered like `adatas`.
    condition_colname: str (default: 'condition')
        Column name in integrated `adata.obs` storing conditions.
    copy: bool (default: False)
        Return a copy instead of modifying `adatas`.

    Returns
    -------
    list of AnnData or None
        If `copy=True`, returns copies of `adatas` with observation
        annotations transferred. Otherwise, updates each AnnData object in
        place and returns None.

        Transferred annotations are stored in:

        - `adata.obs`: specific observation annotations for each AnnData object.
    """

    adatas_cp: AnnDataList = []

    if copy:
        adatas = [_adata.copy() for _adata in adatas]

    for _condition, _adata in zip(conditions, adatas):
        adata_obs = cast(DataFrame, adata.obs)
        _adata_obs = cast(DataFrame, _adata.obs)
        _cond = adata_obs[condition_colname] == _condition
        df = adata_obs.loc[_cond][obs]
        _adata.obs = _adata_obs.merge(
            right=df, how="left", left_index=True, right_index=True
        )
        if copy:
            adatas_cp.append(_adata)

    return adatas_cp if copy else None


def transfer_obs_sti(*args: Any, **kwargs: Any) -> Union[AnnData, None]:
    """
    Deprecated alias for `transfer_obs_to_integrated`.
    """

    _warn_deprecated(
        "`bt.sct.pp.transfer_obs_sti`",
        replacement="`bt.sct.pp.transfer_obs_to_integrated`",
        stacklevel=2,
    )
    return transfer_obs_to_integrated(*args, **kwargs)


def transfer_obs_its(*args: Any, **kwargs: Any) -> Union[AnnDataList, None]:
    """
    Deprecated alias for `transfer_obs_to_specific`.
    """

    _warn_deprecated(
        "`bt.sct.pp.transfer_obs_its`",
        replacement="`bt.sct.pp.transfer_obs_to_specific`",
        stacklevel=2,
    )
    return transfer_obs_to_specific(*args, **kwargs)
