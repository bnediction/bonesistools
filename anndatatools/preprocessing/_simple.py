#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    List,
    Tuple,
    Literal
)
from anndata import AnnData
from pandas import DataFrame
from numpy import nan
from .._typing import (
    adata_checker,
    type_checker,
    UnionType,
    DataFrameList,
    Axis,
    Keys,
    Suffixes
)

import pandas as pd

from ..tools import anndata_to_dataframe

@type_checker(dfs=UnionType(DataFrame,List))
def __generate_unique_index_name(
    dfs: Union[DataFrame,DataFrameList]
) -> str:
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

@adata_checker
def set_index(
    adata: AnnData,
    keys: Keys,
    axis: Axis = 0,
    copy: bool = False
) -> Union[AnnData, None]:
    """
    Create a MultiIndex for `adata.obs` or `adata.var` using current index and existing columns.

    Parameters
    ----------
    adata
        AnnData object
    keys
        either a single column key or a list containing an arbitrary combination of column keys
    axis
        whether to update index from adata.var (0 or `var`) or adata.obs (1 or `obs`)
    copy
        return a copy instead of updating `adata` object
        
    Returns
    -------
    Depending on `copy`, update or return AnnData object.
    """

    adata = adata.copy() if copy else adata
    
    if isinstance(keys, str):
        keys = [keys]
    elif not isinstance(keys, List):
        raise TypeError(f"unsupported argument type for keys: '{type(keys)}' instead of '{str}' or '{List}'")
    
    if axis in [0, "obs"]:
        df = adata.obs.copy()
    elif axis in [1, "var"]:
        df = adata.var.copy()
    else:
        raise TypeError(f"unsupported argument type for axis: `{axis}`")
    
    index_name = __generate_unique_index_name(df)
    df[index_name] = df.index
    df = df.set_index([index_name, *keys])

    if axis in [0, "obs"]:
        adata.obs = df
    elif axis in [1, "var"]:
        adata.var = df
    
    return adata if copy else None

@adata_checker(n=2)
def merge(
    left_ad: AnnData,
    right_ad: AnnData,
    axis: Axis = 0,
    suffixes: Suffixes = ("_x", "_y"),
    copy: bool = False
) -> Union[AnnData, None]:
    """
    Merge DataFrame from `adata.obs` or `adata.var` with an index-based join.

    Parameters
    ----------
    left_ad
        AnnData object receiving information
    right_ad
        AnnData object sending information
    axis
        whether to merge adata.var (0 or `var`) or adata.obs (1 or `obs`)
    suffixes
        length-2 sequence where each element is a string indicating the suffix
        to add to overlapping column names in `left_ad` and `right_ad` respectively
    copy
        return a copy instead of updating `left_ad` object
    
    Returns
    -------
    Depending on `copy`, update or return `left_ad` AnnData object.
    """
    
    left_ad = left_ad.copy() if copy else left_ad

    if axis in [0, "obs"]:
        left_df = left_ad.obs.copy()
        right_df = right_ad.obs.copy()
    elif axis in [1, "var"]:
        left_df = left_ad.var.copy()
        right_df = right_ad.var.copy()
    else:
        raise TypeError(f"unsupported argument type for axis: `{axis}`")

    df = left_df.merge(
        right=right_df,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=suffixes
    )

    if axis in [0, "obs"]:
        left_ad.obs = df
    elif axis in [1, "var"]:
        left_ad.var = df

    return left_ad if copy else None

@adata_checker(n=2)
def transfer_layer(
    left_ad: AnnData,
    right_ad: AnnData,
    layers: Keys,
    copy: bool = False
) -> Union[AnnData, None]:
    """
    Transfer layers from `right_ad.layers` to `left_ad.layers` by preserving
    the order of observations and variables.

    Parameters
    ----------
    left_ad
        AnnData object receiving layer-based information
    right_ad
        AnnData object sending layer-based information
    layers
        sequence where each element is a string indicating the layer to add in `left_ad`
    copy
        return a copy instead of updating `left_ad` object
    
    Returns
    -------
    Depending on `copy`, update or return `left_ad` AnnData object.
    """
    
    left_ad = left_ad.copy() if copy else left_ad

    if isinstance(layers, str):
        layers = [layers]
    elif not isinstance(layers, List):
        raise TypeError(f"unsupported argument type for layers: '{type(layers)}' instead of '{str}' or '{List}'")

    for layer in layers:
        df = anndata_to_dataframe(
            adata=right_ad,
            obs=None,
            layer=layer
        )
        left_cols = set(left_ad.var.index)
        left_idx = set(left_ad.obs.index)
        right_cols = set(right_ad.var.index)
        right_idx = set(right_ad.obs.index)

        cols_to_add = list(left_cols.difference(right_cols))
        idx_to_add = list(left_idx.difference(right_idx))
        if cols_to_add:
            df.loc[:,cols_to_add] = nan
        for idx in idx_to_add:
            df.loc[idx,:] = pd.Series(dtype="float64")

        cols_to_remove = list(right_cols.difference(left_cols))
        index_to_remove = list(right_idx.difference(left_idx))
        df.drop(
            columns=cols_to_remove,
            index=index_to_remove,
            inplace=True
        )

        left_ad.layers[layer] = df[left_ad.var.index.tolist()].reindex(left_ad.obs.index.tolist())

    return left_ad if copy else None

@adata_checker
def transfer_obs_sti(
    adata: AnnData,
    adatas: List[AnnData],
    obs: List[str],
    conditions: List[str],
    condition_colname: str = "condition",
    copy: bool = False
) -> Union[AnnData, None]:
    """
    Transfer observations from specific to integrated dataset, i.e. transfer columns from multiple AnnData `adatas` towards a unique AnnData `adata`.
    This function handles issues whenever there are identical indices in `adatas.obs`.

    Parameters
    ----------
    adata
        AnnData object receiving information (specific datasets)
    adatas
        AnnData objects sending information (integrated dataset)
    obs
        column names in specific `adata.obs` to transfer
    conditions
        conditions related to AnnData objects (ordered w.r.t `adatas`)
    condition_colname
        column name in integrated `adata.obs` related to conditions
    copy
        return a copy instead of updating 'adata' object
        
    Returns
    -------
    Depending on `copy`, update or return AnnData object with new observations.
    """
    
    adata = adata.copy() if copy else adata
    
    all_samples_df = []
    for _condition, _adata in zip(conditions, adatas):
        _df = _adata.obs.loc[:,obs].copy()
        _df[condition_colname] = _condition
        all_samples_df.append(_df)
        del _df

    all_samples_df = pd.concat(all_samples_df)
    one_sample_df = adata.obs.copy()

    index_name = __generate_unique_index_name([all_samples_df, one_sample_df])

    all_samples_df[index_name] = all_samples_df.index
    all_samples_df = all_samples_df.set_index([index_name, condition_colname])

    one_sample_df[index_name] = one_sample_df.index
    one_sample_df = one_sample_df.set_index([index_name, condition_colname])

    merge_df = one_sample_df.merge(
        right=all_samples_df,
        how="left",
        left_index=True,
        right_index=True
    )
    merge_df.reset_index(
        level=(condition_colname,),
        inplace=True
    )
    merge_df.index.name = None

    adata.obs = merge_df

    return adata if copy else None

@adata_checker
def transfer_obs_its(
    adata: AnnData,
    adatas: List[AnnData],
    obs: List[str],
    conditions: List[str],
    condition_colname: str = "condition",
    copy: bool = False
) -> Union[AnnData, None]:
    """
    Transfer observations from integrated to specific datasets, i.e. transfer columns from a unique AnnData `adata` towards multiple AnnData `adatas`.

    Parameters
    ----------
    adata
        AnnData object receiving information (specific datasets)
    adatas
        AnnData objects sending information (integrated dataset)
    obs
        column names in integrated `adata.obs` to transfer
    conditions
        conditions related to AnnData objects (ordered w.r.t `adatas`)
    condition_colname
        column name in integrated `adata.obs` related to conditions
    copy
        return a copy instead of updating 'adata' object
        
    Returns
    -------
    Depending on `copy`, update AnnData objects or return a list of AnnData objects with new observations.
    """

    if copy:
        for _adata in adatas:
            _adata = _adata.copy()
        adatas_cp = []

    for _condition, _adata in zip(conditions, adatas):
        _cond = adata.obs[condition_colname] == _condition
        df = adata.obs.loc[_cond][obs]
        _adata.obs = _adata.obs.merge(
            right=df,
            how="left",
            left_index=True,
            right_index=True
        )
        if copy:
            adatas_cp.append(_adata)
    
    return adatas_cp if copy else None
