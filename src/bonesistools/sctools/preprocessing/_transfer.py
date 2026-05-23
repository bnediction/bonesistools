#!/usr/bin/env python

from typing import List, Union

import pandas as pd
from anndata import AnnData
from numpy import nan
from pandas import DataFrame

from .._typing import (
    AnnDataList,
    DataFrameList,
    Keys,
    UnionType,
    anndata_checker,
    type_checker,
)
from ..tools import anndata_to_dataframe


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
        Updated AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    TypeError
        If `layers` is neither a string nor a list.
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
        df = anndata_to_dataframe(adata=right_ad, obs=None, layer=layer)
        left_cols = set(left_ad.var.index)
        left_idx = set(left_ad.obs.index)
        right_cols = set(right_ad.var.index)
        right_idx = set(right_ad.obs.index)

        cols_to_add = list(left_cols.difference(right_cols))
        idx_to_add = list(left_idx.difference(right_idx))
        if cols_to_add:
            df.loc[:, cols_to_add] = nan
        df = pd.concat(
            [df, pd.DataFrame(data=nan, columns=df.columns, index=idx_to_add)]
        )

        cols_to_remove = list(right_cols.difference(left_cols))
        index_to_remove = list(right_idx.difference(left_idx))
        df.drop(columns=cols_to_remove, index=index_to_remove, inplace=True)

        left_ad.layers[layer] = df[left_ad.var.index.tolist()].reindex(
            left_ad.obs.index.tolist()
        )

    return left_ad if copy else None


@anndata_checker
def transfer_obs_sti(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Transfer observations from specific to integrated dataset, i.e. transfer columns
    from multiple AnnData 'adatas' towards a unique AnnData 'adata'.
    This function handles issues whenever there are identical indices between AnnData objects.

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
        Updated integrated AnnData object if `copy=True`; otherwise None.
    """

    adata = adata.copy() if copy else adata
    obs = [obs] if isinstance(obs, str) else obs

    all_samples_df = []
    for _condition, _adata in zip(conditions, adatas):
        _df = _adata.obs.loc[:, obs].copy()
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
        right=all_samples_df, how="left", left_index=True, right_index=True
    )
    merge_df.reset_index(level=(condition_colname,), inplace=True)
    merge_df.index.name = None

    adata.obs = merge_df

    return adata if copy else None


@anndata_checker
def transfer_obs_its(
    adata: AnnData,
    adatas: AnnDataList,
    obs: Keys,
    conditions: Keys,
    condition_colname: str = "condition",
    copy: bool = False,
) -> Union[AnnDataList, None]:
    """
    Transfer observations from integrated to specific datasets,
    i.e. transfer columns from a unique AnnData 'adata' towards multiple AnnData 'adatas'.

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
        Updated AnnData objects if `copy=True`; otherwise None.
    """

    if copy:
        adatas = [_adata.copy() for _adata in adatas]
        adatas_cp = []

    for _condition, _adata in zip(conditions, adatas):
        _cond = adata.obs[condition_colname] == _condition
        df = adata.obs.loc[_cond][obs]
        _adata.obs = _adata.obs.merge(
            right=df, how="left", left_index=True, right_index=True
        )
        if copy:
            adatas_cp.append(_adata)

    return adatas_cp if copy else None
