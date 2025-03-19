#!/usr/bin/env python

from typing import (
    Union,
    List
)

import anndata as ad
import pandas as pd

def __generate_unique_index_name(
    dfs: Union[List[pd.DataFrame],pd.DataFrame]
) -> str:
    dfs = [dfs] if isinstance(dfs, pd.DataFrame) else dfs
    column_names = set()
    for df in dfs:
        column_names.update(set(df.columns))
    index_name = "index"
    _i = 0
    while index_name in column_names:
        index_name = f"{index_name}_{_i}"
        _i += 1
    return index_name

def transfer_obs_sti(
    adata: ad.AnnData,
    adatas: List[ad.AnnData],
    obs: List[str],
    conditions: List[str],
    condition_colname: str = "condition",
    copy: bool = False
) -> Union[ad.AnnData, None]:
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

def transfer_obs_its(
    adata: ad.AnnData,
    adatas: List[ad.AnnData],
    obs: List[str],
    conditions: List[str],
    condition_colname: str = "condition",
    copy: bool = False
) -> Union[ad.AnnData, None]:
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
