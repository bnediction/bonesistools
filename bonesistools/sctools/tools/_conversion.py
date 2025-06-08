#!/usr/bin/env python

from typing import Optional, Union, Sequence
from .._typing import (
    anndata_checker,
    Keys
)

from anndata import AnnData
from pandas import DataFrame

import numpy as np
from scipy.sparse import issparse

@anndata_checker
def anndata_to_dataframe(
    adata: AnnData,
    obs: Optional[Keys] = None,
    layer: Optional[str] = None,
    is_log: bool = False
) -> DataFrame:
    """
    Convert Anndata instance into Dataframe instance.
    
    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    obs: Keys (optional, default: None)
        Any key or key set in anndata.obs.
        If specified, add adata.obs.loc[obs] to dataframe.
    layer: str (optional, default: None)
        Any key in adata.layers.
        If provided, use adata.layers[layer] for expression values instead of adata.X.
    is_log: bool (default: False)
        Boolean value specifying if the counts are logarithmized.
        If value parameter is 'True', perform an exponential transformation.
        If counts are still logarithmized but user want to keep logarithmized counts,
        please specify 'False' to the value parameter.
    
    Returns
    -------
    Return DataFrame providing information about counts and optionnaly other additionnal chosen information.
    """

    if layer:
        if issparse(adata.layers[layer]):
            counts_df = DataFrame(adata.layers[layer].toarray(), index=adata.obs.index, columns = adata.var.index)
        else:
            counts_df = DataFrame(adata.layers[layer], index=adata.obs.index, columns = adata.var.index)
    elif issparse(adata.X):
        counts_df = DataFrame(adata.X.toarray(), index=adata.obs.index, columns = adata.var.index)
    else:
        counts_df = DataFrame(adata.X, index=adata.obs.index, columns = adata.var.index)
    
    if is_log:
        if "log1p" in adata.uns_keys() and adata.uns["log1p"].get('base') is not None:
            counts_df = np.expm1(counts_df * np.log(adata.uns['log1p']['base']))
        else:
            counts_df = np.expm1(counts_df)
    
    if obs is not None:
        counts_df.loc[:,obs] = adata.obs.loc[:,[obs]] if isinstance (obs, str) else adata.obs.loc[:,obs]
    
    return counts_df
