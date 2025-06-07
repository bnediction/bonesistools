#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Union, Optional
from .._typing import (
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
    Axis
)

import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

from ...databases.ncbi import GeneSynonyms

@anndata_or_mudata_checker
def gene_synonyms_conversion(
    scdata: ScData, # type: ignore
    axis: Axis="var",
    input_type: str="genename",
    output_type: str="referencename",
    keep_if_missing: bool=True,
    copy: bool=False
) -> Union[ScData, None]: # type: ignore
    """
    Replace gene names with theirs gene aliases into 'scdata' object.

    Parameters
    ----------
    scdata
        AnnData or MuData object where names are expected being standardized
    axis
        whether to rename labels from scdata.var (0 or 'obs') or scdata.obs (1 or 'var')
    input_type
        genename|geneid|ensemblid|<database>
    output_type
        referencename|geneid|ensemblid|<database>
    copy
        return a copy instead of updating 'scdata' object
        
    Returns
    -------
    Depending on 'copy', update or return AnnData or MuData object with gene aliases.
    """

    scdata = scdata.copy() if copy else scdata

    if axis in [0, "obs"]:
        GeneSynonyms()(
            scdata.obs,
            axis="index",
            input_type=input_type,
            output_type=output_type,
            keep_if_missing=keep_if_missing,
            copy=False
        )
    elif axis in [1, "var"]:
        GeneSynonyms()(
            scdata.var,
            axis="index",
            input_type=input_type,
            output_type=output_type,
            keep_if_missing=keep_if_missing,
            copy=False
        )
    else:
        raise TypeError(f"unsupported argument type for 'axis': {axis}")
    
    return scdata if copy else None

@anndata_checker
def set_ncbi_reference_name(
    scdata: ScData, # type: ignore
    axis: Axis="var",
    input_type: str="genename",
    keep_if_missing: bool=True,
    copy: bool = False
) -> Union[ScData, None]: # type: ignore
    """
    Replace gene names with theirs reference gene names into 'scdata' object.

    Parameters
    ----------
    scdata
        AnnData or MuData object where names are expected being standardized
    annotations
        whether to rename labels from scdata.var (0 or 'obs') or scdata.obs (1 or 'var')
    input_type
        genename|geneid|ensemblid|<database>
    copy
        return a copy instead of updating 'scdata' object
        
    Returns
    -------
    Depending on 'copy', update or return AnnData or MuData object with reference gene names.
    """

    scdata = scdata.copy() if copy else scdata

    gene_synonyms_conversion(
        scdata=scdata,
        axis=axis,
        input_type=input_type,
        output_type="referencename",
        keep_if_missing=keep_if_missing,
        copy=False
    )

    return scdata if copy else None

@anndata_checker
def var_names_merge_duplicates(
    adata: ad.AnnData,
    var_names_column: Optional[str]=None
) -> Union[ad.AnnData, None]:
    """
    Merge the duplicated index names in adata.var by summing the counts between each duplicated index elements.

    Parameters
    ----------
    adata
        AnnData object where names are expected being standardized
    var_names_column
        if specify, give a priority order for which row in adata.var is saved when rows are merged.
        For instance, if multiple rows have same index but different cell values,
        the row with cell value at 'var_names_column' column is considered as the main information over others.
    
    Returns
    -------
    return AnnData object with duplicated index names in adata.var merged together.
    """

    if var_names_column is None:
        var_names = "copy_var_names"
        adata.var[var_names] = list(adata.var.index)
    else:
        var_names = var_names_column

    obs = pd.DataFrame(adata.obs.index).set_index(0)
    adatas = list()
    duplicated_var_names = {adata.var_names[idx] for idx, value in enumerate(adata.var_names.duplicated()) if value}

    for var_name in duplicated_var_names:
        adata_spec = adata[:,adata.var.index == var_name]
        adata = adata[:,adata.var.index != var_name]

        X = csr_matrix(adata_spec.X.sum(axis=1))
        if var_name in list(adata_spec.var[var_names]):
            filter = adata_spec.var[var_names] == var_name
            var = adata_spec.var[filter].iloc[:1]
        else:
            var = adata_spec.var.iloc[:1]
        adata_spec = ad.AnnData(X=X, var=var, obs=obs)
        adatas.append(adata_spec)

    adatas.append(adata)

    adata = ad.concat(
        adatas=adatas,
        join="outer",
        axis=1,
        merge="same",
        uns_merge="first",
        label=None
    )

    if var_names_column is None:
        adata.var.drop(labels=var_names, axis="columns", inplace=True)

    return adata
