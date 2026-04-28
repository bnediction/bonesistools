#!/usr/bin/env python

from typing import Union, Optional
from .._typing import (
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
    Axis
)
from pandas import DataFrame
from anndata import AnnData
from scipy.sparse import csr_matrix

import anndata as ad

from ...databases.ncbi import (
    InputIdentifierType,
    OutputIdentifierType,
    GeneSynonyms,
    support_legacy_gene_synonyms_args
)

@support_legacy_gene_synonyms_args
@anndata_or_mudata_checker
def convert_gene_identifiers(
    scdata: ScData, # type: ignore
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    copy: bool = False
) -> Union[ScData, None]: # type: ignore
    f"""
    Replace gene identifiers into the desired identifier format.

    Parameters
    ----------
    scdata: ad.AnnData | md.MuData
        Unimodal or multimodal annotated data matrix.
        The stored gene identifiers are converted into the desired identifier format.
    axis: Axis (default: 'var')
        whether to rename labels from scdata.var (0 or 'obs') or scdata.obs (1 or 'var').
    input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
        Gene identifier input format.
    output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
        Gene identifier output format.
    copy: bool (default: False)
        Return a copy instead of updating ScData object.
    
    Returns
    -------
    Depending on 'copy', update 'scdata' or return ScData object.
    """

    scdata = scdata.copy() if copy else scdata

    if axis in [0, "obs"]:
        GeneSynonyms()(
            scdata.obs,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False
        )
    elif axis in [1, "var"]:
        GeneSynonyms()(
            scdata.var,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False
        )
    else:
        raise TypeError(f"unsupported argument type for 'axis': {axis}")
    
    return scdata if copy else None

@support_legacy_gene_synonyms_args
@anndata_checker
def standardize_gene_identifiers(
    scdata: ScData, # type: ignore
    axis: Axis = "var",
    copy: bool = False
) -> Union[ScData, None]: # type: ignore
    """
    Standardize gene names by converting them into their corresponding official names.

    Parameters
    ----------
    scdata: ad.AnnData | md.MuData
        Unimodal or multimodal annotated data matrix.
        The stored gene names are standardized.
    axis: Axis (default: 'var')
        whether to rename labels from scdata.var (0 or 'obs') or scdata.obs (1 or 'var').
    copy: bool (default: False)
        Return a copy instead of updating ScData object.
        
    Returns
    -------
    Depending on 'copy', update 'scdata' or return ScData object.
    """

    return convert_gene_identifiers(
        scdata=scdata,
        axis=axis,
        input_identifier_type="name",
        output_identifier_type="official_name",
        copy=copy
    )

@anndata_checker
def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None
) -> Union[AnnData, None]:
    """
    Merge the duplicated index names in adata.var by summing the counts between each duplicated index elements.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix where duplicated variable names are merged.
    var_names_column: str (optional, default: None)
        If specify, give a priority order for which row in adata.var is saved when rows are merged.
        For instance, if multiple rows have same index but different cell values,
        the row with cell value at 'var_names_column' column is considered as the main information over others.
    
    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.
    """

    if var_names_column is None:
        var_names = "copy_var_names"
        adata.var[var_names] = list(adata.var.index)
    else:
        var_names = var_names_column

    obs = DataFrame(adata.obs.index).set_index(0)
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
        adata_spec = AnnData(X=X, var=var, obs=obs)
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
