#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Union, Optional

import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

from databases.genesyn import GeneSynonyms

def gene_synonyms_conversion(
    adata: ad.AnnData,
    annotations: str="var",
    in_alias_type: str="genename",
    out_alias_type: str="referencename",
    keep_if_missing: bool=True,
    copy: bool = False
) -> Union[ad.AnnData, None]:
    """
    Replace gene names with theirs gene aliases into 'adata' object.

    Parameters
    ----------
    adata
        AnnData object where names are expected being standardized
    annotations
        whether to rename labels from adata.var (0 or 'var') or adata.obs (1 or 'obs')
    in_alias_type
        genename|geneid|ensemblid|<database>
    out_alias_type
        referencename|geneid|ensemblid|<database>
    copy
        return a copy instead of updating 'adata' object
        
    Returns
    -------
    Depending on 'copy', update or return AnnData object with gene aliases.
    """

    adata = adata.copy() if copy else adata

    if annotations in [0, "var"]:
        GeneSynonyms()(
            adata.var,
            axis="index",
            in_alias_type=in_alias_type,
            out_alias_type=out_alias_type,
            keep_if_missing=keep_if_missing,
            copy=False
        )
    elif annotations in [1, "obs"]:
        GeneSynonyms()(
            adata.obs,
            axis="index",
            in_alias_type=in_alias_type,
            out_alias_type=out_alias_type,
            keep_if_missing=keep_if_missing,
            copy=False
        )
    else:
        raise ValueError("invalid value for 'annotations' (got {annotations}, expected 'obs' or 'var')")
    
    return adata if copy else None

def set_ncbi_reference_name(
    adata: ad.AnnData,
    annotations: str="var",
    in_alias_type: str="genename",
    keep_if_missing: bool=True,
    copy: bool = False
) -> Union[ad.AnnData, None]:
    """
    Replace gene names with theirs reference gene names into 'adata' object.

    Parameters
    ----------
    adata
        AnnData object where names are expected being standardized
    annotations
        whether to rename labels from adata.var (0 or 'var') or adata.obs (1 or 'obs')
    in_alias_type
        genename|geneid|ensemblid|<database>
    copy
        return a copy instead of updating 'adata' object
        
    Returns
    -------
    Depending on 'copy', update or return AnnData object with reference gene names.
    """

    adata = adata.copy() if copy else adata

    gene_synonyms_conversion(
        adata=adata,
        annotations=annotations,
        in_alias_type=in_alias_type,
        out_alias_type="referencename",
        keep_if_missing=keep_if_missing,
        copy=False
    )

    return adata if copy else None

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
