#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Union, Optional

import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

from databases.genesyn import GeneSynonyms

def set_ncbi_reference_name(
    adata: ad.AnnData,
    annotations: str = "var",
    copy: bool = False
) -> Union[ad.AnnData, None]:

    adata = adata.copy() if copy else adata

    if annotations == "obs":
        GeneSynonyms()(adata.obs, axis="index", copy=False)
    elif annotations == "var":
        GeneSynonyms()(adata.var, axis="index", copy=False)
    else:
        raise ValueError("`annotations` must take either value `obs` or `var`")
    
    return adata if copy else None

def var_names_merge_duplicates(
    adata: ad.AnnData,
    var_names_column: Optional[str] = None
) -> Union[ad.AnnData, None]:

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
