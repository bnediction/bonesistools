#!/usr/bin/env python

from typing import Union
from anndata import AnnData
from .._typing import (
    anndata_checker,
    Axis,
)

from ...databases.ncbi import (
    GeneType,
    GeneSynonyms
)    

@anndata_checker
def mitochondrial_genes(
    adata: AnnData,  # type: ignore
    index_type: GeneType = "genename",
    key: str = "mt",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[AnnData, None]: # type: ignore
    """
    Create a column storing whether gene encodes a mitochondrial protein using current index.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    index_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
        Input identifier type of indices.
    key: str (default: 'mt')
        Column name storing whether gene encodes a mitochondrial protein.
    axis: Axis (default: 1)
        Whether to classify from adata.var (0 or 'obs') or adata.obs (1 or 'var').
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.
    
    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.
    """

    adata = adata.copy() if copy else adata

    genesyn = GeneSynonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        adata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        adata.var[key] = False
    else:
        raise TypeError(f"unsupported argument type for 'axis': {axis}")
    
    mt_id = set()
    for k,v in genesyn.gene_aliases_mapping["genename"].items():
        if k.startswith("mt-"):
            mt_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        geneid = genesyn.get_geneid(
            gene=index,
            gene_type=index_type
        )
        if geneid in mt_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None

@anndata_checker
def ribosomal_genes(
    adata: AnnData,  # type: ignore
    index_type: GeneType = "genename",
    key: str = "rps",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[AnnData, None]: # type: ignore
    """
    Create a column storing whether gene encodes a ribosomal protein using current index.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    index_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
        Input identifier type of indices.
    key: str (default: 'rps')
        Column name storing whether gene encodes a ribosomal protein.
    axis: Axis (default: 1)
        Whether to classify from adata.var (0 or 'obs') or adata.obs (1 or 'var').
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.
    
    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.
    """

    adata = adata.copy() if copy else adata

    genesyn = GeneSynonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        adata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        adata.var[key] = False
    else:
        raise TypeError(f"unsupported argument type for 'axis': {axis}")
    
    rps_id = set()
    for k,v in genesyn.gene_aliases_mapping["genename"].items():
        if k.startswith(("Rps","Rpl","Mrp")):
            rps_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        geneid = genesyn.get_geneid(
            gene=index,
            gene_type=index_type
        )
        if geneid in rps_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None
