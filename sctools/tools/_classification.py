#!/usr/bin/env python

from typing import Optional, Union, Sequence
from .._typing import (
    anndata_or_mudata_checker,
    ScData,
    Axis,
)

from ...databases.ncbi import GeneSynonyms

import pandas as pd

@anndata_or_mudata_checker
def mitochondrial_genes(
    scdata: ScData,  # type: ignore
    index_type: str = "genename",
    key: str = "mt",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[ScData, None]: # type: ignore
    """
    Create a column storing whether gene encodes a mitochondrial protein using current index.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    genesyn
        GeneSynonyms object
    index_type
        genename|geneid|ensemblid|<database>
    key
        column name storing whether gene encodes a mitochondrial protein
    axis
        whether to classify from scdata.var (0 or 'obs') or scdata.obs (1 or 'var')
    copy
        return a copy instead of updating 'scdata' object
    
    Returns
    -------
    Depending on 'copy', update or return AnnData or MuData object.
    """

    scdata = scdata.copy() if copy else scdata

    genesyn = GeneSynonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        scdata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        scdata.var[key] = False
    else:
        raise ValueError(f"invalid value for 'annotations' (got {axis}, expected 'obs' or 'var')")
    
    mt_id = set()
    for k,v in genesyn.gene_aliases_mapping["genename"].items():
        if k.startswith("mt-"):
            mt_id.add(v.value.decode())

    for index in eval(f"scdata.{axis}.index"):
        geneid = genesyn.get_geneid(
            alias=index,
            alias_type=index_type
        )
        if geneid in mt_id:
            exec(f"scdata.{axis}.at['{index}','{key}'] = True")

    return scdata if copy else None

@anndata_or_mudata_checker
def ribosomal_genes(
    scdata: ScData,  # type: ignore
    index_type: str = "genename",
    key: str = "rps",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[ScData, None]: # type: ignore
    """
    Create a column storing whether gene encodes a ribosomal protein using current index.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    genesyn
        GeneSynonyms object
    index_type
        genename|geneid|ensemblid|<database>
    key
        column name storing whether gene encodes a ribosomal protein
    axis
        whether to classify from scdata.var (0 or 'obs') or scdata.obs (1 or 'var')
    copy
        return a copy instead of updating 'scdata' object
    
    Returns
    -------
    Depending on 'copy', update or return AnnData or MuData object.
    """

    scdata = scdata.copy() if copy else scdata

    genesyn = GeneSynonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        scdata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        scdata.var[key] = False
    else:
        raise ValueError(f"invalid value for 'annotations' (got {axis}, expected 'obs' or 'var')")
    
    rps_id = set()
    for k,v in genesyn.gene_aliases_mapping["genename"].items():
        if k.startswith(("Rps","Rpl","Mrp")):
            rps_id.add(v.value.decode())

    for index in eval(f"scdata.{axis}.index"):
        geneid = genesyn.get_geneid(
            alias=index,
            alias_type=index_type
        )
        if geneid in rps_id:
            exec(f"scdata.{axis}.at['{index}','{key}'] = True")

    return scdata if copy else None
