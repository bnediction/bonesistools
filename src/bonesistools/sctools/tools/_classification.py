#!/usr/bin/env python

from typing import Union

from anndata import AnnData

from ...databases.ncbi import genesyn as create_gene_synonyms
from ...databases.ncbi._typing import InputIdentifierType
from .._typing import (
    Axis,
    anndata_checker,
)


@anndata_checker
def mitochondrial_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[AnnData, None]:  # type: ignore
    """
    Annotate genes encoding mitochondrial proteins.

    The annotation is stored as a Boolean column in `adata.var` by default.
    Set `axis="obs"` to annotate observation names instead.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    index_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
        Input identifier type of the selected index.
    key: str (default: 'mt')
        Column name storing whether each gene encodes a mitochondrial protein.
    axis: {0, 1, "obs", "var"} (default: 1)
        If 0 or `"obs"`, annotate `adata.obs`. If 1 or `"var"`, annotate
        `adata.var`.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        Annotated AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    ValueError
        If `axis` is not 0, 1, `"obs"` or `"var"`.
    """

    adata = adata.copy() if copy else adata

    genesyn = create_gene_synonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        adata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        adata.var[key] = False
    else:
        raise ValueError(
            f"invalid argument value for 'axis': "
            f"expected 0, 1, 'obs' or 'var' but received {axis!r}"
        )

    mt_id = set()
    for k, v in genesyn.gene_aliases_mapping["name"].items():
        if k.startswith("mt-"):
            mt_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        gene_id = genesyn.get_gene_id(gene=index, input_identifier_type=index_type)
        if gene_id in mt_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None


@anndata_checker
def ribosomal_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: Axis = 1,
    copy: bool = False,
) -> Union[AnnData, None]:  # type: ignore
    """
    Annotate genes encoding ribosomal proteins.

    The annotation is stored as a Boolean column in `adata.var` by default.
    Set `axis="obs"` to annotate observation names instead.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    index_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
        Input identifier type of the selected index.
    key: str (default: 'rps')
        Column name storing whether each gene encodes a ribosomal protein.
    axis: {0, 1, "obs", "var"} (default: 1)
        If 0 or `"obs"`, annotate `adata.obs`. If 1 or `"var"`, annotate
        `adata.var`.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        Annotated AnnData object if `copy=True`; otherwise None.

    Raises
    ------
    ValueError
        If `axis` is not 0, 1, `"obs"` or `"var"`.
    """

    adata = adata.copy() if copy else adata

    genesyn = create_gene_synonyms()
    if axis in [0, "obs"]:
        axis = "obs"
        adata.obs[key] = False
    elif axis in [1, "var"]:
        axis = "var"
        adata.var[key] = False
    else:
        raise ValueError(
            f"invalid argument value for 'axis': "
            f"expected 0, 1, 'obs' or 'var' but received {axis!r}"
        )

    rps_id = set()
    for k, v in genesyn.gene_aliases_mapping["name"].items():
        if k.startswith(("Rps", "Rpl", "Mrp")):
            rps_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        gene_id = genesyn.get_gene_id(gene=index, input_identifier_type=index_type)
        if gene_id in rps_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None
