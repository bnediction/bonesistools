#!/usr/bin/env python

from typing import Optional, overload

from anndata import AnnData

from ..._compat import Literal
from ...databases.ncbi import genesyn as create_gene_synonyms
from ...databases.ncbi._typing import InputIdentifierType
from .._typing import (
    AnnDataAxisWithInteger,
    anndata_checker,
)
from .._validation import _as_anndata_axis


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: Literal[False] = False,
) -> None: ...


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    *,
    copy: Literal[True],
) -> AnnData: ...


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def mitochondrial_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
) -> Optional[AnnData]:  # type: ignore
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
    axis: {"obs", "var"} (default: "var")
        If `"obs"`, annotate `adata.obs`. If `"var"`, annotate `adata.var`.
        Deprecated values 0 and 1 are still accepted as aliases for `"obs"`
        and `"var"`.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with mitochondrial gene
        annotations added. Otherwise, updates `adata` in place and returns
        None.

        Mitochondrial gene annotations are stored in:

        - `adata.obs[key]`: Boolean annotation if `axis="obs"`;
        - `adata.var[key]`: Boolean annotation if `axis="var"`.

    """

    adata = adata.copy() if copy else adata

    genesyn = create_gene_synonyms()
    axis = _as_anndata_axis(axis, allow_integer=True)
    if axis == "obs":
        adata.obs[key] = False
    elif axis == "var":
        adata.var[key] = False

    mt_id = set()
    for k, v in genesyn.gene_aliases_mapping["name"].items():
        if k.startswith("mt-"):
            mt_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        gene_id = genesyn.get_gene_id(gene=index, input_identifier_type=index_type)
        if gene_id in mt_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: Literal[False] = False,
) -> None: ...


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    *,
    copy: Literal[True],
) -> AnnData: ...


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
) -> Optional[AnnData]: ...


@anndata_checker
def ribosomal_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
) -> Optional[AnnData]:  # type: ignore
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
    axis: {"obs", "var"} (default: "var")
        If `"obs"`, annotate `adata.obs`. If `"var"`, annotate `adata.var`.
        Deprecated values 0 and 1 are still accepted as aliases for `"obs"`
        and `"var"`.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with ribosomal gene
        annotations added. Otherwise, updates `adata` in place and returns
        None.

        Ribosomal gene annotations are stored in:

        - `adata.obs[key]`: Boolean annotation if `axis="obs"`;
        - `adata.var[key]`: Boolean annotation if `axis="var"`.

    """

    adata = adata.copy() if copy else adata

    genesyn = create_gene_synonyms()
    axis = _as_anndata_axis(axis, allow_integer=True)
    if axis == "obs":
        adata.obs[key] = False
    elif axis == "var":
        adata.var[key] = False

    rps_id = set()
    for k, v in genesyn.gene_aliases_mapping["name"].items():
        if k.startswith(("Rps", "Rpl", "Mrp")):
            rps_id.add(v.value.decode())

    for index in eval(f"adata.{axis}.index"):
        gene_id = genesyn.get_gene_id(gene=index, input_identifier_type=index_type)
        if gene_id in rps_id:
            exec(f"adata.{axis}.at['{index}','{key}'] = True")

    return adata if copy else None
