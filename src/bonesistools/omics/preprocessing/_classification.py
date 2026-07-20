#!/usr/bin/env python

from typing import AbstractSet, Optional, cast, overload

import pandas as pd
from anndata import AnnData

from ..._compat import Literal
from ..._warnings import _rename_deprecated_arguments
from ...resources.ncbi import identifiers as create_identifiers
from ...resources.ncbi._identifiers import GeneIdentifiers
from ...resources.ncbi._typing import InputIdentifierType
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
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> None: ...


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    *,
    copy: Literal[True],
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> AnnData: ...


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> Optional[AnnData]: ...


@_rename_deprecated_arguments(genesyn="identifiers")
@anndata_checker
def mitochondrial_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
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
    organism: str (default: "mouse")
        Organism used to resolve gene identifiers and mitochondrial aliases.
    identifiers: GeneIdentifiers, optional
        Existing NCBI gene information converter. If provided, `organism` is
        ignored. If `None`, a converter is created for `organism`.

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

    axis = _as_anndata_axis(axis, allow_integer=True)
    adata = adata.copy() if copy else adata
    gene_identifiers = (
        create_identifiers(organism=organism) if identifiers is None else identifiers
    )

    mt_id = {
        gene_id
        for gene_id, gene in gene_identifiers._iter_gene_identifiers()
        if gene.chromosome == "MT"
    }

    annotations = cast(pd.DataFrame, adata.obs if axis == "obs" else adata.var)
    _annotate_gene_id_set(
        annotations,
        key=key,
        gene_identifiers=gene_identifiers,
        index_type=index_type,
        target_gene_ids=mt_id,
    )

    return adata if copy else None


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: Literal[False] = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> None: ...


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    *,
    copy: Literal[True],
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> AnnData: ...


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> Optional[AnnData]: ...


@_rename_deprecated_arguments(genesyn="identifiers")
@anndata_checker
def ribosomal_genes(
    adata: AnnData,  # type: ignore
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
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
    organism: str (default: "mouse")
        Organism used to resolve gene identifiers and ribosomal official
        symbols.
    identifiers: GeneIdentifiers, optional
        Existing NCBI gene information converter. If provided, `organism` is
        ignored. If `None`, a converter is created for `organism`.

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

    axis = _as_anndata_axis(axis, allow_integer=True)
    adata = adata.copy() if copy else adata
    gene_identifiers = (
        create_identifiers(organism=organism) if identifiers is None else identifiers
    )

    rps_id = {
        gene_id
        for gene_id, gene in gene_identifiers._iter_gene_identifiers()
        if getattr(gene, "symbol", "")
        .lower()
        .startswith(("rps", "rpl", "mrps", "mrpl"))
    }

    annotations = cast(pd.DataFrame, adata.obs if axis == "obs" else adata.var)
    _annotate_gene_id_set(
        annotations,
        key=key,
        gene_identifiers=gene_identifiers,
        index_type=index_type,
        target_gene_ids=rps_id,
    )

    return adata if copy else None


def _annotate_gene_id_set(
    annotations: pd.DataFrame,
    *,
    key: str,
    gene_identifiers: GeneIdentifiers,
    index_type: InputIdentifierType,
    target_gene_ids: AbstractSet[str],
) -> None:
    """Annotate whether index entries resolve to selected NCBI gene IDs."""

    convert_sequence = getattr(gene_identifiers, "convert_sequence", None)
    if convert_sequence is None:
        gene_ids = [
            gene_identifiers.get_gene_id(gene=gene, input_type=index_type)
            for gene in annotations.index
        ]
    else:
        gene_ids = convert_sequence(
            annotations.index.tolist(),
            input_type=index_type,
            output_type="gene_id",
            keep_if_missing=False,
        )
    annotations[key] = [gene_id in target_gene_ids for gene_id in gene_ids]
