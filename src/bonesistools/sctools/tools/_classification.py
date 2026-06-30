#!/usr/bin/env python

from typing import Any, Optional, overload

from anndata import AnnData

from ..._compat import Literal
from ..._warnings import _warn_deprecated
from ...databases.ncbi import genesyn as create_gene_synonyms
from ...databases.ncbi._typing import InputIdentifierType
from .._typing import AnnDataAxisWithInteger


@overload
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: Literal[False] = False,
    *,
    organism: str = "mouse",
    genesyn: Optional[Any] = None,
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
    genesyn: Optional[Any] = None,
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
    genesyn: Optional[Any] = None,
) -> Optional[AnnData]: ...


def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    genesyn: Optional[Any] = None,
) -> Optional[AnnData]:
    """
    Deprecated alias for `bt.sct.pp.mitochondrial_genes`.
    """

    from ..preprocessing import _classification as _preprocessing_classification

    _warn_deprecated(
        "`bt.sct.tl.mitochondrial_genes`",
        replacement="`bt.sct.pp.mitochondrial_genes`",
        stacklevel=2,
    )
    original_create_gene_synonyms = _preprocessing_classification.create_gene_synonyms
    try:
        _preprocessing_classification.create_gene_synonyms = create_gene_synonyms
        return _preprocessing_classification.mitochondrial_genes(
            adata,
            index_type=index_type,
            key=key,
            axis=axis,
            copy=copy,
            organism=organism,
            genesyn=genesyn,
        )
    finally:
        _preprocessing_classification.create_gene_synonyms = (
            original_create_gene_synonyms
        )


@overload
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: Literal[False] = False,
    *,
    organism: str = "mouse",
    genesyn: Optional[Any] = None,
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
    genesyn: Optional[Any] = None,
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
    genesyn: Optional[Any] = None,
) -> Optional[AnnData]: ...


def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    genesyn: Optional[Any] = None,
) -> Optional[AnnData]:
    """
    Deprecated alias for `bt.sct.pp.ribosomal_genes`.
    """

    from ..preprocessing import _classification as _preprocessing_classification

    _warn_deprecated(
        "`bt.sct.tl.ribosomal_genes`",
        replacement="`bt.sct.pp.ribosomal_genes`",
        stacklevel=2,
    )
    original_create_gene_synonyms = _preprocessing_classification.create_gene_synonyms
    try:
        _preprocessing_classification.create_gene_synonyms = create_gene_synonyms
        return _preprocessing_classification.ribosomal_genes(
            adata,
            index_type=index_type,
            key=key,
            axis=axis,
            copy=copy,
            organism=organism,
            genesyn=genesyn,
        )
    finally:
        _preprocessing_classification.create_gene_synonyms = (
            original_create_gene_synonyms
        )
