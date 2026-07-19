#!/usr/bin/env python

from typing import Optional, overload

from anndata import AnnData

from ..._compat import Literal
from ..._warnings import _rename_deprecated_arguments, _warn_deprecated
from ...resources.ncbi import identifiers as create_identifiers
from ...resources.ncbi._identifiers import GeneIdentifiers
from ...resources.ncbi._typing import InputIdentifierType
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
def mitochondrial_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "mt",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> Optional[AnnData]:
    """
    Deprecated alias for `bt.omics.pp.mitochondrial_genes`.
    """

    from ..preprocessing import _classification as _preprocessing_classification

    _warn_deprecated(
        "`bt.omics.tl.mitochondrial_genes`",
        replacement="`bt.omics.pp.mitochondrial_genes`",
        stacklevel=2,
    )
    original_create_identifiers = _preprocessing_classification.create_identifiers
    try:
        _preprocessing_classification.create_identifiers = create_identifiers
        return _preprocessing_classification.mitochondrial_genes(
            adata,
            index_type=index_type,
            key=key,
            axis=axis,
            copy=copy,
            organism=organism,
            identifiers=identifiers,
        )
    finally:
        _preprocessing_classification.create_identifiers = original_create_identifiers


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
def ribosomal_genes(
    adata: AnnData,
    index_type: InputIdentifierType = "name",
    key: str = "rps",
    axis: AnnDataAxisWithInteger = "var",
    copy: bool = False,
    *,
    organism: str = "mouse",
    identifiers: Optional[GeneIdentifiers] = None,
) -> Optional[AnnData]:
    """
    Deprecated alias for `bt.omics.pp.ribosomal_genes`.
    """

    from ..preprocessing import _classification as _preprocessing_classification

    _warn_deprecated(
        "`bt.omics.tl.ribosomal_genes`",
        replacement="`bt.omics.pp.ribosomal_genes`",
        stacklevel=2,
    )
    original_create_identifiers = _preprocessing_classification.create_identifiers
    try:
        _preprocessing_classification.create_identifiers = create_identifiers
        return _preprocessing_classification.ribosomal_genes(
            adata,
            index_type=index_type,
            key=key,
            axis=axis,
            copy=copy,
            organism=organism,
            identifiers=identifiers,
        )
    finally:
        _preprocessing_classification.create_identifiers = original_create_identifiers
