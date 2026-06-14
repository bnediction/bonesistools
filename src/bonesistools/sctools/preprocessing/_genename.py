#!/usr/bin/env python

from __future__ import annotations

from typing import Optional, Union, overload

from ..._compat import Literal
from ...databases.ncbi import genesyn as create_gene_synonyms
from ...databases.ncbi._genesyn import support_legacy_gene_synonyms_args
from ...databases.ncbi._typing import (
    GeneSynonymsLike,
    InputIdentifierType,
    OutputIdentifierType,
)
from .._typing import (
    Axis,
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
)
from .._validation import _as_anndata_axis


@overload
def convert_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: Literal[False] = False,
) -> None:
    ...


@overload
def convert_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: Literal[True] = True,
) -> ScData:
    ...


@overload
def convert_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: bool = False,
) -> Union[ScData, None]:
    ...


@support_legacy_gene_synonyms_args
@anndata_or_mudata_checker
def convert_gene_identifiers(
    scdata: ScData,  # type: ignore
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: bool = False,
) -> Union[ScData, None]:  # type: ignore
    """
    Convert gene identifiers stored in `scdata.obs` or `scdata.var`.

    Gene identifiers are converted through the bundled NCBI GeneSynonyms
    resource. By default, variable names are interpreted as gene names and
    converted to official gene symbols.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
        The stored gene identifiers are converted into the requested identifier
        type.
    axis: {"obs", "var"} (default: "var")
        If `"obs"`, convert `scdata.obs.index`. If `"var"`, convert
        `scdata.var.index`. Deprecated values 0 and 1 are still accepted as
        aliases for `"obs"` and `"var"`.
    input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
        (default: 'name')
        Input gene identifier type. Valid database-specific values are listed
        in `databases`.
    output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
        'ensembl_id' | <database> (default: 'official_name')
        Output gene identifier type. Valid database-specific values are listed
        in `databases`.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert gene identifiers. If None, create
        a default GeneSynonyms instance.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Returns
    -------
    AnnData, MuData or None
        If `copy=True`, returns a copy of `scdata` with gene identifiers
        converted. Otherwise, updates `scdata` in place and returns None.

        Converted identifiers are stored in:

        - `scdata.obs.index`: converted identifiers if `axis="obs"`;
        - `scdata.var.index`: converted identifiers if `axis="var"`.

    """

    scdata = scdata.copy() if copy else scdata
    genesyn = create_gene_synonyms() if genesyn is None else genesyn

    axis = _as_anndata_axis(axis, allow_integer=True)
    if axis == "obs":
        genesyn(
            scdata.obs,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False,
        )
    elif axis == "var":
        genesyn(
            scdata.var,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False,
        )
    return scdata if copy else None


@overload
def standardize_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: Literal[False] = False,
) -> None:
    ...


@overload
def standardize_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: Literal[True] = True,
) -> ScData:
    ...


@overload
def standardize_gene_identifiers(
    scdata: ScData,
    axis: Axis = "var",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: bool = False,
) -> Union[ScData, None]:
    ...


@support_legacy_gene_synonyms_args
@anndata_checker
def standardize_gene_identifiers(
    scdata: ScData,  # type: ignore
    axis: Axis = "var",
    genesyn: Optional[GeneSynonymsLike] = None,
    copy: bool = False,
) -> Union[ScData, None]:  # type: ignore
    """
    Standardize gene names by converting them to official names.

    This is a convenience wrapper around `convert_gene_identifiers` with
    `input_identifier_type="name"` and
    `output_identifier_type="official_name"`.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
        The stored gene names are standardized.
    axis: {"obs", "var"} (default: "var")
        If `"obs"`, standardize `scdata.obs.index`. If `"var"`, standardize
        `scdata.var.index`. Deprecated values 0 and 1 are still accepted as
        aliases for `"obs"` and `"var"`.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert gene identifiers. If None, create
        a default GeneSynonyms instance.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Returns
    -------
    AnnData, MuData or None
        If `copy=True`, returns a copy of `scdata` with gene identifiers
        standardized. Otherwise, updates `scdata` in place and returns None.

        Standardized identifiers are stored in:

        - `scdata.obs.index`: standardized identifiers if `axis="obs"`;
        - `scdata.var.index`: standardized identifiers if `axis="var"`.
    """

    return convert_gene_identifiers(
        scdata=scdata,
        axis=axis,
        input_identifier_type="name",
        output_identifier_type="official_name",
        genesyn=genesyn,
        copy=copy,
    )
