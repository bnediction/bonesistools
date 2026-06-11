#!/usr/bin/env python

from __future__ import annotations

from typing import Optional, Union

from ...databases.ncbi import GeneSynonyms
from ...databases.ncbi import genesyn as create_gene_synonyms
from ...databases.ncbi._genesyn import support_legacy_gene_synonyms_args
from ...databases.ncbi._typing import (
    InputIdentifierType,
    OutputIdentifierType,
)
from .._typing import (
    Axis,
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
)


@support_legacy_gene_synonyms_args
@anndata_or_mudata_checker
def convert_gene_identifiers(
    scdata: ScData,  # type: ignore
    axis: Axis = "var",
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    genesyn: Optional[GeneSynonyms] = None,
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
    axis: {0, 1, "obs", "var"} (default: "var")
        If 0 or `"obs"`, convert `scdata.obs.index`. If 1 or `"var"`,
        convert `scdata.var.index`.
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
        Converted object if `copy=True`; otherwise None.

    Raises
    ------
    ValueError
        If `axis` is not 0, 1, `"obs"` or `"var"`.
    """

    scdata = scdata.copy() if copy else scdata
    genesyn = create_gene_synonyms() if genesyn is None else genesyn

    if axis in [0, "obs"]:
        genesyn(
            scdata.obs,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False,
        )
    elif axis in [1, "var"]:
        genesyn(
            scdata.var,
            axis="index",
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False,
        )
    else:
        raise ValueError(
            f"invalid argument value for 'axis': "
            f"expected 0, 1, 'obs' or 'var' but received {axis!r}"
        )

    return scdata if copy else None


@support_legacy_gene_synonyms_args
@anndata_checker
def standardize_gene_identifiers(
    scdata: ScData,  # type: ignore
    axis: Axis = "var",
    genesyn: Optional[GeneSynonyms] = None,
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
    axis: {0, 1, "obs", "var"} (default: "var")
        If 0 or `"obs"`, standardize `scdata.obs.index`. If 1 or `"var"`,
        standardize `scdata.var.index`.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert gene identifiers. If None, create
        a default GeneSynonyms instance.
    copy: bool (default: False)
        Return a copy instead of modifying `scdata`.

    Returns
    -------
    AnnData, MuData or None
        Standardized object if `copy=True`; otherwise None.
    """

    return convert_gene_identifiers(
        scdata=scdata,
        axis=axis,
        input_identifier_type="name",
        output_identifier_type="official_name",
        genesyn=genesyn,
        copy=copy,
    )
