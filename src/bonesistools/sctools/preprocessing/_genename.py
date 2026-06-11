#!/usr/bin/env python

from __future__ import annotations

import warnings
from typing import Any, Optional, Union, cast

import anndata as ad
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import csr_matrix

from ...databases.ncbi import GeneSynonyms
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
    genesyn = GeneSynonyms() if genesyn is None else genesyn

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


@anndata_checker
def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None,
    copy: bool = True,
) -> Union[AnnData, None]:
    """
    Merge duplicated variable names by summing counts across duplicate genes.

    The resulting AnnData object contains one variable per variable name. Counts
    from duplicated variables are summed across columns. If `var_names_column`
    is provided, that column is used to choose which metadata row is kept for
    each duplicated variable.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix containing duplicated variable names.
    var_names_column: str, optional
        Column used to prioritize which `adata.var` row is kept when duplicate
        variable names are merged.
    copy: bool (default: True)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        AnnData object with duplicated variable names merged if `copy=True`;
        otherwise None.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        source = adata.copy() if copy else adata

    if var_names_column is None:
        var_names = "copy_var_names"
        adata_var = cast(DataFrame, source.var)
        adata_var[var_names] = list(adata_var.index)
    else:
        var_names = var_names_column

    adata_obs = cast(DataFrame, source.obs)
    obs = DataFrame(adata_obs.index).set_index(0)
    adatas = list()
    duplicated_var_names = {
        str(source.var_names[idx])
        for idx, value in enumerate(source.var_names.duplicated())
        if value
    }

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )

        for var_name in duplicated_var_names:
            adata_spec = source[:, source.var.index == var_name]
            source = source[:, source.var.index != var_name]

            X = csr_matrix(cast(Any, adata_spec.X).sum(axis=1))
            adata_spec_var = cast(DataFrame, adata_spec.var)
            if var_name in list(adata_spec_var[var_names]):
                filter = adata_spec_var[var_names] == var_name
                var = adata_spec_var[filter].iloc[:1]
            else:
                var = adata_spec_var.iloc[:1]
            adata_spec = AnnData(X=X, var=var, obs=obs)
            adatas.append(adata_spec)

        adatas.append(source)

        merged = ad.concat(
            adatas=adatas,
            join="outer",
            axis=1,
            merge="same",
            uns_merge="first",
            label=None,
        )

    if var_names_column is None:
        adata_var = cast(DataFrame, merged.var)
        adata_var.drop(labels=var_names, axis="columns", inplace=True)

    if copy:
        return merged

    adata._init_as_actual(merged)
    return None
