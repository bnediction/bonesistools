#!/usr/bin/env python

from __future__ import annotations

from typing import (
    Any,
    Optional,
    cast,
)

import numpy as np
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import issparse

from .._typing import Keys, anndata_checker


@anndata_checker
def anndata_to_dataframe(
    adata: AnnData,
    obs: Optional[Keys] = None,
    layer: Optional[str] = None,
    is_log: bool = False,
) -> DataFrame:
    """
    Convert an AnnData expression matrix to a DataFrame.

    The returned DataFrame is indexed by observations and uses variable names
    as expression columns. Observation annotations can optionally be appended.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    obs: str or sequence of str, optional
        Observation columns from `adata.obs` to append to the output DataFrame.
    layer: str, optional
        Layer to use instead of `adata.X`.
    is_log: bool (default: False)
        If True, back-transform log1p values with `expm1`.

    Returns
    -------
    DataFrame
        Expression values with optional observation annotations.
    """

    matrix = adata.layers[layer] if layer else adata.X

    if issparse(matrix):
        matrix = cast(Any, matrix).toarray()

    matrix = cast(np.ndarray, matrix)

    adata_obs = cast(DataFrame, adata.obs)
    adata_var = cast(DataFrame, adata.var)

    counts_df = DataFrame(matrix, index=adata_obs.index, columns=adata_var.index)

    if is_log:
        if "log1p" in adata.uns and adata.uns["log1p"].get("base") is not None:
            matrix = np.expm1(counts_df * np.log(adata.uns["log1p"]["base"]))
        else:
            matrix = np.expm1(counts_df)

        matrix = cast(np.ndarray, matrix)
        counts_df = DataFrame(matrix, index=adata_obs.index, columns=adata_var.index)

    if obs is not None:
        if isinstance(obs, str):
            obs_df = adata_obs.loc[:, [obs]]
        else:
            obs_df = adata_obs.loc[:, list(obs)]

        counts_df = counts_df.join(obs_df)

    return counts_df
