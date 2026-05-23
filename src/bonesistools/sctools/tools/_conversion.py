#!/usr/bin/env python

from typing import (
    Optional,
    Sequence,
    Union,
)
from .._typing import anndata_checker, Keys

from anndata import AnnData
from pandas import DataFrame

import numpy as np
from scipy.sparse import issparse


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

    if layer:
        if issparse(adata.layers[layer]):
            counts_df = DataFrame(
                adata.layers[layer].toarray(),
                index=adata.obs.index,
                columns=adata.var.index,
            )
        else:
            counts_df = DataFrame(
                adata.layers[layer], index=adata.obs.index, columns=adata.var.index
            )
    elif issparse(adata.X):
        counts_df = DataFrame(
            adata.X.toarray(), index=adata.obs.index, columns=adata.var.index
        )
    else:
        counts_df = DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)

    if is_log:
        if "log1p" in adata.uns_keys() and adata.uns["log1p"].get("base") is not None:
            counts_df = np.expm1(counts_df * np.log(adata.uns["log1p"]["base"]))
        else:
            counts_df = np.expm1(counts_df)

    if obs is not None:
        counts_df.loc[:, obs] = (
            adata.obs.loc[:, [obs]] if isinstance(obs, str) else adata.obs.loc[:, obs]
        )

    return counts_df
