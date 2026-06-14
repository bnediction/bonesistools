#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Optional, Union, cast, overload

import numpy as np
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import issparse

from ..._compat import Literal
from .._dependencies import require_sklearn
from .._typing import Keys, anndata_checker


def __linear_regress_out_feature(
    interest: np.ndarray,
    regressors: np.ndarray,
    linear_regression,
    intercept: bool = False,
    n_jobs: int = 1,
):
    regression_model = linear_regression(fit_intercept=False, n_jobs=n_jobs)
    regression_model.fit(regressors, interest)
    prediction = regression_model.predict(regressors)

    if intercept:
        intercept = regression_model.coef_[0][0]
        predicted = interest - prediction + intercept
    else:
        predicted = interest - prediction

    return predicted[:, 0]


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    *,
    copy: Literal[False],
    n_jobs: int = 1,
) -> None:
    ...


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    copy: Literal[True] = True,
    n_jobs: int = 1,
) -> AnnData:
    ...


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    copy: bool = False,
    n_jobs: int = 1,
) -> Union[AnnData, None]:
    ...


@require_sklearn
@anndata_checker
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    copy: bool = False,
    n_jobs: int = 1,
) -> Union[AnnData, None]:
    """
    Regress out unwanted sources of variation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    keys: Keys
        Keys in `adata.obs` used as regressors.
    layer: str, optional
        Layer to use instead of `adata.X`.
    intercept: bool (default: False)
        If True, preserve the fitted intercept after regressing out covariates.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    n_jobs: int (default: 1)
        Number of allocated processors.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with corrected values added.
        Otherwise, updates `adata` in place and returns None.

        Corrected values are stored in:

        - `adata.X`: corrected matrix if `layer=None`;
        - `adata.layers[layer]`: corrected layer otherwise.
    """

    from sklearn.linear_model import LinearRegression

    adata = adata.copy() if copy else adata

    if layer is None:
        counts = cast(Any, adata.X).copy()
    else:
        counts = cast(Any, adata.layers[layer]).copy()

    if issparse(counts):
        counts = cast(Any, counts).toarray()

    counts = np.asarray(counts)
    adata_obs = cast(DataFrame, adata.obs)
    regressors = adata_obs[[keys]] if isinstance(keys, str) else adata_obs[keys]
    regressors.insert(0, "ones", 1.0)
    regressors = regressors.to_numpy()

    for i in range(adata.n_vars):
        interest = counts[:, i].reshape(-1, 1)
        counts[:, i] = __linear_regress_out_feature(
            interest,
            regressors,
            LinearRegression,
            intercept=intercept,
            n_jobs=n_jobs,
        )

    if layer is None:
        adata.X = counts
    else:
        adata.layers[layer] = counts

    return adata if copy else None
