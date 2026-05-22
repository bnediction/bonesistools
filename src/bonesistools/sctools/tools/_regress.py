#!/usr/bin/env python

from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse

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
        If specified, regress values from `adata.layers[layer]` instead of
        `adata.X`.
    intercept: bool (default: False)
        If true, add intercept as regressor on which to regress on.
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.
    n_jobs: int (default: 1)
        Number of jobs for parallel computation.

    Returns
    -------
    Depending on 'copy', update 'adata' or return AnnData object.
    """

    from sklearn.linear_model import LinearRegression

    adata = adata.copy() if copy else adata

    if layer is None:
        counts = adata.X.copy()
    else:
        counts = adata.layers[layer].copy()

    if issparse(counts):
        counts = counts.toarray()
    regressors = adata.obs[[keys]] if isinstance(keys, str) else adata.obs[keys]
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
