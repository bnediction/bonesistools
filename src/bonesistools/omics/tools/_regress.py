#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Optional, Union, cast, overload

import numpy as np
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import issparse

from ..._compat import Literal
from ..._validation import _as_positive_integer
from .._typing import Keys, anndata_checker


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    *,
    copy: Literal[False],
    n_jobs: int = 1,
) -> None: ...


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    copy: Literal[True] = True,
    n_jobs: int = 1,
) -> AnnData: ...


@overload
def regress_out(
    adata: AnnData,
    keys: Keys,
    layer: Optional[str] = None,
    intercept: bool = False,
    copy: bool = False,
    n_jobs: int = 1,
) -> Union[AnnData, None]: ...


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
        If `True`, preserve the fitted intercept after regressing out covariates.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    n_jobs: int (default: 1)
        Number of BLAS/LAPACK threads allocated to the least-squares solve and
        matrix multiplication.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with corrected values added.
        Otherwise, updates `adata` in place and returns None.

        Corrected values are stored in:

        - `adata.X`: corrected matrix if `layer=None`;
        - `adata.layers[layer]`: corrected layer otherwise.
    """

    from threadpoolctl import threadpool_limits

    n_jobs = _as_positive_integer(n_jobs, "n_jobs")

    adata = adata.copy() if copy else adata

    if layer is None:
        counts = cast(Any, adata.X).copy()
    else:
        counts = cast(Any, adata.layers[layer]).copy()

    if issparse(counts):
        counts = cast(Any, counts).toarray()

    counts = np.asarray(counts)
    if not np.issubdtype(counts.dtype, np.floating):
        counts = counts.astype(np.float32)

    adata_obs = cast(DataFrame, adata.obs)
    regressors = adata_obs[[keys]] if isinstance(keys, str) else adata_obs[list(keys)]
    regressors = regressors.to_numpy()
    ones = np.ones((regressors.shape[0], 1), dtype=regressors.dtype)
    regressors = np.concatenate([ones, regressors], axis=1)

    with threadpool_limits(limits=n_jobs):
        coefficients = np.linalg.lstsq(regressors, counts, rcond=None)[0]
        prediction = regressors @ coefficients
    if intercept:
        counts = counts - prediction + coefficients[0, :]
    else:
        counts = counts - prediction

    if layer is None:
        adata.X = counts
    else:
        adata.layers[layer] = counts

    return adata if copy else None
