#!/usr/bin/env python

from typing import Any, Sequence, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import bonesistools as bt


def _expected_regress_out(
    matrix: np.ndarray,
    obs: pd.DataFrame,
    keys: Sequence[str],
    intercept: bool = False,
) -> np.ndarray:

    regressors = obs[list(keys)].to_numpy()
    regressors = np.concatenate(
        [np.ones((regressors.shape[0], 1), dtype=regressors.dtype), regressors],
        axis=1,
    )
    coefficients = np.linalg.lstsq(regressors, matrix, rcond=None)[0]
    prediction = regressors @ coefficients
    if intercept:
        return matrix - prediction + coefficients[0, :]
    return matrix - prediction


def _make_regression_adata(
    sparse_input: bool = False,
    dtype: Any = np.float32,
) -> ad.AnnData:

    rng = np.random.default_rng(0)
    n_obs = 40
    n_vars = 25
    score = rng.normal(size=n_obs).astype(dtype)
    batch = rng.normal(size=n_obs).astype(dtype)
    obs = pd.DataFrame(
        {"score": score, "batch": batch},
        index=[f"c{i}" for i in range(n_obs)],
    )

    base = rng.normal(size=(n_obs, n_vars)).astype(dtype)
    weights = rng.normal(size=(2, n_vars)).astype(dtype)
    matrix = (
        base
        + score[:, None] * weights[0, :]
        + batch[:, None] * weights[1, :]
        + dtype(2.0)
    )

    if sparse_input:
        matrix = matrix.copy()
        matrix[rng.random(size=matrix.shape) < 0.5] = 0
        layer = sparse.csr_matrix(matrix)
    else:
        layer = matrix

    adata = ad.AnnData(
        X=np.zeros((n_obs, n_vars), dtype=dtype),
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_vars)]),
    )
    adata.layers["counts"] = layer
    return adata


def _dense_counts(adata: ad.AnnData) -> np.ndarray:

    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        return cast(Any, counts).toarray()
    return np.asarray(counts)


def test_regress_out_matches_vectorized_least_squares_without_intercept():
    adata = _make_regression_adata(dtype=np.float64)
    original = _dense_counts(adata).copy()
    expected = _expected_regress_out(
        original,
        cast(pd.DataFrame, adata.obs),
        keys=["score", "batch"],
        intercept=False,
    )

    bt.omics.tl.regress_out(
        adata,
        keys=["score", "batch"],
        layer="counts",
        intercept=False,
    )

    np.testing.assert_allclose(_dense_counts(adata), expected)


def test_regress_out_matches_vectorized_least_squares_with_intercept():
    adata = _make_regression_adata(dtype=np.float64)
    original = _dense_counts(adata).copy()
    expected = _expected_regress_out(
        original,
        cast(pd.DataFrame, adata.obs),
        keys=["score", "batch"],
        intercept=True,
    )

    bt.omics.tl.regress_out(
        adata,
        keys=["score", "batch"],
        layer="counts",
        intercept=True,
    )

    np.testing.assert_allclose(_dense_counts(adata), expected)


def test_regress_out_sparse_input_is_densified_and_matches_dense_result():
    sparse_adata = _make_regression_adata(sparse_input=True)
    dense_original = _dense_counts(sparse_adata).copy()
    expected = _expected_regress_out(
        dense_original,
        cast(pd.DataFrame, sparse_adata.obs),
        keys=["score", "batch"],
        intercept=True,
    )

    bt.omics.tl.regress_out(
        sparse_adata,
        keys=["score", "batch"],
        layer="counts",
        intercept=True,
    )

    assert isinstance(sparse_adata.layers["counts"], np.ndarray)
    np.testing.assert_allclose(_dense_counts(sparse_adata), expected, rtol=1e-6)


def test_regress_out_is_strictly_reproducible_on_identical_inputs():
    adata = _make_regression_adata(sparse_input=True)
    first = adata.copy()
    second = adata.copy()

    bt.omics.tl.regress_out(
        first,
        keys=["score", "batch"],
        layer="counts",
        intercept=True,
    )
    bt.omics.tl.regress_out(
        second,
        keys=["score", "batch"],
        layer="counts",
        intercept=True,
    )

    np.testing.assert_array_equal(_dense_counts(first), _dense_counts(second))


def test_regress_out_limits_threadpool_with_n_jobs(monkeypatch):
    import threadpoolctl

    calls = []

    class FakeThreadpoolLimits:
        def __init__(self, limits):
            self.limits = limits
            calls.append(("init", limits))

        def __enter__(self):
            calls.append(("enter", self.limits))
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            calls.append(("exit", self.limits))

    monkeypatch.setattr(threadpoolctl, "threadpool_limits", FakeThreadpoolLimits)

    adata = _make_regression_adata()
    bt.omics.tl.regress_out(
        adata,
        keys=["score", "batch"],
        layer="counts",
        n_jobs=3,
    )

    assert calls == [("init", 3), ("enter", 3), ("exit", 3)]


def test_regress_out_validates_n_jobs():
    adata = _make_regression_adata()

    with pytest.raises(ValueError, match="invalid argument value for 'n_jobs'"):
        bt.omics.tl.regress_out(
            adata,
            keys=["score", "batch"],
            layer="counts",
            n_jobs=0,
        )


def test_regress_out_copy_preserves_original_layer():
    adata = _make_regression_adata()
    original = _dense_counts(adata).copy()

    copied = bt.omics.tl.regress_out(
        adata,
        keys=["score", "batch"],
        layer="counts",
        intercept=True,
        copy=True,
    )

    np.testing.assert_array_equal(_dense_counts(adata), original)
    assert copied is not None
    assert not np.array_equal(_dense_counts(copied), original)
