#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.preprocessing import _qc


def _qc_adata(sparse: bool = False) -> ad.AnnData:

    X = np.array(
        [
            [1.0, 0.0, 3.0, 2.0],
            [0.0, 5.0, 0.0, 1.0],
            [2.0, 2.0, 2.0, 2.0],
        ],
        dtype=np.float32,
    )
    return ad.AnnData(
        X=csr_matrix(X) if sparse else X,
        obs=pd.DataFrame(index=["c1", "c2", "c3"]),
        var=pd.DataFrame(
            {"mito": [True, False, True, False]},
            index=["g1", "g2", "g3", "g4"],
        ),
    )


def _assert_expected_qc_metrics(adata: ad.AnnData) -> None:

    obs = cast(pd.DataFrame, adata.obs)
    var = cast(pd.DataFrame, adata.var)
    expected_obs = pd.DataFrame(
        {
            "n_features": [3, 2, 4],
            "log1p_n_features": np.log1p([3, 2, 4]),
            "total": np.array([6.0, 6.0, 8.0], dtype=np.float32),
            "log1p_total": np.log1p(np.array([6.0, 6.0, 8.0], dtype=np.float32)),
            "pct_top1_features": [50.0, 500.0 / 6.0, 25.0],
            "pct_top2_features": [500.0 / 6.0, 100.0, 50.0],
            "pct_top4_features": [100.0, 100.0, 100.0],
            "total_mito": np.array([4.0, 0.0, 4.0], dtype=np.float32),
            "log1p_total_mito": np.log1p(np.array([4.0, 0.0, 4.0], dtype=np.float32)),
            "pct_mito": np.array([200.0 / 3.0, 0.0, 50.0], dtype=np.float32),
        },
        index=["c1", "c2", "c3"],
    )
    expected_var = pd.DataFrame(
        {
            "n_barcodes": [2, 2, 2, 3],
            "mean": np.array([1.0, 7.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0]),
            "log1p_mean": np.log1p(np.array([1.0, 7.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0])),
            "variance": np.array([1.0, 19.0 / 3.0, 7.0 / 3.0, 1.0 / 3.0]),
            "log1p_variance": np.log1p(
                np.array([1.0, 19.0 / 3.0, 7.0 / 3.0, 1.0 / 3.0])
            ),
            "median": np.array([1.0, 2.0, 2.0, 2.0]),
            "log1p_median": np.log1p(np.array([1.0, 2.0, 2.0, 2.0])),
            "mad": np.array([1.0, 2.0, 1.0, 0.0]),
            "log1p_mad": np.log1p(np.array([1.0, 2.0, 1.0, 0.0])),
            "pct_dropout": [100.0 / 3.0, 100.0 / 3.0, 100.0 / 3.0, 0.0],
            "total": np.array([3.0, 7.0, 5.0, 5.0], dtype=np.float32),
            "log1p_total": np.log1p(np.array([3.0, 7.0, 5.0, 5.0], dtype=np.float32)),
        },
        index=["g1", "g2", "g3", "g4"],
    )

    pd.testing.assert_frame_equal(
        obs.loc[:, expected_obs.columns],
        expected_obs,
        check_dtype=False,
        rtol=1e-6,
        atol=1e-6,
    )
    pd.testing.assert_frame_equal(
        var.loc[:, expected_var.columns],
        expected_var,
        check_dtype=False,
        rtol=1e-6,
        atol=1e-6,
    )


@pytest.mark.parametrize("sparse", [False, True])
def test_qc_matches_expected_metrics_for_dense_and_sparse_inputs(sparse):
    adata = _qc_adata(sparse=sparse)

    result = bt.sct.pp.qc(
        adata,
        qc_vars=["mito"],
        percent_top=[1, 2, 4],
    )

    assert result is None
    _assert_expected_qc_metrics(adata)


def test_qc_auto_backend_uses_python_fallback_when_numba_is_unavailable(monkeypatch):
    adata = _qc_adata(sparse=True)
    monkeypatch.setattr(_qc, "_median_mad_sparse_csc_numba", None)
    monkeypatch.setattr(_qc, "_percent_top_sparse_csr_numba", None)

    bt.sct.pp.qc(
        adata,
        qc_vars=["mito"],
        percent_top=[1, 2, 4],
        backend="auto",
    )

    _assert_expected_qc_metrics(adata)


def test_qc_numba_backend_requires_numba(monkeypatch):
    adata = _qc_adata(sparse=True)
    monkeypatch.setattr(_qc, "_numba_njit", None)

    with pytest.raises(ImportError, match="backend='numba' requires"):
        bt.sct.pp.qc(adata, backend="numba")


def test_qc_python_and_numba_backends_agree_within_tolerance():
    if _qc._numba_njit is None:
        pytest.skip("Numba is not installed")

    python_adata = _qc_adata(sparse=True)
    numba_adata = _qc_adata(sparse=True)

    bt.sct.pp.qc(
        python_adata,
        qc_vars=["mito"],
        percent_top=[1, 2, 4],
        backend="python",
    )
    bt.sct.pp.qc(
        numba_adata,
        qc_vars=["mito"],
        percent_top=[1, 2, 4],
        backend="numba",
    )

    pd.testing.assert_frame_equal(
        cast(pd.DataFrame, numba_adata.obs),
        cast(pd.DataFrame, python_adata.obs),
        check_dtype=False,
        rtol=1e-5,
        atol=1e-8,
    )
    pd.testing.assert_frame_equal(
        cast(pd.DataFrame, numba_adata.var),
        cast(pd.DataFrame, python_adata.var),
        check_dtype=False,
        rtol=1e-5,
        atol=1e-8,
    )


def test_qc_uses_expression_layer_and_can_disable_log_and_top_metrics():
    adata = _qc_adata()
    adata.layers["doubled"] = cast(Any, adata.X).copy() * 2

    bt.sct.pp.qc(
        adata,
        expression="doubled",
        qc_vars="mito",
        percent_top=None,
        log1p=False,
    )

    assert adata.obs["total"].tolist() == [12.0, 12.0, 16.0]
    assert adata.obs["total_mito"].tolist() == [8.0, 0.0, 8.0]
    assert adata.var["total"].tolist() == [6.0, 14.0, 10.0, 10.0]
    assert "log1p_total" not in adata.obs
    assert "pct_top1_features" not in adata.obs


def test_qc_copy_returns_modified_copy_without_touching_input():
    adata = _qc_adata()

    copied = bt.sct.pp.qc(
        adata,
        qc_vars=["mito"],
        percent_top=[1],
        copy=True,
    )

    assert "total" not in adata.obs
    assert "total" in copied.obs
    assert "total" in copied.var


def test_qc_handles_zero_count_observations():
    adata = ad.AnnData(
        X=np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["zero", "nonzero"]),
        var=pd.DataFrame({"mito": [True, False]}, index=["g1", "g2"]),
    )

    bt.sct.pp.qc(adata, qc_vars=["mito"], percent_top=[1])

    obs = cast(pd.DataFrame, adata.obs)
    assert np.isnan(cast(float, obs.loc["zero", "pct_top1_features"]))
    assert np.isnan(cast(float, obs.loc["zero", "pct_mito"]))
    assert cast(float, obs.loc["nonzero", "pct_top1_features"]) == 100.0
    assert cast(float, obs.loc["nonzero", "pct_mito"]) == 100.0


def test_qc_sparse_median_and_mad_include_implicit_zeros():
    adata = ad.AnnData(
        X=csr_matrix(
            np.array(
                [
                    [0.0, 0.0, 10.0],
                    [0.0, 4.0, 0.0],
                    [6.0, 0.0, 0.0],
                    [8.0, 0.0, 0.0],
                ],
                dtype=np.float32,
            )
        ),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )

    bt.sct.pp.qc(adata, percent_top=None)

    np.testing.assert_allclose(adata.var["median"].to_numpy(), [3.0, 0.0, 0.0])
    np.testing.assert_allclose(adata.var["mad"].to_numpy(), [3.0, 0.0, 0.0])


def test_qc_sparse_median_and_mad_numba_path_matches_python_fallback(monkeypatch):
    if _qc._median_mad_sparse_csc_numba is None:
        pytest.skip("Numba is not installed")

    matrix = csr_matrix(
        np.array(
            [
                [0.0, 0.0, 10.0, 0.0],
                [0.0, 4.0, 0.0, 2.0],
                [6.0, 0.0, 0.0, 0.0],
                [8.0, 0.0, 0.0, 2.0],
                [0.0, 0.0, 0.0, 5.0],
            ],
            dtype=np.float32,
        )
    )

    monkeypatch.setattr(_qc, "_median_mad_sparse_csc_numba", None)
    fallback_median, fallback_mad = _qc._median_mad_per_feature(matrix)
    monkeypatch.undo()

    numba_median, numba_mad = _qc._median_mad_per_feature(matrix)

    np.testing.assert_array_equal(numba_median, fallback_median)
    np.testing.assert_array_equal(numba_mad, fallback_mad)


def test_qc_validates_arguments():
    adata = _qc_adata()

    with pytest.raises(IndexError, match="Positions outside range of features"):
        bt.sct.pp.qc(adata, percent_top=[5])

    with pytest.raises(TypeError, match="unsupported argument type for 'qc_vars'"):
        bt.sct.pp.qc(adata, qc_vars=cast(Any, 1), percent_top=None)

    with pytest.raises(KeyError):
        bt.sct.pp.qc(adata, qc_vars=["missing"], percent_top=None)

    adata.var["bad"] = ["yes", "no", "yes", "no"]
    with pytest.raises(TypeError, match="unsupported column dtype"):
        bt.sct.pp.qc(adata, qc_vars=["bad"], percent_top=None)

    with pytest.raises(ValueError, match="unsupported backend"):
        bt.sct.pp.qc(adata, backend=cast(Any, "bad"))
