#!/usr/bin/env python

from typing import Any, Optional, cast

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt


def _expected_scale(
    matrix: np.ndarray,
    *,
    zero_center: bool = True,
    max_value: Optional[float] = None,
) -> np.ndarray:

    expected = matrix.astype(float, copy=True)
    means = expected.mean(axis=0)
    ddof = 1 if expected.shape[0] > 1 else 0
    stds = expected.std(axis=0, ddof=ddof)
    stds[(stds == 0) | ~np.isfinite(stds)] = 1
    if zero_center:
        expected = expected - means
    expected = expected / stds
    if max_value is not None:
        expected = np.clip(expected, -max_value, max_value)
    return expected


def test_scale_dense_x_in_place_adds_variable_statistics():

    values = np.array(
        [
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            [3.0, 6.0, 9.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    bt.sct.pp.scale(adata)

    expected = _expected_scale(values)
    np.testing.assert_allclose(np.asarray(cast(Any, adata.X)), expected)
    np.testing.assert_allclose(np.asarray(cast(Any, adata.X)).mean(axis=0), 0)
    np.testing.assert_allclose(np.asarray(cast(Any, adata.X)).std(axis=0, ddof=1), 1)
    np.testing.assert_allclose(adata.var["mean"].to_numpy(), values.mean(axis=0))
    np.testing.assert_allclose(adata.var["std"].to_numpy(), values.std(axis=0, ddof=1))


def test_scale_dense_expression_layer_in_place():

    values = np.array(
        [
            [1.0, 5.0],
            [2.0, 7.0],
            [4.0, 11.0],
        ]
    )
    adata = ad.AnnData(X=np.zeros(values.shape))
    adata.layers["log1p"] = values.copy()

    bt.sct.pp.scale(adata, expression="log1p")

    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["log1p"])),
        _expected_scale(values),
    )
    np.testing.assert_allclose(adata.var["mean"].to_numpy(), values.mean(axis=0))
    np.testing.assert_allclose(adata.var["std"].to_numpy(), values.std(axis=0, ddof=1))


def test_scale_dense_x_to_new_layer_preserves_source():

    values = np.array(
        [
            [1.0, 5.0],
            [2.0, 7.0],
            [4.0, 11.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    bt.sct.pp.scale(adata, key_added="scaled")

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.X)), values)
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["scaled"])),
        _expected_scale(values),
    )
    assert "scaled" not in adata.obsm


def test_scale_dense_expression_layer_to_new_layer_preserves_source():

    values = np.array(
        [
            [1.0, 5.0],
            [2.0, 7.0],
            [4.0, 11.0],
        ]
    )
    adata = ad.AnnData(X=np.zeros(values.shape))
    adata.layers["log1p"] = values.copy()

    bt.sct.pp.scale(adata, expression="log1p", key_added="scaled")

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.layers["log1p"])), values)
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["scaled"])),
        _expected_scale(values),
    )


def test_scale_dense_representation_in_place_does_not_update_var():

    values = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [4.0, 40.0],
        ]
    )
    adata = ad.AnnData(X=np.zeros((values.shape[0], 3)))
    adata.obsm["X_pca"] = values.copy()

    bt.sct.pp.scale(adata, representation="X_pca")

    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.obsm["X_pca"])),
        _expected_scale(values),
    )
    assert "mean" not in adata.var
    assert "std" not in adata.var


def test_scale_dense_representation_to_new_obsm_preserves_source():

    values = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [4.0, 40.0],
        ]
    )
    adata = ad.AnnData(X=np.zeros((values.shape[0], 3)))
    adata.obsm["X_pca"] = values.copy()

    bt.sct.pp.scale(
        adata,
        representation="X_pca",
        key_added="X_pca_scaled",
    )

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.obsm["X_pca"])), values)
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.obsm["X_pca_scaled"])),
        _expected_scale(values),
    )
    assert "X_pca_scaled" not in adata.layers


def test_scale_rejects_expression_and_representation_together():

    adata = ad.AnnData(X=np.ones((2, 2)))
    adata.layers["log1p"] = np.ones((2, 2))
    adata.obsm["X_pca"] = np.ones((2, 2))

    with pytest.raises(ValueError):
        bt.sct.pp.scale(
            adata,
            expression="log1p",
            representation="X_pca",
        )


def test_scale_zero_center_false_divides_by_standard_deviation_only():

    values = np.array(
        [
            [1.0, 2.0],
            [2.0, 4.0],
            [4.0, 8.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    bt.sct.pp.scale(adata, zero_center=False)

    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.X)),
        _expected_scale(values, zero_center=False),
    )


def test_scale_max_value_clips_scaled_values():

    values = np.array(
        [
            [0.0, 0.0],
            [0.0, 10.0],
            [0.0, 100.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    bt.sct.pp.scale(adata, max_value=1.0)

    scaled = np.asarray(cast(Any, adata.X))
    np.testing.assert_allclose(
        scaled,
        _expected_scale(values, max_value=1.0),
    )
    assert np.max(np.abs(scaled)) <= 1.0


def test_scale_constant_columns_do_not_create_nan_or_inf():

    values = np.array(
        [
            [1.0, 5.0, 2.0],
            [1.0, 7.0, 2.0],
            [1.0, 9.0, 2.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    bt.sct.pp.scale(adata)

    scaled = np.asarray(cast(Any, adata.X))
    expected = _expected_scale(values)
    np.testing.assert_allclose(scaled, expected)
    assert np.isfinite(scaled).all()
    np.testing.assert_array_equal(scaled[:, 0], np.zeros(values.shape[0]))
    np.testing.assert_array_equal(scaled[:, 2], np.zeros(values.shape[0]))
    np.testing.assert_allclose(
        adata.var["std"].to_numpy(),
        np.array([1.0, values[:, 1].std(ddof=1), 1.0]),
    )


def test_scale_sparse_zero_center_false_remains_sparse():

    values = np.array(
        [
            [0.0, 2.0, 0.0],
            [3.0, 0.0, 4.0],
            [0.0, 0.0, 0.0],
        ]
    )
    adata = ad.AnnData(X=csr_matrix(values))

    bt.sct.pp.scale(adata, zero_center=False)

    assert sparse.issparse(adata.X)
    np.testing.assert_allclose(
        cast(Any, adata.X).toarray(),
        _expected_scale(values, zero_center=False),
    )


def test_scale_sparse_zero_center_true_densifies_and_scales():

    values = np.array(
        [
            [0.0, 2.0, 0.0],
            [3.0, 0.0, 4.0],
            [0.0, 0.0, 0.0],
        ]
    )
    adata = ad.AnnData(X=csr_matrix(values))

    bt.sct.pp.scale(adata, zero_center=True)

    assert not sparse.issparse(adata.X)
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.X)),
        _expected_scale(values),
    )


def test_scale_copy_returns_new_anndata_and_leaves_original_unchanged():

    values = np.array(
        [
            [1.0, 2.0],
            [2.0, 4.0],
            [4.0, 8.0],
        ]
    )
    adata = ad.AnnData(X=values.copy())

    copied = bt.sct.pp.scale(adata, copy=True)

    assert copied is not adata
    np.testing.assert_array_equal(np.asarray(cast(Any, adata.X)), values)
    np.testing.assert_allclose(
        np.asarray(cast(Any, copied.X)),
        _expected_scale(values),
    )


def test_scale_rejects_invalid_argument_types(mini_adata):

    with pytest.raises(TypeError):
        bt.sct.pp.scale(mini_adata, expression=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.scale(mini_adata, representation=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.scale(mini_adata, key_added=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.scale(mini_adata, zero_center=cast(Any, "yes"))

    with pytest.raises(ValueError):
        bt.sct.pp.scale(mini_adata, max_value=0)
