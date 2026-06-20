#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt


def _scbolt_reference_normalize(
    adata: ad.AnnData,
    layer: str,
    target_sum: float = 1e4,
) -> None:

    matrix = adata.layers[layer]

    if sparse.issparse(matrix):
        matrix = matrix.tocsr(copy=True)
        if not np.issubdtype(matrix.dtype, np.floating):
            matrix = matrix.astype(np.float32)
        counts = np.asarray(matrix.sum(axis=1)).ravel()
        counts = counts / target_sum
        counts = counts + (counts == 0)
        matrix.data = np.true_divide(
            matrix.data,
            np.repeat(counts, np.diff(matrix.indptr)),
        )
        adata.layers[layer] = matrix
        return

    matrix = np.asarray(matrix)
    if not np.issubdtype(matrix.dtype, np.floating):
        matrix = matrix.astype(np.float32)
    counts = np.asarray(matrix.sum(axis=1)).ravel()
    counts = counts / target_sum
    counts = counts + (counts == 0)
    np.true_divide(matrix, counts[:, None], out=matrix)
    adata.layers[layer] = matrix


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)
    return cast(csr_matrix, matrix)


def _assert_same_sparse_matrix(left: Any, right: Any) -> None:

    left_csr = _as_csr(left)
    right_csr = _as_csr(right)

    np.testing.assert_array_equal(left_csr.indptr, right_csr.indptr)
    np.testing.assert_array_equal(left_csr.indices, right_csr.indices)
    np.testing.assert_array_equal(left_csr.data, right_csr.data)


def test_normalize_sparse_integer_x_remains_sparse_float32():

    counts = csr_matrix(
        np.array(
            [
                [1, 2, 0],
                [0, 0, 0],
                [3, 0, 1],
            ],
            dtype=np.int64,
        )
    )
    adata = ad.AnnData(X=counts)

    bt.sct.pp.normalize(adata, target_sum=10)

    normalized = _as_csr(adata.X)
    assert normalized.dtype == np.float32
    np.testing.assert_allclose(
        normalized.toarray(),
        np.array(
            [
                [10 / 3, 20 / 3, 0],
                [0, 0, 0],
                [7.5, 0, 2.5],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_dense_integer_x_uses_same_row_scaling():

    adata = ad.AnnData(
        X=np.array(
            [
                [1, 2, 0],
                [0, 0, 0],
                [3, 0, 1],
            ],
            dtype=np.int64,
        )
    )

    bt.sct.pp.normalize(adata, target_sum=10)

    normalized = np.asarray(cast(Any, adata.X))
    assert normalized.dtype == np.float32
    np.testing.assert_allclose(
        normalized,
        np.array(
            [
                [10 / 3, 20 / 3, 0],
                [0, 0, 0],
                [7.5, 0, 2.5],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_expression_overwrites_selected_layer_only():

    adata = ad.AnnData(X=np.ones((2, 2), dtype=np.float32))
    adata.layers["counts"] = csr_matrix(
        np.array(
            [
                [1, 3],
                [2, 0],
            ],
            dtype=np.int64,
        )
    )
    original_x = np.asarray(cast(Any, adata.X)).copy()

    bt.sct.pp.normalize(adata, target_sum=8, expression="counts")

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.X)), original_x)
    normalized = _as_csr(adata.layers["counts"])
    np.testing.assert_allclose(
        normalized.toarray(),
        np.array(
            [
                [2, 6],
                [8, 0],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_key_added_preserves_input_layer():

    counts = csr_matrix(
        np.array(
            [
                [1, 1],
                [2, 0],
            ],
            dtype=np.int64,
        )
    )
    adata = ad.AnnData(X=np.zeros((2, 2), dtype=np.float32))
    adata.layers["counts"] = counts.copy()

    bt.sct.pp.normalize(
        adata,
        target_sum=4,
        expression="counts",
        key_added="normalized",
    )

    _assert_same_sparse_matrix(adata.layers["counts"], counts)
    normalized = _as_csr(adata.layers["normalized"])
    np.testing.assert_allclose(
        normalized.toarray(),
        np.array(
            [
                [2, 2],
                [4, 0],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_key_added_preserves_x_when_expression_is_none():

    counts = np.array(
        [
            [1, 1],
            [2, 0],
        ],
        dtype=np.int64,
    )
    adata = ad.AnnData(X=counts.copy())

    bt.sct.pp.normalize(adata, target_sum=4, key_added="normalized")

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.X)), counts)
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["normalized"])),
        np.array(
            [
                [2, 2],
                [4, 0],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_target_sum_sets_nonzero_row_sums():

    adata = ad.AnnData(
        X=csr_matrix(
            np.array(
                [
                    [1, 2, 0],
                    [0, 0, 0],
                    [3, 0, 1],
                ],
                dtype=np.int64,
            )
        )
    )

    bt.sct.pp.normalize(adata, target_sum=25)

    row_sums = np.asarray(_as_csr(adata.X).sum(axis=1)).ravel()
    np.testing.assert_allclose(row_sums, np.array([25, 0, 25], dtype=np.float32))


def test_normalize_copy_returns_new_anndata_and_leaves_original_unchanged():

    counts = csr_matrix(
        np.array(
            [
                [1, 1],
                [0, 2],
            ],
            dtype=np.int64,
        )
    )
    adata = ad.AnnData(X=counts.copy())

    copied = bt.sct.pp.normalize(adata, target_sum=6, copy=True)

    assert copied is not adata
    _assert_same_sparse_matrix(adata.X, counts)
    normalized = _as_csr(copied.X)
    np.testing.assert_allclose(
        normalized.toarray(),
        np.array(
            [
                [3, 3],
                [0, 6],
            ],
            dtype=np.float32,
        ),
    )


def test_normalize_matches_scbolt_reference_helper_for_sparse_layer():

    counts = csr_matrix(
        np.array(
            [
                [0, 3, 1],
                [0, 0, 0],
                [5, 2, 3],
                [7, 0, 0],
            ],
            dtype=np.int64,
        )
    )
    reference = ad.AnnData(X=np.zeros(counts.shape, dtype=np.float32))
    reference.layers["counts"] = counts.copy()
    candidate = reference.copy()

    _scbolt_reference_normalize(reference, layer="counts", target_sum=1e4)
    bt.sct.pp.normalize(candidate, expression="counts", target_sum=1e4)

    _assert_same_sparse_matrix(candidate.layers["counts"], reference.layers["counts"])


def test_normalize_matches_scbolt_reference_helper_for_dense_layer():

    counts = np.array(
        [
            [0, 3, 1],
            [0, 0, 0],
            [5, 2, 3],
            [7, 0, 0],
        ],
        dtype=np.int64,
    )
    reference = ad.AnnData(X=np.zeros(counts.shape, dtype=np.float32))
    reference.layers["counts"] = counts.copy()
    candidate = reference.copy()

    _scbolt_reference_normalize(reference, layer="counts", target_sum=1e4)
    bt.sct.pp.normalize(candidate, expression="counts", target_sum=1e4)

    np.testing.assert_array_equal(
        np.asarray(cast(Any, candidate.layers["counts"])),
        np.asarray(cast(Any, reference.layers["counts"])),
    )


def test_normalize_rejects_invalid_arguments(mini_adata):

    with pytest.raises(ValueError):
        bt.sct.pp.normalize(mini_adata, target_sum=0)

    with pytest.raises(TypeError):
        bt.sct.pp.normalize(mini_adata, expression=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.normalize(mini_adata, key_added=cast(Any, object()))
