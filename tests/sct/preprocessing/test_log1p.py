#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)
    return cast(csr_matrix, matrix)


def _assert_same_sparse_matrix(left: Any, right: Any) -> None:

    left_csr = _as_csr(left)
    right_csr = _as_csr(right)

    np.testing.assert_array_equal(left_csr.indptr, right_csr.indptr)
    np.testing.assert_array_equal(left_csr.indices, right_csr.indices)
    np.testing.assert_array_equal(left_csr.data, right_csr.data)


def test_log1p_dense_x_in_place():

    counts = np.array(
        [
            [0.0, 1.0, 3.0],
            [4.0, 0.0, 8.0],
        ],
        dtype=np.float64,
    )
    adata = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(adata)

    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.X)),
        np.log1p(counts),
    )


def test_log1p_dense_expression_layer_in_place():

    values = np.array(
        [
            [1.0, 3.0],
            [0.0, 8.0],
        ],
        dtype=np.float64,
    )
    adata = ad.AnnData(X=np.zeros_like(values))
    adata.layers["normalized"] = values.copy()

    bt.sct.pp.log1p(adata, expression="normalized")

    np.testing.assert_array_equal(np.asarray(cast(Any, adata.X)), np.zeros_like(values))
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["normalized"])),
        np.log1p(values),
    )


def test_log1p_key_added_preserves_source_layer():

    values = np.array(
        [
            [1.0, 3.0],
            [0.0, 8.0],
        ],
        dtype=np.float64,
    )
    adata = ad.AnnData(X=np.zeros_like(values))
    adata.layers["normalized"] = values.copy()

    bt.sct.pp.log1p(adata, expression="normalized", key_added="log1p")

    np.testing.assert_array_equal(
        np.asarray(cast(Any, adata.layers["normalized"])),
        values,
    )
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata.layers["log1p"])),
        np.log1p(values),
    )


def test_log1p_sparse_input_remains_sparse_and_transforms_nonzero_values():

    counts = csr_matrix(
        np.array(
            [
                [0, 1, 3],
                [4, 0, 8],
            ],
            dtype=np.int64,
        )
    )
    adata = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(adata)

    transformed = _as_csr(adata.X)
    assert transformed.dtype == np.float32
    np.testing.assert_array_equal(transformed.indptr, counts.indptr)
    np.testing.assert_array_equal(transformed.indices, counts.indices)
    np.testing.assert_allclose(transformed.data, np.log1p(counts.data.astype(float)))
    np.testing.assert_allclose(transformed.toarray(), np.log1p(counts.toarray()))


def test_log1p_dense_integer_input_returns_float_values():

    counts = np.array(
        [
            [0, 1],
            [3, 8],
        ],
        dtype=np.int64,
    )
    adata = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(adata)

    transformed = np.asarray(cast(Any, adata.X))
    assert transformed.dtype == np.float32
    np.testing.assert_allclose(transformed, np.log1p(counts))


def test_log1p_base_two_and_ten():

    counts = np.array(
        [
            [0.0, 1.0],
            [3.0, 8.0],
        ],
        dtype=np.float64,
    )
    adata_base2 = ad.AnnData(X=counts.copy())
    adata_base10 = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(adata_base2, base=2)
    bt.sct.pp.log1p(adata_base10, base=10)

    np.testing.assert_allclose(
        np.asarray(cast(Any, adata_base2.X)),
        np.log1p(counts) / np.log(2),
    )
    np.testing.assert_allclose(
        np.asarray(cast(Any, adata_base10.X)),
        np.log1p(counts) / np.log(10),
    )


def test_log1p_copy_returns_new_anndata_and_leaves_original_unchanged():

    counts = csr_matrix(
        np.array(
            [
                [0, 1],
                [3, 8],
            ],
            dtype=np.int64,
        )
    )
    adata = ad.AnnData(X=counts.copy())

    copied = bt.sct.pp.log1p(adata, copy=True)

    assert copied is not adata
    _assert_same_sparse_matrix(adata.X, counts)
    np.testing.assert_allclose(_as_csr(copied.X).toarray(), np.log1p(counts.toarray()))


def test_log1p_max_memory_chunking_preserves_dense_results():

    counts = np.array(
        [
            [0, 1, 3, 8],
            [2, 4, 6, 10],
            [0, 0, 5, 7],
        ],
        dtype=np.int64,
    )
    unchunked = ad.AnnData(X=counts.copy())
    chunked = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(unchunked, max_memory=None)
    bt.sct.pp.log1p(chunked, max_memory=16)

    np.testing.assert_array_equal(
        np.asarray(cast(Any, chunked.X)),
        np.asarray(cast(Any, unchunked.X)),
    )


def test_log1p_zero_entries_remain_zero():

    counts = np.array(
        [
            [0.0, 0.0],
            [0.0, 3.0],
        ]
    )
    adata = ad.AnnData(X=counts.copy())

    bt.sct.pp.log1p(adata)

    transformed = np.asarray(cast(Any, adata.X))
    assert transformed[0, 0] == 0
    assert transformed[0, 1] == 0
    assert transformed[1, 0] == 0
    assert transformed[1, 1] == np.log1p(3)


def test_log1p_rejects_invalid_base_values(mini_adata):

    with pytest.raises(ValueError):
        bt.sct.pp.log1p(mini_adata, base=0)

    with pytest.raises(ValueError):
        bt.sct.pp.log1p(mini_adata, base=-2)

    with pytest.raises(ValueError):
        bt.sct.pp.log1p(mini_adata, base=1)


def test_log1p_rejects_invalid_expression_key_added_and_memory(mini_adata):

    with pytest.raises(TypeError):
        bt.sct.pp.log1p(mini_adata, expression=cast(Any, object()))

    with pytest.raises(TypeError):
        bt.sct.pp.log1p(mini_adata, key_added=cast(Any, object()))

    with pytest.raises(ValueError):
        bt.sct.pp.log1p(mini_adata, max_memory="1XB")
