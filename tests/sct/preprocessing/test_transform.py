#!/usr/bin/env python

from typing import Any, Optional, cast

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


def _scbolt_reference_normalize(
    adata: ad.AnnData,
    layer: str,
    target_sum: float = 1e4,
) -> None:

    matrix = adata.layers[layer]

    if sparse.issparse(matrix):
        matrix = cast(csr_matrix, cast(Any, matrix).tocsr(copy=True))
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
    reference = ad.AnnData(X=np.zeros((4, 3), dtype=np.float32))
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
