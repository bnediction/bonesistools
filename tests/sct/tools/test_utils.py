#!/usr/bin/env python

from typing import Any, cast

import numpy as np
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _utils


def test_get_expression_selects_x_layer_and_raw(mini_adata):
    mini_adata.raw = mini_adata
    expression_mtx = cast(np.ndarray, mini_adata.X)
    counts_mtx = cast(np.ndarray, mini_adata.layers["counts"])
    raw_mtx = cast(np.ndarray, mini_adata.raw.X)

    expression_mtx_copy = cast(np.ndarray, bt.sct.tl.get_expression(mini_adata))

    assert np.array_equal(expression_mtx_copy, expression_mtx)
    assert expression_mtx_copy is not expression_mtx
    assert np.array_equal(
        cast(np.ndarray, bt.sct.tl.get_expression(mini_adata, copy=False)),
        expression_mtx,
    )
    assert np.array_equal(
        cast(
            np.ndarray,
            bt.sct.tl.get_expression(
                mini_adata,
                layer="counts",
                copy=False,
            ),
        ),
        counts_mtx,
    )
    assert np.array_equal(
        cast(
            np.ndarray,
            bt.sct.tl.get_expression(
                mini_adata,
                use_raw=True,
                copy=False,
            ),
        ),
        raw_mtx,
    )


def test_get_expression_rejects_raw_and_layer(mini_adata):
    with pytest.raises(ValueError, match="invalid argument combination"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            layer="counts",
        )


def test_get_expression_selects_variable_subset(mini_adata):
    raw_adata = mini_adata.copy()
    raw_adata.var_names = ["raw_g1", "raw_g2", "raw_g3"]
    mini_adata.raw = raw_adata
    mini_adata.var["selected"] = [True, False, True]
    expression_mtx = cast(np.ndarray, mini_adata.X)
    raw_mtx = cast(np.ndarray, mini_adata.raw.X)

    subset_by_column = cast(
        np.ndarray,
        bt.sct.tl.get_expression(
            mini_adata,
            var_subset="selected",
            copy=False,
        ),
    )
    subset_by_names = cast(
        np.ndarray,
        bt.sct.tl.get_expression(
            mini_adata,
            var_subset=["g2", "g1"],
            copy=False,
        ),
    )
    raw_subset = cast(
        np.ndarray,
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset=["raw_g3", "raw_g1"],
            copy=False,
        ),
    )

    assert np.array_equal(subset_by_column, expression_mtx[:, [0, 2]])
    assert np.array_equal(subset_by_names, expression_mtx[:, [0, 1]])
    assert np.array_equal(raw_subset, raw_mtx[:, [0, 2]])


def test_get_expression_validates_raw_variable_subset(mini_adata):
    raw_adata = mini_adata.copy()
    raw_adata.var_names = ["raw_g1", "raw_g2", "raw_g3"]
    raw_adata.var["selected"] = [True, False, True]
    raw_adata.var["not_bool"] = ["yes", "no", "yes"]
    mini_adata.raw = raw_adata

    with pytest.raises(KeyError, match="column 'missing' not found"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset="missing",
        )

    with pytest.raises(TypeError, match="unsupported column dtype"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset="not_bool",
        )

    with pytest.raises(ValueError, match="expected at least one variable name"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset=[],
        )

    with pytest.raises(TypeError, match="unsupported element type"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset=cast(Any, ["raw_g1", 1]),
        )

    with pytest.raises(KeyError, match="variable\\(s\\) not found"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset=["raw_g1", "missing"],
        )

    with pytest.raises(TypeError, match="unsupported argument type"):
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            var_subset=cast(Any, object()),
        )


def test_get_expression_with_gene_names_resolves_sources_and_subsets(mini_adata):
    raw_adata = mini_adata.copy()
    raw_adata.var_names = ["raw_g1", "raw_g2", "raw_g3"]
    raw_adata.var["selected"] = [True, False, True]
    mini_adata.raw = raw_adata
    mini_adata.layers["log"] = cast(np.ndarray, mini_adata.X).copy()
    mini_adata.var["selected"] = [False, True, True]

    expression_mtx, gene_names = _utils._get_expression_with_gene_names(
        mini_adata,
        expression="log",
        var_subset=["g3", "g1"],
    )
    raw_mtx, raw_names = _utils._get_expression_with_gene_names(
        mini_adata,
        expression="raw.X",
        var_subset="selected",
    )

    assert expression_mtx.shape == (mini_adata.n_obs, 2)
    assert gene_names.tolist() == ["g1", "g3"]
    assert raw_mtx.shape == (mini_adata.n_obs, 2)
    assert raw_names.tolist() == ["raw_g1", "raw_g3"]


def test_get_expression_with_gene_names_validates_subset(mini_adata):
    mini_adata.var["empty"] = [False, False, False]
    mini_adata.var["not_bool"] = ["yes", "no", "yes"]

    with pytest.raises(ValueError, match="adata.raw is required"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression="raw.X",
            var_subset=None,
        )

    with pytest.raises(KeyError, match="column 'missing' not found"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset="missing",
        )

    with pytest.raises(TypeError, match="unsupported column dtype"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset="not_bool",
        )

    with pytest.raises(ValueError, match="selects no variables"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset="empty",
        )

    with pytest.raises(TypeError, match="unsupported argument type"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset=cast(Any, object()),
        )

    with pytest.raises(ValueError, match="expected at least one variable name"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset=[],
        )

    with pytest.raises(TypeError, match="unsupported element type"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset=cast(Any, ["g1", 1]),
        )

    with pytest.raises(KeyError, match="variable\\(s\\) not found"):
        _utils._get_expression_with_gene_names(
            mini_adata,
            expression=None,
            var_subset=["g1", "missing"],
        )


def test_as_dense_matrix_chunk_supports_sparse_and_rejects_non_matrix():
    sparse_chunk = _utils._as_dense_matrix_chunk(
        csr_matrix([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]]),
        start=1,
        end=3,
    )

    np.testing.assert_array_equal(sparse_chunk, np.array([[0.0, 2.0], [3.0, 0.0]]))

    with pytest.raises(ValueError, match="expected a two-dimensional matrix"):
        _utils._as_dense_matrix_chunk(np.array([1.0, 2.0, 3.0]), start=0, end=2)


def test_get_representation_truncates_dimensions(mini_adata):
    pca_mtx = cast(np.ndarray, mini_adata.obsm["X_pca"])
    full_representation_mtx = bt.sct.tl.get_representation(mini_adata, use_rep=None)
    representation_mtx = cast(
        np.ndarray,
        bt.sct.tl.get_representation(
            mini_adata,
            use_rep="X_pca",
            n_components=2,
        ),
    )

    assert full_representation_mtx is pca_mtx
    assert representation_mtx.shape == (4, 2)
    assert np.array_equal(representation_mtx, pca_mtx[:, :2])


def test_get_representation_reports_missing_key(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found in scdata.obsm"):
        bt.sct.tl.get_representation(mini_adata, use_rep="missing")

    pca_missing = mini_adata.copy()
    del pca_missing.obsm["X_pca"]

    with pytest.raises(KeyError, match="bonesistools.sct.tl.pca"):
        bt.sct.tl.get_representation(pca_missing, use_rep=None)


def test_get_pairwise_selects_obsp_and_varp(mini_adata):
    mini_adata.varp["correlations"] = csr_matrix(
        [
            [1.0, 0.5, 0.0],
            [0.5, 1.0, 0.2],
            [0.0, 0.2, 1.0],
        ]
    )

    obs_pairwise_mtx = bt.sct.tl.get_pairwise(mini_adata, "connectivities")
    var_pairwise_mtx = bt.sct.tl.get_pairwise(
        mini_adata,
        "correlations",
        axis="var",
    )

    assert obs_pairwise_mtx is mini_adata.obsp["connectivities"]
    assert var_pairwise_mtx is mini_adata.varp["correlations"]


def test_get_pairwise_rejects_invalid_axis(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.get_pairwise(
            mini_adata,
            "connectivities",
            axis=cast(Any, "bad"),
        )


def test_get_distances_and_connectivities_from_default_neighbors(mini_adata):
    assert _utils._get_distances(mini_adata) is mini_adata.obsp["distances"]
    assert _utils._get_connectivities(mini_adata) is mini_adata.obsp["connectivities"]


def test_get_distances_and_connectivities_from_explicit_keys(mini_adata):
    mini_adata.obsp["custom_distances"] = mini_adata.obsp["distances"].copy()
    mini_adata.obsp["custom_connectivities"] = mini_adata.obsp["connectivities"].copy()
    mini_adata.uns["custom_neighbors"] = {
        "distances_key": "custom_distances",
        "connectivities_key": "custom_connectivities",
    }

    assert (
        _utils._get_distances(
            mini_adata,
            obsp="custom_distances",
        )
        is mini_adata.obsp["custom_distances"]
    )
    assert (
        _utils._get_connectivities(
            mini_adata,
            obsp="custom_connectivities",
        )
        is mini_adata.obsp["custom_connectivities"]
    )
    assert (
        _utils._get_distances(
            mini_adata,
            neighbors_key="custom_neighbors",
        )
        is mini_adata.obsp["custom_distances"]
    )
    assert (
        _utils._get_connectivities(
            mini_adata,
            neighbors_key="custom_neighbors",
        )
        is mini_adata.obsp["custom_connectivities"]
    )


@pytest.mark.parametrize(
    "getter,obsp_key,neighbors_key,error_message",
    [
        (
            _utils._get_distances,
            "distances",
            "neighbors",
            "'obsp' and 'neighbors_key' cannot be both specified",
        ),
        (
            _utils._get_connectivities,
            "connectivities",
            "neighbors",
            "'obsp' and 'neighbors_key' cannot be both specified",
        ),
    ],
)
def test_get_distances_and_connectivities_reject_ambiguous_keys(
    mini_adata,
    getter,
    obsp_key,
    neighbors_key,
    error_message,
):
    with pytest.raises(ValueError, match=error_message):
        getter(mini_adata, obsp=obsp_key, neighbors_key=neighbors_key)


@pytest.mark.parametrize(
    "getter,error_message",
    [
        (_utils._get_distances, "distances not found"),
        (_utils._get_connectivities, "connectivities not found"),
    ],
)
def test_get_distances_and_connectivities_report_missing_defaults(
    mini_adata,
    getter,
    error_message,
):
    mini_adata.uns.clear()

    with pytest.raises(KeyError, match=error_message):
        getter(mini_adata)
