#!/usr/bin/env python

import numpy as np
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _utils


def test_get_expression_selects_x_layer_and_raw(mini_adata):
    mini_adata.raw = mini_adata

    expression_mtx_copy = bt.sct.tl.get_expression(mini_adata)

    assert np.array_equal(expression_mtx_copy, mini_adata.X)
    assert expression_mtx_copy is not mini_adata.X
    assert np.array_equal(
        bt.sct.tl.get_expression(mini_adata, copy=False),
        mini_adata.X,
    )
    assert np.array_equal(
        bt.sct.tl.get_expression(
            mini_adata,
            layer="counts",
            copy=False,
        ),
        mini_adata.layers["counts"],
    )
    assert np.array_equal(
        bt.sct.tl.get_expression(
            mini_adata,
            use_raw=True,
            copy=False,
        ),
        mini_adata.raw.X,
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

    subset_by_column = bt.sct.tl.get_expression(
        mini_adata,
        var_subset="selected",
        copy=False,
    )
    subset_by_names = bt.sct.tl.get_expression(
        mini_adata,
        var_subset=["g2", "g1"],
        copy=False,
    )
    raw_subset = bt.sct.tl.get_expression(
        mini_adata,
        use_raw=True,
        var_subset=["raw_g3", "raw_g1"],
        copy=False,
    )

    assert np.array_equal(subset_by_column, mini_adata.X[:, [0, 2]])
    assert np.array_equal(subset_by_names, mini_adata.X[:, [0, 1]])
    assert np.array_equal(raw_subset, mini_adata.raw.X[:, [0, 2]])


def test_get_representation_truncates_dimensions(mini_adata):
    full_representation_mtx = bt.sct.tl.get_representation(mini_adata, use_rep=None)
    representation_mtx = bt.sct.tl.get_representation(
        mini_adata,
        use_rep="X_pca",
        n_components=2,
    )

    assert full_representation_mtx is mini_adata.obsm["X_pca"]
    assert representation_mtx.shape == (4, 2)
    assert np.array_equal(representation_mtx, mini_adata.obsm["X_pca"][:, :2])


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
