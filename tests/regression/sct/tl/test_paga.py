#!/usr/bin/env python

from typing import cast

import numpy as np
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _graph as graph_tools
from tests.regression.sct.toy_data import make_nestorowa_hvg_adata


def test_paga_stores_connectivities_tree_and_metadata(mini_adata):
    result = bt.sct.tl.paga(mini_adata, groupby="cluster")

    expected_connectivities = np.array([[0.0, 0.75], [0.75, 0.0]])

    assert result is None
    assert np.allclose(
        mini_adata.uns["paga"]["connectivities"].toarray(),
        expected_connectivities,
    )
    assert mini_adata.uns["paga"]["connectivities_tree"].nnz == 1
    assert np.allclose(mini_adata.uns["paga"]["connectivities_tree"].data, [0.75])
    assert mini_adata.uns["paga"]["params"] == {
        "method": "paga",
        "groupby": "cluster",
        "neighbors_key": "neighbors",
        "obsp": None,
    }


def test_paga_uses_custom_obsp_with_copy(mini_adata):
    mini_adata.obsp["custom_distances"] = mini_adata.obsp["distances"].copy()

    copied = bt.sct.tl.paga(
        mini_adata,
        groupby="cluster",
        neighbors_key=None,
        obsp="custom_distances",
        key_added="custom_paga",
        copy=True,
    )

    assert copied is not None
    assert "custom_paga" not in mini_adata.uns
    assert np.allclose(
        copied.uns["custom_paga"]["connectivities"].toarray(),
        np.array([[0.0, 0.75], [0.75, 0.0]]),
    )
    assert copied.uns["custom_paga"]["params"] == {
        "method": "paga",
        "groupby": "cluster",
        "neighbors_key": None,
        "obsp": "custom_distances",
    }


def test_paga_validates_inputs(mini_adata):
    with pytest.raises(ValueError, match="cannot be both specified"):
        bt.sct.tl.paga(
            mini_adata,
            groupby="cluster",
            neighbors_key="neighbors",
            obsp="distances",
        )

    with pytest.raises(KeyError, match="key 'missing' not found"):
        bt.sct.tl.paga(mini_adata, groupby="cluster", neighbors_key="missing")

    missing_labels = mini_adata.copy()
    missing_labels.obs["cluster"] = missing_labels.obs["cluster"].astype(object)
    missing_labels.obs.loc["c1", "cluster"] = np.nan

    with pytest.raises(ValueError, match="non-missing group labels"):
        bt.sct.tl.paga(missing_labels, groupby="cluster")


def test_paga_matches_scanpy_v1_2_on_nestorowa():
    sc = pytest.importorskip("scanpy")

    adata = make_nestorowa_hvg_adata()
    bt.sct.tl.pca(
        adata,
        n_components=20,
        seed=10,
        copy=False,
    )
    bt.sct.tl.neighbors(
        adata,
        n_neighbors=15,
        representation="X_pca",
        n_pcs=20,
        backend="exact",
        metric="euclidean",
        n_jobs=1,
        copy=False,
    )

    scanpy_adata = adata.copy()
    bonesis_adata = adata.copy()

    sc.tl.paga(
        scanpy_adata,
        groups="clusters",
        neighbors_key="neighbors",
        model="v1.2",
        copy=False,
    )
    bt.sct.tl.paga(
        bonesis_adata,
        groupby="clusters",
        neighbors_key="neighbors",
        copy=False,
    )

    for key in ["connectivities", "connectivities_tree"]:
        scanpy_mtx = cast(csr_matrix, scanpy_adata.uns["paga"][key]).tocsr()
        bonesis_mtx = cast(csr_matrix, bonesis_adata.uns["paga"][key]).tocsr()

        assert np.array_equal(bonesis_mtx.indptr, scanpy_mtx.indptr)
        assert np.array_equal(bonesis_mtx.indices, scanpy_mtx.indices)
        assert np.array_equal(bonesis_mtx.data, scanpy_mtx.data)


def test_paga_dense_and_sparse_group_count_paths_match(mini_adata, monkeypatch):
    dense_adata = mini_adata.copy()
    sparse_adata = mini_adata.copy()

    bt.sct.tl.paga(dense_adata, groupby="cluster")
    monkeypatch.setattr(graph_tools, "_DENSE_PAGA_GROUP_PAIRS_LIMIT", 0)
    bt.sct.tl.paga(sparse_adata, groupby="cluster")

    for key in ["connectivities", "connectivities_tree"]:
        dense_mtx = cast(csr_matrix, dense_adata.uns["paga"][key]).tocsr()
        sparse_mtx = cast(csr_matrix, sparse_adata.uns["paga"][key]).tocsr()

        assert np.array_equal(dense_mtx.indptr, sparse_mtx.indptr)
        assert np.array_equal(dense_mtx.indices, sparse_mtx.indices)
        assert np.array_equal(dense_mtx.data, sparse_mtx.data)
