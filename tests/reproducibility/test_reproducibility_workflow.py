#!/usr/bin/env python

import os
from typing import cast

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def _run_clustering_workflow():

    adata = bt.sct.datasets.nestorowa()

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
        metric="euclidean",
        n_jobs=1,
        copy=False,
    )
    bt.sct.tl.leiden(
        adata,
        neighbors_key="neighbors",
        resolution=0.35,
        key_added="cluster",
        seed=10,
        copy=False,
    )
    bt.sct.tl.umap(
        adata,
        neighbors_key="neighbors",
        n_components=2,
        min_dist=0.5,
        spread=2.0,
        seed=10,
        n_jobs=1,
        copy=False,
    )
    bt.sct.tl.tsne(
        adata,
        representation="X_pca",
        n_pcs=20,
        n_components=2,
        max_iter=250,
        perplexity=30.0,
        metric="euclidean",
        seed=10,
        n_jobs=1,
        copy=False,
    )
    bt.sct.tl.spectral(
        adata,
        neighbors_key="neighbors",
        n_components=2,
        eigen_solver="arpack",
        seed=10,
        n_jobs=1,
        copy=False,
    )

    return adata


def _assert_same_sparse_matrix(left: csr_matrix, right: csr_matrix) -> None:

    assert left.shape == right.shape
    assert np.array_equal(left.indptr, right.indptr)
    assert np.array_equal(left.indices, right.indices)
    assert np.array_equal(left.data, right.data)


def _assert_same_value(left, right) -> None:

    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        assert isinstance(left, np.ndarray)
        assert isinstance(right, np.ndarray)
        assert np.array_equal(left, right)
    elif isinstance(left, pd.Series) or isinstance(right, pd.Series):
        assert isinstance(left, pd.Series)
        assert isinstance(right, pd.Series)
        assert left.equals(right)
    elif isinstance(left, dict) or isinstance(right, dict):
        assert isinstance(left, dict)
        assert isinstance(right, dict)
        assert left.keys() == right.keys()
        for key in left:
            _assert_same_value(left[key], right[key])
    elif isinstance(left, (list, tuple)) or isinstance(right, (list, tuple)):
        assert type(left) is type(right)
        assert len(left) == len(right)
        for left_value, right_value in zip(left, right):
            _assert_same_value(left_value, right_value)
    else:
        assert left == right


def test_clustering_workflow_is_reproducible_across_runs():

    first = _run_clustering_workflow()
    second = _run_clustering_workflow()

    assert np.array_equal(first.obsm["X_pca"], second.obsm["X_pca"])
    assert np.array_equal(first.varm["PCs"], second.varm["PCs"])
    _assert_same_value(first.uns["pca"], second.uns["pca"])

    first_distances = cast(csr_matrix, first.obsp["distances"])
    second_distances = cast(csr_matrix, second.obsp["distances"])
    first_connectivities = cast(csr_matrix, first.obsp["connectivities"])
    second_connectivities = cast(csr_matrix, second.obsp["connectivities"])

    _assert_same_sparse_matrix(first_distances, second_distances)
    _assert_same_sparse_matrix(
        first_connectivities,
        second_connectivities,
    )
    _assert_same_value(first.uns["neighbors"], second.uns["neighbors"])
    assert first_connectivities.shape == (first.n_obs, first.n_obs)

    assert first.obs["cluster"].equals(second.obs["cluster"])
    _assert_same_value(first.uns["cluster"], second.uns["cluster"])

    assert np.array_equal(first.obsm["X_umap"], second.obsm["X_umap"])
    _assert_same_value(first.uns["X_umap"], second.uns["X_umap"])

    assert np.array_equal(first.obsm["X_tsne"], second.obsm["X_tsne"])
    _assert_same_value(first.uns["X_tsne"], second.uns["X_tsne"])

    assert np.array_equal(first.obsm["X_se"], second.obsm["X_se"])
    _assert_same_value(first.uns["X_se"], second.uns["X_se"])
