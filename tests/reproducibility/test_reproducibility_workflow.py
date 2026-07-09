#!/usr/bin/env python

import os
import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
from umap.umap_ import fuzzy_simplicial_set

import bonesistools as bt
from tests.regression.sct.toy_data import make_nestorowa_hvg_adata

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def _load_nestorowa_hvg() -> ad.AnnData:

    return make_nestorowa_hvg_adata()


def _run_clustering_workflow():

    adata = _load_nestorowa_hvg()

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
        n_iter=250,
        perplexity=30.0,
        metric="euclidean",
        seed=10,
        n_jobs=1,
        copy=False,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(
                "Graph is not fully connected, spectral embedding may not "
                "work as expected."
            ),
            category=UserWarning,
            module="sklearn.manifold._spectral_embedding",
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


def test_neighbors_connectivities_match_umap_reference_on_nestorowa():

    adata = _load_nestorowa_hvg()
    bt.sct.tl.pca(
        adata,
        n_components=20,
        seed=10,
        copy=False,
    )

    representation_mtx = cast(np.ndarray, adata.obsm["X_pca"])[:, :20]
    neighbors_model = NearestNeighbors(
        n_neighbors=15,
        metric="euclidean",
        n_jobs=1,
    )
    neighbors_model.fit(representation_mtx)
    knn_distances, knn_indices = neighbors_model.kneighbors(representation_mtx)
    reference_graph = cast(
        Any,
        fuzzy_simplicial_set(
            X=representation_mtx,
            n_neighbors=15,
            random_state=np.random.RandomState(0),
            metric="euclidean",
            knn_indices=knn_indices,
            knn_dists=knn_distances,
            set_op_mix_ratio=1.0,
            local_connectivity=1.0,
            apply_set_operations=True,
            verbose=False,
        )[0],
    )
    reference_connectivities = cast(
        csr_matrix,
        reference_graph.tocsr(),
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
    connectivities = cast(csr_matrix, adata.obsp["connectivities"])

    assert np.array_equal(connectivities.indptr, reference_connectivities.indptr)
    assert np.array_equal(connectivities.indices, reference_connectivities.indices)
    np.testing.assert_allclose(
        connectivities.data,
        reference_connectivities.data,
        rtol=3e-5,
        atol=2e-6,
    )


@pytest.mark.parametrize("backend", ["exact", "pynndescent"])
def test_neighbors_are_reproducible_on_nestorowa(backend):

    if backend == "pynndescent":
        pytest.importorskip("pynndescent")

    first = _load_nestorowa_hvg()
    second = _load_nestorowa_hvg()

    for adata in [first, second]:
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
            backend=backend,
            metric="euclidean",
            seed=10,
            n_jobs=1,
            copy=False,
        )

    _assert_same_sparse_matrix(
        cast(csr_matrix, first.obsp["distances"]),
        cast(csr_matrix, second.obsp["distances"]),
    )
    _assert_same_sparse_matrix(
        cast(csr_matrix, first.obsp["connectivities"]),
        cast(csr_matrix, second.obsp["connectivities"]),
    )
    _assert_same_value(first.uns["neighbors"], second.uns["neighbors"])


def test_qc_numba_backend_is_reproducible_on_nestorowa():

    pytest.importorskip("numba")

    first = _load_nestorowa_hvg()
    second = _load_nestorowa_hvg()
    first.X = csr_matrix(first.X)
    second.X = csr_matrix(second.X)

    for adata in [first, second]:
        bt.sct.pp.qc(
            adata,
            percent_top=[10, 20, 50],
            backend="numba",
            copy=False,
        )

    assert cast(pd.DataFrame, first.obs).equals(cast(pd.DataFrame, second.obs))
    assert cast(pd.DataFrame, first.var).equals(cast(pd.DataFrame, second.var))
