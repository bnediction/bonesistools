#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple, cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.omics.preprocessing import _hvg

GOLDEN_DIR = Path(__file__).parent
PBMC3K_PATH = GOLDEN_DIR / "pbmc3k.h5ad"
EXPECTED_DIR = GOLDEN_DIR / "expected"

QC_OBS_COLUMNS = (
    "n_features",
    "log1p_n_features",
    "total",
    "log1p_total",
    "pct_top50_features",
    "pct_top100_features",
    "pct_top200_features",
    "pct_top500_features",
)
QC_VAR_COLUMNS = (
    "n_barcodes",
    "mean",
    "log1p_mean",
    "variance",
    "log1p_variance",
    "median",
    "log1p_median",
    "mad",
    "log1p_mad",
    "pct_dropout",
    "total",
    "log1p_total",
)
HVG_N_FEATURES = 2000
PCA_N_COMPONENTS = 50
N_NEIGHBORS = 15
EMBEDDING_N_COMPONENTS = 2
TSNE_N_ITER = 300
KNNSC_CELLS_PER_CLUSTER = 50
KNNSC_N_NEIGHBORS = 15


def load_pbmc3k() -> ad.AnnData:

    return ad.read_h5ad(PBMC3K_PATH)


def run_omics_workflow() -> Dict[str, Dict[str, Any]]:

    adata = load_pbmc3k()
    outputs: Dict[str, Dict[str, Any]] = {}

    bt.omics.pp.qc(
        adata,
        percent_top=(50, 100, 200, 500),
        backend="python",
    )
    outputs["qc"] = _collect_qc(adata)

    bt.omics.pp.hvg(
        adata,
        n_features=HVG_N_FEATURES,
        method="loess",
        key_added="highly_variable",
    )
    outputs["hvg_loess"] = _collect_hvg(adata, "highly_variable")

    if "binning" in _hvg.HVG_METHODS:
        bt.omics.pp.hvg(
            adata,
            n_features=HVG_N_FEATURES,
            method=cast(Any, "binning"),
            key_added="highly_variable_binning",
        )
        outputs["hvg_binning"] = _collect_hvg(adata, "highly_variable_binning")

    bt.omics.pp.normalize(
        adata,
        target_sum=1e4,
        key_added="normalized",
    )
    bt.omics.pp.log1p(
        adata,
        expression="normalized",
        key_added="log1p",
    )

    bt.omics.tl.pca(
        adata,
        n_components=PCA_N_COMPONENTS,
        layer="log1p",
        var_subset="highly_variable",
        zero_center=True,
        svd_solver="arpack",
        seed=0,
        n_jobs=1,
    )
    outputs["pca"] = _collect_pca(adata)

    bt.omics.tl.neighbors(
        adata,
        n_neighbors=N_NEIGHBORS,
        representation="X_pca",
        n_pcs=PCA_N_COMPONENTS,
        backend="exact",
        connectivity_method="fuzzy",
        metric="euclidean",
        seed=0,
        n_jobs=1,
    )
    outputs["neighbors"] = _collect_neighbors(adata)

    bt.omics.tl.louvain(
        adata,
        resolution=1.0,
        neighbors_key="neighbors",
        key_added="knnsc_clusters",
        seed=0,
    )

    bt.omics.tl.neighbors(
        adata,
        n_neighbors=N_NEIGHBORS,
        representation="X_pca",
        n_pcs=PCA_N_COMPONENTS,
        backend="exact",
        connectivity_method="binary",
        metric="euclidean",
        key_added="spectral_neighbors",
        distances_key="spectral_distances",
        connectivities_key="spectral_connectivities",
        seed=0,
        n_jobs=1,
    )
    bt.omics.tl.spectral(
        adata,
        n_components=EMBEDDING_N_COMPONENTS,
        neighbors_key="spectral_neighbors",
        eigen_solver="arpack",
        key_added="X_spectral",
        seed=0,
        n_jobs=1,
    )
    outputs["spectral"] = _collect_embedding(adata, "X_spectral")

    bt.omics.tl.umap(
        adata,
        n_components=EMBEDDING_N_COMPONENTS,
        neighbors_key="neighbors",
        n_iter=500,
        init_pos="spectral",
        key_added="X_umap",
        seed=0,
        n_jobs=1,
    )
    outputs["umap"] = _collect_embedding(adata, "X_umap")

    knnsc_adata = _subset_knnsc_adata(adata)
    knnsc = bt.omics.tl.KNNSC()
    knnsc.fit(
        knnsc_adata,
        cluster_key="knnsc_clusters",
        representation="X_pca",
        n_components=EMBEDDING_N_COMPONENTS,
        n_neighbors=KNNSC_N_NEIGHBORS,
        n_jobs=1,
    )
    outputs["knnsc"] = _collect_knnsc(knnsc)

    bt.omics.tl.tsne(
        adata,
        representation="X_pca",
        n_pcs=PCA_N_COMPONENTS,
        n_components=EMBEDDING_N_COMPONENTS,
        n_iter=TSNE_N_ITER,
        perplexity=30.0,
        key_added="X_tsne",
        seed=0,
        n_jobs=1,
    )
    outputs["tsne"] = _collect_embedding(adata, "X_tsne")

    return outputs


def save_expected(
    outputs: Mapping[str, Mapping[str, Any]],
    output_dir: Path = EXPECTED_DIR,
) -> None:

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, values in outputs.items():
        np.savez_compressed(output_dir / f"{name}.npz", **values)


def load_expected(
    name: str,
    expected_dir: Path = EXPECTED_DIR,
) -> Dict[str, np.ndarray]:

    with np.load(expected_dir / f"{name}.npz") as expected:
        return {key: expected[key] for key in expected.files}


def available_expected_names(expected_dir: Path = EXPECTED_DIR) -> Tuple[str, ...]:

    return tuple(sorted(path.stem for path in expected_dir.glob("*.npz")))


def _collect_qc(adata: ad.AnnData) -> Dict[str, Any]:

    obs = cast(pd.DataFrame, adata.obs)
    var = cast(pd.DataFrame, adata.var)

    return {
        "obs_columns": np.asarray(QC_OBS_COLUMNS),
        "obs_values": obs.loc[:, QC_OBS_COLUMNS].to_numpy(dtype=np.float64),
        "var_columns": np.asarray(QC_VAR_COLUMNS),
        "var_values": var.loc[:, QC_VAR_COLUMNS].to_numpy(dtype=np.float64),
    }


def _collect_hvg(
    adata: ad.AnnData,
    key: str,
) -> Dict[str, Any]:

    var = cast(pd.DataFrame, adata.var)
    selected_names = var.index[var[key].to_numpy(dtype=bool)].astype(str).to_numpy()

    return {
        "mask": var[key].to_numpy(dtype=bool),
        "rank": var[f"{key}_rank"].to_numpy(dtype=np.float64),
        "score": var[f"{key}_score"].to_numpy(dtype=np.float64),
        "selected_names": np.asarray(selected_names, dtype=str),
    }


def _collect_pca(adata: ad.AnnData) -> Dict[str, Any]:

    var = cast(pd.DataFrame, adata.var)
    hvg_mask = var["highly_variable"].to_numpy(dtype=bool)

    return {
        "embedding": np.asarray(adata.obsm["X_pca"]),
        "loadings_hvg": np.asarray(adata.varm["PCs"])[hvg_mask, :],
        "variance": np.asarray(adata.uns["pca"]["variance"]),
        "variance_ratio": np.asarray(adata.uns["pca"]["variance_ratio"]),
        "hvg_indices": np.flatnonzero(hvg_mask),
    }


def _collect_neighbors(adata: ad.AnnData) -> Dict[str, Any]:

    distances = _as_csr(adata.obsp["distances"])
    connectivities = _as_csr(adata.obsp["connectivities"])

    return {
        "distances_data": distances.data,
        "distances_indices": distances.indices,
        "distances_indptr": distances.indptr,
        "distances_shape": np.asarray(distances.shape, dtype=np.int64),
        "connectivities_data": connectivities.data,
        "connectivities_indices": connectivities.indices,
        "connectivities_indptr": connectivities.indptr,
        "connectivities_shape": np.asarray(connectivities.shape, dtype=np.int64),
    }


def _collect_embedding(
    adata: ad.AnnData,
    key: str,
) -> Dict[str, Any]:

    return {
        "embedding": np.asarray(adata.obsm[key]),
    }


def _subset_knnsc_adata(adata: ad.AnnData) -> ad.AnnData:

    obs = cast(pd.DataFrame, adata.obs)
    labels = cast(pd.Series, obs["knnsc_clusters"])
    indices = []
    for cluster in labels.cat.categories:
        cluster_indices = obs.index[labels == cluster]
        indices.extend(cluster_indices[:KNNSC_CELLS_PER_CLUSTER])

    subset = adata[indices, :].copy()
    subset.obs["knnsc_clusters"] = cast(
        pd.Series,
        subset.obs["knnsc_clusters"],
    ).cat.remove_unused_categories()

    return subset


def _collect_knnsc(knnsc: bt.omics.tl.KNNSC) -> Dict[str, Any]:

    shortest_path_lengths = knnsc.shortest_path_lengths_df

    return {
        "obs_names": np.asarray(shortest_path_lengths.index.astype(str), dtype=str),
        "cluster_names": np.asarray(
            shortest_path_lengths.columns.astype(str),
            dtype=str,
        ),
        "shortest_path_lengths": shortest_path_lengths.to_numpy(dtype=np.float64),
    }


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)

    return cast(csr_matrix, matrix)


def assert_array_keys(
    result: Mapping[str, np.ndarray],
    expected: Mapping[str, np.ndarray],
) -> None:

    assert set(result) == set(expected)


def assert_equal_arrays(
    result: Mapping[str, np.ndarray],
    expected: Mapping[str, np.ndarray],
    keys: Iterable[str],
) -> None:

    for key in keys:
        np.testing.assert_array_equal(result[key], expected[key])


def assert_close_arrays(
    result: Mapping[str, np.ndarray],
    expected: Mapping[str, np.ndarray],
    keys: Iterable[str],
    *,
    rtol: float,
    atol: float,
    equal_nan: bool = False,
) -> None:

    for key in keys:
        np.testing.assert_allclose(
            result[key],
            expected[key],
            rtol=rtol,
            atol=atol,
            equal_nan=equal_nan,
        )
