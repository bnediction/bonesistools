#!/usr/bin/env python

import numpy as np
from scipy.sparse import csr_matrix

import bonesistools as bt


def test_anndata_to_dataframe_uses_layer_obs_and_log_base(mini_adata):
    mini_adata.uns["log1p"] = {"base": 2}
    mini_adata.layers["log2_counts"] = np.log1p(mini_adata.layers["counts"]) / np.log(2)

    df = bt.sct.tl.anndata_to_dataframe(
        mini_adata,
        obs=["cluster", "batch"],
        layer="log2_counts",
        is_log=True,
    )

    assert df.columns.tolist() == ["g1", "g2", "g3", "cluster", "batch"]
    assert np.allclose(df.loc[:, ["g1", "g2", "g3"]], mini_adata.layers["counts"])
    assert df["batch"].tolist() == ["b1", "b2", "b1", "b2"]


def test_anndata_to_dataframe_handles_sparse_x_sparse_layer_and_obs_string(mini_adata):
    sparse_adata = mini_adata.copy()
    sparse_adata.X = csr_matrix(sparse_adata.X)

    x_df = bt.sct.tl.anndata_to_dataframe(sparse_adata)

    assert np.allclose(x_df, mini_adata.X)

    sparse_adata.layers["sparse_log_counts"] = csr_matrix(
        np.log1p(mini_adata.layers["counts"])
    )
    layer_df = bt.sct.tl.anndata_to_dataframe(
        sparse_adata,
        obs="cluster",
        layer="sparse_log_counts",
        is_log=True,
    )

    assert layer_df.columns.tolist() == ["g1", "g2", "g3", "cluster"]
    assert np.allclose(
        layer_df.loc[:, ["g1", "g2", "g3"]],
        mini_adata.layers["counts"],
    )
    assert layer_df["cluster"].tolist() == ["A", "A", "B", "B"]
