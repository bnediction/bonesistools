#!/usr/bin/env python

import pandas as pd
from scipy import io, sparse

import bonesistools as bt


def test_to_csv_writes_dense_matrix(mini_adata, tmp_path):
    filename = tmp_path / "matrix"

    bt.sct.tl.to_csv(mini_adata, filename)

    written = pd.read_csv(tmp_path / "matrix.csv", index_col=0)
    assert written.shape == mini_adata.X.shape
    assert written.iloc[0, 0] == mini_adata.X[0, 0]


def test_to_mtx_writes_layer_matrix(mini_adata, tmp_path):
    filename = tmp_path / "counts"

    bt.sct.tl.to_mtx(mini_adata, filename, layer="counts")

    written = io.mmread(tmp_path / "counts.mtx")
    assert written.shape == mini_adata.layers["counts"].shape


def test_to_npz_writes_sparse_matrix(mini_adata, tmp_path):
    mini_adata.layers["sparse_counts"] = sparse.csr_matrix(mini_adata.layers["counts"])
    filename = tmp_path / "sparse_counts"

    bt.sct.tl.to_npz(mini_adata, filename, layer="sparse_counts")

    written = sparse.load_npz(tmp_path / "sparse_counts.npz")
    assert written.shape == mini_adata.layers["sparse_counts"].shape
    assert written.nnz == mini_adata.layers["sparse_counts"].nnz


def test_to_csv_or_mtx_selects_writer_from_density(mini_adata, tmp_path):
    bt.sct.tl.to_csv_or_mtx(mini_adata, tmp_path / "dense")

    mini_adata.layers["sparse_counts"] = sparse.csr_matrix(mini_adata.layers["counts"])
    bt.sct.tl.to_csv_or_mtx(mini_adata, tmp_path / "sparse", layer="sparse_counts")

    assert (tmp_path / "dense.csv").exists()
    assert (tmp_path / "sparse.mtx").exists()


def test_to_csv_or_npz_selects_writer_from_density(mini_adata, tmp_path):
    bt.sct.tl.to_csv_or_npz(mini_adata, tmp_path / "dense")

    mini_adata.layers["sparse_counts"] = sparse.csr_matrix(mini_adata.layers["counts"])
    bt.sct.tl.to_csv_or_npz(mini_adata, tmp_path / "sparse", layer="sparse_counts")

    assert (tmp_path / "dense.csv").exists()
    assert (tmp_path / "sparse.npz").exists()


def test_writers_cover_layer_and_sparse_x_branches(mini_adata, tmp_path):
    mini_adata.layers["dense_counts"] = mini_adata.layers["counts"].copy()

    bt.sct.tl.to_csv(mini_adata, tmp_path / "dense_layer", layer="dense_counts")
    bt.sct.tl.to_csv_or_mtx(
        mini_adata,
        tmp_path / "auto_dense_layer_mtx",
        layer="dense_counts",
    )
    bt.sct.tl.to_csv_or_npz(
        mini_adata,
        tmp_path / "auto_dense_layer_npz",
        layer="dense_counts",
    )

    mini_adata.X = sparse.csr_matrix(mini_adata.X)

    bt.sct.tl.to_mtx(mini_adata, tmp_path / "sparse_x_mtx")
    bt.sct.tl.to_npz(mini_adata, tmp_path / "sparse_x_npz")
    bt.sct.tl.to_csv_or_mtx(mini_adata, tmp_path / "auto_sparse_x_mtx")
    bt.sct.tl.to_csv_or_npz(mini_adata, tmp_path / "auto_sparse_x_npz")

    assert (tmp_path / "dense_layer.csv").exists()
    assert (tmp_path / "auto_dense_layer_mtx.csv").exists()
    assert (tmp_path / "auto_dense_layer_npz.csv").exists()
    assert (tmp_path / "sparse_x_mtx.mtx").exists()
    assert (tmp_path / "sparse_x_npz.npz").exists()
    assert (tmp_path / "auto_sparse_x_mtx.mtx").exists()
    assert (tmp_path / "auto_sparse_x_npz.npz").exists()
