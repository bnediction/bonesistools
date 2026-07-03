#!/usr/bin/env python

from pathlib import Path
from typing import Any, cast

import anndata as ad
from scipy import sparse
from scipy.sparse import csr_matrix

_GOLDEN_DIR = Path(__file__).parent


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)

    return cast(csr_matrix, matrix)


def test_pbmc3k_golden_dataset_is_loadable_and_raw():
    adata = ad.read_h5ad(_GOLDEN_DIR / "pbmc3k.h5ad")
    counts = _as_csr(adata.X)

    assert adata.shape == (2700, 32738)
    assert adata.obs_names[:3].tolist() == [
        "AAACATACAACCAC-1",
        "AAACATTGAGCTAC-1",
        "AAACATTGATCAGC-1",
    ]
    assert adata.var_names[:3].tolist() == ["MIR1302-10", "FAM138A", "OR4F5"]
    assert counts.nnz == 2286884
    assert str(counts.dtype) == "float64"
    assert list(adata.obs.columns) == []
    assert list(adata.var.columns) == ["Accession", "symbol"]
    assert adata.uns["pbmc3k"]["license"] == "CC BY 4.0"
