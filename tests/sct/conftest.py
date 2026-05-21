#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


@pytest.fixture
def mini_adata():
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 0.0, 3.0],
                [2.0, 1.0, 0.0],
                [0.0, 3.0, 1.0],
                [4.0, 0.0, 2.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(["A", "A", "B", "B"]),
                "batch": ["b1", "b2", "b1", "b2"],
                "score": [0.1, 0.7, 0.2, 0.8],
            },
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(
            {"kind": ["keep", "drop", "keep"]},
            index=["g1", "g2", "g3"],
        ),
    )
    adata.layers["counts"] = adata.X.copy() + 1.0
    adata.obsm["X_pca"] = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.2, 0.1, 1.1],
            [2.0, 2.0, 0.0],
            [2.2, 2.1, 0.1],
        ]
    )
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {
            "n_neighbors": 3,
            "n_pcs": 2,
            "use_rep": "X_pca",
        },
    }
    adata.obsp["distances"] = csr_matrix(
        [
            [0.0, 1.0, 2.0, 0.0],
            [1.0, 0.0, 0.0, 2.0],
            [2.0, 0.0, 0.0, 1.0],
            [0.0, 2.0, 1.0, 0.0],
        ]
    )
    adata.obsp["connectivities"] = csr_matrix(
        [
            [1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 1.0],
            [0.0, 1.0, 1.0, 1.0],
        ]
    )
    return adata
