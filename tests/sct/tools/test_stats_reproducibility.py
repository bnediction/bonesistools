#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd

import bonesistools as bt


def test_wilcoxon_tests_max_memory_preserves_results():
    base = ad.AnnData(
        X=np.array(
            [
                [1.0, 0.0, 2.0],
                [3.0, 0.0, 2.0],
                [2.0, 5.0, 2.0],
                [4.0, 10.0, 9.0],
                [5.0, 20.0, 8.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B", "C"])},
            index=["c1", "c2", "c3", "c4", "c5"],
        ),
        var=pd.DataFrame(index=["Gata1", "Klf1", "Tal1"]),
    )
    repeats = 10_001
    adata = ad.AnnData(
        X=np.tile(base.X, (1, repeats)),
        obs=base.obs.copy(),
        var=pd.DataFrame(
            index=[
                f"{gene}_{repeat}"
                for repeat in range(repeats)
                for gene in base.var_names
            ]
        ),
    )

    whole = bt.sct.tl.wilcoxon_tests(
        adata,
        obs="cluster",
        groups=["A"],
        correction=None,
        max_memory="100GB",
    )
    memory_limited = bt.sct.tl.wilcoxon_tests(
        adata,
        obs="cluster",
        groups=["A"],
        correction=None,
        max_memory="1MB",
    )

    np.testing.assert_allclose(
        memory_limited["statistics"].to_numpy(),
        whole["statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        memory_limited["u_statistics"].to_numpy(),
        whole["u_statistics"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        memory_limited["sum_ranks"].to_numpy(),
        whole["sum_ranks"].to_numpy(),
        equal_nan=True,
    )
