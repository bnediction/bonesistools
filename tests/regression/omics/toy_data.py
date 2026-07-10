#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd


def make_nestorowa_hvg_adata() -> ad.AnnData:

    rng = np.random.default_rng(2020)
    n_obs = 1200
    n_vars = 60
    labels = np.array(["HSC", "LMPP", "MPP", "CMP", "MEP", "GMP"])
    group_ids = np.arange(n_obs) % len(labels)
    X = rng.poisson(1.0, size=(n_obs, n_vars)).astype(np.float32)

    for group_id in range(len(labels)):
        group_mask = group_ids == group_id
        start = group_id * 5
        stop = start + 5
        X[group_mask, start:stop] += group_id + 1

    obs = pd.DataFrame(
        {
            "label": pd.Categorical(labels[group_ids], categories=labels),
            "clusters": pd.Categorical(group_ids.astype(str)),
        },
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {
            "ensembl": [f"ENSMUSG{i:011d}" for i in range(n_vars)],
            "symbol": [f"Gene{i}" for i in range(n_vars)],
        },
        index=[f"Gene{i}" for i in range(n_vars)],
    )

    return ad.AnnData(X=X, obs=obs, var=var)
