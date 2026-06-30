#!/usr/bin/env python

import numpy as np
import pytest

import bonesistools as bt


def test_pairwise_distances_returns_or_stores_matrix(
    mini_adata,
    expected_mini_pca2_distances,
):
    distances = bt.sct.tl.pairwise_distances(
        mini_adata,
        representation="X_pca",
        n_components=2,
    )

    assert distances.shape == (4, 4)
    assert np.allclose(distances, expected_mini_pca2_distances)

    result = bt.sct.tl.pairwise_distances(
        mini_adata,
        representation="X_pca",
        key_added="custom_distances",
    )
    assert result is None
    assert "custom_distances" in mini_adata.obsp
    assert mini_adata.uns["custom_distances"]["representation"] == "X_pca"


def test_barycenters_returns_cluster_means(
    mini_adata,
    expected_mini_cluster_barycenters,
):
    barycenters = bt.sct.tl.barycenters(
        mini_adata,
        obs="cluster",
        representation="X_pca",
    )

    assert sorted(barycenters) == ["A", "B"]
    assert np.allclose(barycenters["A"], expected_mini_cluster_barycenters["A"])
    assert np.allclose(barycenters["B"], expected_mini_cluster_barycenters["B"])


def test_barycenters_requires_categorical_obs(mini_adata):
    mini_adata.obs["plain"] = ["a", "a", "b", "b"]

    with pytest.raises(AttributeError, match="has no attribute 'cat'"):
        bt.sct.tl.barycenters(mini_adata, obs="plain")
