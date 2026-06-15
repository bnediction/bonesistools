#!/usr/bin/env python

from typing import Any, cast

import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt


def _set_two_component_graph(adata, key):
    adata.obsp[key] = csr_matrix(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )


def test_leiden_stores_clusters_and_metadata(mini_adata):
    _set_two_component_graph(mini_adata, "connectivities")

    result = bt.sct.tl.leiden(
        mini_adata,
        resolution=1.0,
        neighbors_key="neighbors",
        key_added="leiden",
        seed=10,
    )

    assert result is None
    clusters = mini_adata.obs["leiden"]
    assert isinstance(clusters.dtype, pd.CategoricalDtype)
    assert clusters.loc["c1"] == clusters.loc["c2"]
    assert clusters.loc["c3"] == clusters.loc["c4"]
    assert clusters.loc["c1"] != clusters.loc["c3"]
    assert mini_adata.uns["leiden"] == {
        "params": {
            "method": "leiden",
            "resolution": 1.0,
            "neighbors_key": "neighbors",
            "obsp": None,
            "directed": True,
            "weighted": True,
            "n_iterations": "auto",
            "seed": 10,
        }
    }


def test_leiden_obsp_copy_and_unweighted(mini_adata):
    _set_two_component_graph(mini_adata, "custom_connectivities")

    copied = bt.sct.tl.leiden(
        mini_adata,
        resolution=1.0,
        neighbors_key=None,
        obsp="custom_connectivities",
        weighted=False,
        directed=False,
        n_iterations=2,
        key_added="clusters",
        seed=10,
        copy=True,
    )

    assert "clusters" not in mini_adata.obs
    assert copied is not None
    assert copied.obs["clusters"].loc["c1"] == copied.obs["clusters"].loc["c2"]
    assert copied.obs["clusters"].loc["c3"] == copied.obs["clusters"].loc["c4"]
    assert copied.obs["clusters"].loc["c1"] != copied.obs["clusters"].loc["c3"]
    assert copied.uns["clusters"]["params"]["obsp"] == "custom_connectivities"
    assert copied.uns["clusters"]["params"]["directed"] is False
    assert copied.uns["clusters"]["params"]["weighted"] is False
    assert copied.uns["clusters"]["params"]["n_iterations"] == 2


def test_leiden_validates_graph_source_and_arguments(mini_adata):
    with pytest.raises(ValueError):
        bt.sct.tl.leiden(mini_adata, neighbors_key="neighbors", obsp="connectivities")

    with pytest.raises(KeyError):
        bt.sct.tl.leiden(mini_adata, neighbors_key="missing")

    with pytest.raises(ValueError):
        bt.sct.tl.leiden(mini_adata, resolution=0)

    with pytest.raises(ValueError):
        bt.sct.tl.leiden(mini_adata, n_iterations=0)

    with pytest.raises(ValueError):
        bt.sct.tl.leiden(mini_adata, n_iterations=-1)

    with pytest.raises(ValueError):
        bt.sct.tl.leiden(mini_adata, n_iterations=cast(Any, "bad"))

    with pytest.raises(TypeError):
        bt.sct.tl.leiden(mini_adata, directed=cast(Any, "bad"))
