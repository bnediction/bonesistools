#!/usr/bin/env python

import numpy as np
import pandas as pd
import pytest

import bonesistools as bt


def test_filter_obs_and_filter_var_subset_in_place(mini_adata):
    bt.sct.pp.filter_obs(mini_adata, "score", lambda values: values > 0.15)
    bt.sct.pp.filter_var(mini_adata, "kind", lambda values: values == "keep")

    assert mini_adata.obs_names.tolist() == ["c2", "c3", "c4"]
    assert mini_adata.var_names.tolist() == ["g1", "g3"]


def test_filter_obs_and_filter_var_validate_inputs(mini_adata):
    with pytest.raises(KeyError, match="key 'missing' not found in adata.obs"):
        bt.sct.pp.filter_obs(mini_adata, "missing", lambda values: values)

    with pytest.raises(TypeError, match="unsupported argument type for 'function'"):
        bt.sct.pp.filter_var(mini_adata, "kind", "not callable")


def test_regress_out_updates_layer_or_returns_copy(mini_adata):
    original_layer = mini_adata.layers["counts"].copy()

    copied = bt.sct.pp.regress_out(
        mini_adata,
        keys="score",
        layer="counts",
        copy=True,
    )

    assert np.array_equal(mini_adata.layers["counts"], original_layer)
    assert copied.layers["counts"].shape == original_layer.shape
    assert not np.allclose(copied.layers["counts"], original_layer)


def test_merge_obs_and_var(mini_adata):
    right = mini_adata.copy()
    right.obs = pd.DataFrame(
        {"annotation": ["x", "y", "z", "w"]},
        index=mini_adata.obs_names,
    )
    right.var = pd.DataFrame(
        {"symbol": ["G1", "G2", "G3"]},
        index=mini_adata.var_names,
    )

    merged_obs = bt.sct.pp.merge(mini_adata, right, axis="obs", copy=True)
    merged_var = bt.sct.pp.merge(mini_adata, right, axis="var", copy=True)

    assert merged_obs.obs["annotation"].tolist() == ["x", "y", "z", "w"]
    assert merged_var.var["symbol"].tolist() == ["G1", "G2", "G3"]


def test_transfer_layer_aligns_obs_and_var(mini_adata):
    right = mini_adata[[1, 3], [2, 0]].copy()
    right.layers["counts"] = np.array([[10.0, 11.0], [20.0, 21.0]])

    bt.sct.pp.transfer_layer(mini_adata, right, layers="counts")

    transferred = pd.DataFrame(
        mini_adata.layers["counts"],
        index=mini_adata.obs_names,
        columns=mini_adata.var_names,
    )
    assert transferred.loc["c2", "g3"] == 10.0
    assert transferred.loc["c4", "g1"] == 21.0
    assert pd.isna(transferred.loc["c1", "g1"])
    assert pd.isna(transferred.loc["c2", "g2"])


def test_transfer_layer_rejects_invalid_layers(mini_adata):
    with pytest.raises(TypeError, match="unsupported argument type for 'layers'"):
        bt.sct.pp.transfer_layer(mini_adata, mini_adata.copy(), layers=("counts",))
