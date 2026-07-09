#!/usr/bin/env python

from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt


def test_sort_reorders_observations_and_variables_stably():
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ]
        ),
        obs=pd.DataFrame({"batch": ["b", "a", "c"]}, index=["cell2", "cell1", "cell3"]),
        var=pd.DataFrame({"kind": ["z", "x", "y"]}, index=["gene3", "gene1", "gene2"]),
    )

    result = bt.sct.pp.sort(adata, on="both")

    assert result is None
    assert adata.obs_names.tolist() == ["cell1", "cell2", "cell3"]
    assert adata.var_names.tolist() == ["gene1", "gene2", "gene3"]
    np.testing.assert_array_equal(
        adata.X,
        np.array(
            [
                [5.0, 6.0, 4.0],
                [2.0, 3.0, 1.0],
                [8.0, 9.0, 7.0],
            ]
        ),
    )


def test_sort_copy_preserves_input_and_can_sort_one_axis():
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame(index=["cell2", "cell1"]),
        var=pd.DataFrame(index=["gene2", "gene1"]),
    )

    copied = bt.sct.pp.sort(adata, on="obs", copy=True)

    assert copied is not None
    assert adata.obs_names.tolist() == ["cell2", "cell1"]
    assert adata.var_names.tolist() == ["gene2", "gene1"]
    assert copied.obs_names.tolist() == ["cell1", "cell2"]
    assert copied.var_names.tolist() == ["gene2", "gene1"]


def test_sort_anndata_is_deprecated():
    adata = ad.AnnData(
        X=np.array([[1.0], [2.0]]),
        obs=pd.DataFrame(index=["cell2", "cell1"]),
        var=pd.DataFrame(index=["gene1"]),
    )

    with pytest.warns(FutureWarning, match="bt.sct.pp.sort_anndata"):
        result = bt.sct.pp.sort_anndata(adata)

    assert result is None
    assert adata.obs_names.tolist() == ["cell1", "cell2"]


def test_filter_obs_and_filter_var_subset_in_place(mini_adata):
    bt.sct.pp.filter_obs(mini_adata, "score", lambda values: values > 0.15)
    bt.sct.pp.filter_var(mini_adata, "kind", lambda values: values == "keep")

    assert mini_adata.obs_names.tolist() == ["c2", "c3", "c4"]
    assert mini_adata.var_names.tolist() == ["g1", "g3"]


def test_filter_obs_and_filter_var_validate_inputs(mini_adata):
    with pytest.raises(KeyError):
        bt.sct.pp.filter_obs(mini_adata, "missing", lambda values: values)

    with pytest.raises(TypeError):
        bt.sct.pp.filter_obs(mini_adata, cast(Any, 1), lambda values: values)

    with pytest.raises(TypeError):
        bt.sct.pp.filter_obs(mini_adata, "score", cast(Any, "not callable"))

    with pytest.raises(KeyError):
        bt.sct.pp.filter_var(mini_adata, "missing", lambda values: values)

    with pytest.raises(TypeError):
        bt.sct.pp.filter_var(mini_adata, cast(Any, 1), lambda values: values)

    with pytest.raises(TypeError):
        bt.sct.pp.filter_var(mini_adata, "kind", cast(Any, "not callable"))


def test_regress_out_removes_linear_covariate_and_preserves_copy():
    score = np.array([-1.0, 0.0, 1.0, 2.0])
    g1 = 5 + 2 * score
    g2 = 4 - score

    adata = ad.AnnData(
        X=np.zeros((4, 2)),
        obs=pd.DataFrame(
            {"score": score},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.layers["counts"] = np.column_stack([g1, g2])
    original_layer = adata.layers["counts"].copy()

    copied = bt.sct.tl.regress_out(
        adata,
        keys="score",
        layer="counts",
        intercept=True,
        copy=True,
    )

    assert np.array_equal(adata.layers["counts"], original_layer)
    assert np.allclose(
        copied.layers["counts"],
        np.array(
            [
                [5.0, 4.0],
                [5.0, 4.0],
                [5.0, 4.0],
                [5.0, 4.0],
            ]
        ),
    )

    adata.layers["sparse_counts"] = csr_matrix(np.column_stack([g1, g2]))
    sparse_copied = bt.sct.tl.regress_out(
        adata,
        keys="score",
        layer="sparse_counts",
        intercept=True,
        copy=True,
    )

    assert np.allclose(sparse_copied.layers["sparse_counts"], copied.layers["counts"])


def test_preprocessing_regress_out_is_deprecated(mini_adata):
    with pytest.warns(FutureWarning, match="bt.sct.pp.regress_out"):
        bt.sct.pp.regress_out(mini_adata, keys="score")


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

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.pp.merge(mini_adata, right, axis=cast(Any, "bad"))


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
