#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import bonesistools as bt


def _dense_list(value: object) -> list:
    matrix = cast(Any, value)
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).tolist()


def _duplicated_adata(**kwargs: Any) -> ad.AnnData:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        return ad.AnnData(**kwargs)


def test_merge_duplicate_vars_returns_independent_copy_without_duplicates():
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]]),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame({"tag": ["first", "second"]}, index=["g1", "g2"]),
        uns={"nested": {"value": 1}},
    )
    adata.obsm["coords"] = np.array([[1.0, 2.0], [3.0, 4.0]])
    adata.obsp["distances"] = np.array([[0.0, 1.0], [1.0, 0.0]])

    merged = bt.omics.pp.merge_duplicate_vars(adata, copy=True)

    assert merged is not adata
    assert merged.var_names.tolist() == ["g1", "g2"]
    assert _dense_list(merged.X) == [[1.0, 2.0], [3.0, 4.0]]

    merged.uns["nested"]["value"] = 2
    merged_coords = cast(np.ndarray, merged.obsm["coords"])
    merged_distances = cast(np.ndarray, merged.obsp["distances"])
    coords = cast(np.ndarray, adata.obsm["coords"])
    distances = cast(np.ndarray, adata.obsp["distances"])

    merged_coords[0, 0] = 10.0
    merged_distances[0, 1] = 5.0

    assert adata.uns["nested"]["value"] == 1
    np.testing.assert_array_equal(coords, np.array([[1.0, 2.0], [3.0, 4.0]]))
    np.testing.assert_array_equal(distances, np.array([[0.0, 1.0], [1.0, 0.0]]))


def test_merge_duplicate_vars_preserves_axis_mappings_and_metadata():
    adata = _duplicated_adata(
        X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        obs=pd.DataFrame({"condition": ["ctrl", "treated"]}, index=["c1", "c2"]),
        var=pd.DataFrame(
            {"symbol": ["Aars", "Aars", "Myc"], "tag": ["same", "same", "unique"]},
            index=pd.Index(["Aars", "Aars", "Myc"], name="gene"),
        ),
        uns={"pipeline": {"step": "normalization"}},
    )
    adata.obsm["embedding"] = np.array([[0.1, 0.2], [0.3, 0.4]])
    adata.obsp["connectivities"] = np.array([[1.0, 0.5], [0.5, 1.0]])

    merged = bt.omics.pp.merge_duplicate_vars(adata, keep="consensus")

    assert merged.var_names.tolist() == ["Aars", "Myc"]
    assert merged.var_names.name == "gene"
    assert _dense_list(merged.X) == [[3.0, 3.0], [9.0, 6.0]]
    assert cast(pd.DataFrame, merged.obs).equals(cast(pd.DataFrame, adata.obs))
    assert _dense_list(merged.obsm["embedding"]) == [[0.1, 0.2], [0.3, 0.4]]
    assert _dense_list(merged.obsp["connectivities"]) == [[1.0, 0.5], [0.5, 1.0]]
    assert merged.uns == {"pipeline": {"step": "normalization"}}


def test_merge_duplicate_vars_consensus_treats_matching_nan_as_consensus():
    adata = _duplicated_adata(
        X=np.array([[1.0, 2.0, 3.0]]),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(
            {
                "nullable": [np.nan, np.nan, "defined"],
                "conflict": ["left", "right", "stable"],
            },
            index=["g1", "g1", "g2"],
        ),
    )

    merged = bt.omics.pp.merge_duplicate_vars(adata, keep="consensus")
    var = cast(pd.DataFrame, merged.var)

    assert pd.isna(var.loc["g1", "nullable"])
    assert pd.isna(var.loc["g1", "conflict"])
    assert var.loc["g2", "nullable"] == "defined"
    assert var.loc["g2", "conflict"] == "stable"


def test_merge_duplicate_vars_dataframe_varm_mean_and_nan_strategies():
    adata = _duplicated_adata(
        X=np.ones((1, 3)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1", "g2"]),
    )
    adata.varm["numeric"] = pd.DataFrame(
        {"score": [1.0, 3.0, 5.0], "weight": [2.0, 4.0, 6.0]},
        index=["g1", "g1", "g2"],
    )

    merged_mean = bt.omics.pp.merge_duplicate_vars(adata, varm="mean")
    numeric = merged_mean.varm["numeric"]

    assert isinstance(numeric, pd.DataFrame)
    assert numeric.index.tolist() == ["g1", "g2"]
    assert numeric.to_dict("list") == {"score": [2.0, 5.0], "weight": [3.0, 6.0]}

    adata_with_labels = _duplicated_adata(
        X=np.ones((1, 3)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1", "g2"]),
    )
    adata_with_labels.varm["labels"] = pd.DataFrame(
        {"label": ["left", "right", "single"]},
        index=["g1", "g1", "g2"],
    )

    merged_nan = bt.omics.pp.merge_duplicate_vars(adata_with_labels, varm="nan")
    labels = merged_nan.varm["labels"]

    assert isinstance(labels, pd.DataFrame)
    assert pd.isna(labels.loc["g1", "label"])
    assert labels.loc["g2", "label"] == "single"


def test_merge_duplicate_vars_rejects_non_numeric_varm_mean():
    adata_with_array = _duplicated_adata(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1"]),
    )
    adata_with_array.varm["labels"] = np.array([["left"], ["right"]], dtype=object)

    with pytest.raises(TypeError, match="varm='mean' requires numeric values"):
        bt.omics.pp.merge_duplicate_vars(adata_with_array, varm="mean")

    adata_with_dataframe = _duplicated_adata(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1"]),
    )
    adata_with_dataframe.varm["labels"] = pd.DataFrame(
        {"label": ["left", "right"]},
        index=["g1", "g1"],
    )

    with pytest.raises(TypeError, match="varm='mean' requires numeric values"):
        bt.omics.pp.merge_duplicate_vars(adata_with_dataframe, varm="mean")


def test_merge_duplicate_vars_varm_nan_handles_bool_and_object_arrays():
    adata = _duplicated_adata(
        X=np.ones((1, 3)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1", "g2"]),
    )
    adata.varm["flag"] = np.array([[True], [False], [True]])
    adata.varm["label"] = np.array([["left"], ["right"], ["single"]], dtype=object)

    merged = bt.omics.pp.merge_duplicate_vars(adata, varm="nan")
    flag = cast(np.ndarray, merged.varm["flag"])
    label = cast(np.ndarray, merged.varm["label"])

    assert flag.dtype.kind == "f"
    np.testing.assert_array_equal(np.isnan(flag), np.array([[True], [False]]))
    np.testing.assert_allclose(flag[1, 0], 1.0)
    missing_label = cast(float, label[0, 0])
    single_label = cast(str, label[1, 0])

    assert np.isnan(missing_label)
    assert single_label == "single"


def test_merge_duplicate_vars_validates_strategies_and_deprecated_alias():
    adata = _duplicated_adata(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["g1", "g1"]),
    )

    with pytest.raises(ValueError, match="invalid argument value for 'keep'"):
        bt.omics.pp.merge_duplicate_vars(adata, keep=cast(Any, "last"))

    with pytest.raises(ValueError, match="invalid argument value for 'varm'"):
        bt.omics.pp.merge_duplicate_vars(adata, varm=cast(Any, "drop"))

    with pytest.raises(ValueError, match="invalid argument value for 'varp'"):
        bt.omics.pp.merge_duplicate_vars(adata, varp=cast(Any, "keep"))

    with pytest.warns(FutureWarning):
        merged = getattr(bt.omics.pp, "var_names_merge_duplicates")(
            adata,
            var_names_column="symbol",
        )

    assert merged.var_names.tolist() == ["g1"]
    assert _dense_list(merged.X) == [[2.0]]
