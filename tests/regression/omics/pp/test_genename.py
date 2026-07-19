#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.omics.preprocessing import _genename


def test_convert_gene_identifiers_copy_and_axis_validation(
    monkeypatch,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_genename, "create_identifiers", fake_gene_synonyms_cls)

    adata = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["NF-kappaB", "unknown"]),
        var=pd.DataFrame(index=["Tp53", "Myc"]),
    )

    converted_var = bt.omics.pp.convert_gene_identifiers(adata, axis="var", copy=True)
    converted_obs = bt.omics.pp.convert_gene_identifiers(
        adata,
        axis="obs",
        copy=True,
    )

    assert converted_var.var_names.tolist() == ["Trp53", "Myc"]
    assert converted_obs.obs_names.tolist() == ["Nfkb1", "unknown"]
    assert adata.var_names.tolist() == ["Tp53", "Myc"]
    assert adata.obs_names.tolist() == ["NF-kappaB", "unknown"]

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.omics.pp.convert_gene_identifiers(adata, axis=cast(Any, "bad"))


def test_standardize_gene_identifiers_is_deprecated(
    monkeypatch,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_genename, "create_identifiers", fake_gene_synonyms_cls)
    adata = ad.AnnData(
        X=np.ones((1, 1)),
        var=pd.DataFrame(index=["Tp53"]),
    )

    with pytest.warns(FutureWarning, match="bt.omics.pp.standardize_gene_identifiers"):
        result = getattr(bt.omics.pp, "standardize_gene_identifiers")(adata)

    assert result is None
    assert adata.var_names.tolist() == ["Trp53"]


def _dense_list(value: object) -> list:
    matrix = cast(Any, value)
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).tolist()


def test_convert_gene_identifiers_accepts_explicit_identifiers():
    class FakeGeneIdentifiers:
        def __init__(self):
            self.calls = []

        def __call__(self, df, **kwargs):
            self.calls.append((df, kwargs))
            df.index = ["Trp53", "Myc"]

    identifiers = FakeGeneIdentifiers()
    adata = ad.AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["cell"]),
        var=pd.DataFrame(index=["Tp53", "Myc"]),
    )

    result = bt.omics.pp.convert_gene_identifiers(
        adata,
        axis="var",
        identifiers=identifiers,
        copy=False,
    )

    assert result is None
    assert adata.var_names.tolist() == ["Trp53", "Myc"]
    assert len(identifiers.calls) == 1
    called_df, kwargs = identifiers.calls[0]
    assert called_df is adata.var
    assert kwargs == {
        "axis": "index",
        "input_type": "name",
        "output_type": "symbol",
        "copy": False,
    }


def test_convert_gene_identifiers_accepts_deprecated_identifier_arguments(
    fake_gene_synonyms_cls,
):
    adata = ad.AnnData(
        X=np.ones((1, 1)),
        var=pd.DataFrame(index=["Tp53"]),
    )

    with pytest.warns(FutureWarning) as warning_records:
        cast(Any, bt.omics.pp.convert_gene_identifiers)(
            adata,
            genesyn=fake_gene_synonyms_cls(),
            input_identifier_type="name",
            output_identifier_type="symbol",
        )

    messages = [str(record.message) for record in warning_records]
    assert any(
        "genesyn" in message and "identifiers" in message for message in messages
    )
    assert any("input_identifier_type" in message for message in messages)
    assert any("output_identifier_type" in message for message in messages)
    assert adata.var_names.tolist() == ["Trp53"]


def test_merge_duplicate_vars_sums_counts_and_var_rows():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=csr_matrix(
                np.array(
                    [
                        [1.0, 2.0, 3.0, 4.0],
                        [5.0, 6.0, 7.0, 8.0],
                    ]
                )
            ),
            obs=pd.DataFrame({"batch": ["a", "b"]}, index=["c1", "c2"]),
            var=pd.DataFrame(
                {
                    "symbol": ["alias", "g1", "x", "y"],
                    "tag": ["first", "preferred", "fallback", "ignored"],
                },
                index=["g1", "g1", "g2", "g2"],
            ),
        )

    with warnings.catch_warnings(record=True) as records:
        warnings.simplefilter("always")
        merged = bt.omics.pp.merge_duplicate_vars(adata)

    assert not any(
        "Variable names are not unique" in str(warning.message) for warning in records
    )
    var = cast(pd.DataFrame, merged.var)
    assert set(merged.var_names) == {"g1", "g2"}
    assert var.loc["g1", "tag"] == "first"
    assert var.loc["g2", "tag"] == "fallback"
    assert cast(pd.DataFrame, merged.obs).equals(cast(pd.DataFrame, adata.obs))
    assert adata.var_names.tolist() == ["g1", "g1", "g2", "g2"]

    ordered = merged[:, ["g1", "g2"]]
    assert _dense_list(ordered.X) == [[3.0, 7.0], [11.0, 15.0]]

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata_without_column = ad.AnnData(
            X=csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]])),
            obs=pd.DataFrame(index=["c1", "c2"]),
            var=pd.DataFrame(
                {"tag": ["first", "second"]},
                index=["g1", "g1"],
            ),
        )

    with warnings.catch_warnings(record=True) as records:
        warnings.simplefilter("always")
        merged_without_column = bt.omics.pp.merge_duplicate_vars(
            adata_without_column,
        )

    assert not any(
        "Variable names are not unique" in str(warning.message) for warning in records
    )
    assert merged_without_column.var_names.tolist() == ["g1"]
    assert "copy_var_names" not in merged_without_column.var
    assert _dense_list(merged_without_column.X) == [[3.0], [7.0]]


def test_merge_duplicate_vars_can_build_consensus_var_rows():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.ones((1, 3)),
            obs=pd.DataFrame(index=["c1"]),
            var=pd.DataFrame(
                {
                    "symbol": ["Aars", "Aars1", "Myc"],
                    "biotype": ["protein_coding", "protein_coding", "protein_coding"],
                },
                index=["Aars", "Aars", "Myc"],
            ),
        )

    merged = bt.omics.pp.merge_duplicate_vars(adata, keep="consensus")

    var = cast(pd.DataFrame, merged.var)
    assert merged.var_names.tolist() == ["Aars", "Myc"]
    assert pd.isna(var.loc["Aars", "symbol"])
    assert var.loc["Aars", "biotype"] == "protein_coding"
    assert var.loc["Myc", "symbol"] == "Myc"
    assert var.loc["Myc", "biotype"] == "protein_coding"


def test_merge_duplicate_vars_sums_layers():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            obs=pd.DataFrame(index=["c1", "c2"]),
            var=pd.DataFrame(index=["g1", "g1", "g2"]),
        )
    adata.layers["dense"] = np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]])
    adata.layers["sparse"] = csr_matrix(
        np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]])
    )

    merged = bt.omics.pp.merge_duplicate_vars(adata)

    assert _dense_list(merged.X) == [[3.0, 3.0], [9.0, 6.0]]
    assert _dense_list(merged.layers["dense"]) == [[30.0, 30.0], [90.0, 60.0]]
    assert _dense_list(merged.layers["sparse"]) == [
        [300.0, 300.0],
        [900.0, 600.0],
    ]


def test_merge_duplicate_vars_merges_varm_with_nan_strategy():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.ones((1, 3)),
            obs=pd.DataFrame(index=["c1"]),
            var=pd.DataFrame(index=["g1", "g1", "g2"]),
        )
        adata.varm["scores"] = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])

    merged = bt.omics.pp.merge_duplicate_vars(adata, varm="nan")

    assert np.isnan(merged.varm["scores"][0]).all()
    assert _dense_list(merged.varm["scores"][1]) == [3.0, 30.0]


def test_merge_duplicate_vars_merges_varm_with_first_strategy():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.ones((1, 3)),
            obs=pd.DataFrame(index=["c1"]),
            var=pd.DataFrame(
                {"symbol": ["alias", "g1", "g2"]},
                index=["g1", "g1", "g2"],
            ),
        )
        adata.varm["scores"] = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])
        adata.varm["df_scores"] = pd.DataFrame(
            {"score": [1.0, 2.0, 3.0]},
            index=["g1", "g1", "g2"],
        )

    merged = bt.omics.pp.merge_duplicate_vars(adata, varm="first")

    assert _dense_list(merged.varm["scores"]) == [[1.0, 10.0], [3.0, 30.0]]
    df_scores = merged.varm["df_scores"]
    assert isinstance(df_scores, pd.DataFrame)
    assert df_scores.index.tolist() == ["g1", "g2"]
    assert df_scores["score"].tolist() == [1.0, 3.0]


def test_merge_duplicate_vars_merges_varm_with_mean_strategy():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.ones((1, 3)),
            obs=pd.DataFrame(index=["c1"]),
            var=pd.DataFrame(index=["g1", "g1", "g2"]),
        )
        adata.varm["scores"] = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])

    merged = bt.omics.pp.merge_duplicate_vars(adata, varm="mean")

    assert _dense_list(merged.varm["scores"]) == [[1.5, 15.0], [3.0, 30.0]]


def test_merge_duplicate_vars_rejects_varp_by_default():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=np.ones((1, 2)),
            obs=pd.DataFrame(index=["c1"]),
            var=pd.DataFrame(index=["g1", "g1"]),
        )
        adata.varp["correlation"] = np.ones((2, 2))

    with pytest.raises(NotImplementedError, match="does not merge `.varp`"):
        bt.omics.pp.merge_duplicate_vars(adata)

    merged = bt.omics.pp.merge_duplicate_vars(adata, varp="drop")

    assert len(merged.varp) == 0
    assert _dense_list(merged.X) == [[2.0]]


def test_merge_duplicate_vars_can_modify_in_place():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        adata = ad.AnnData(
            X=csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]])),
            obs=pd.DataFrame(index=["c1", "c2"]),
            var=pd.DataFrame(
                {"symbol": ["alias", "g1"], "tag": ["first", "preferred"]},
                index=["g1", "g1"],
            ),
        )

    result = bt.omics.pp.merge_duplicate_vars(
        adata,
        copy=False,
    )

    assert result is None
    var = cast(pd.DataFrame, adata.var)
    assert adata.var_names.tolist() == ["g1"]
    assert var.loc["g1", "tag"] == "first"
    assert _dense_list(adata.X) == [[3.0], [7.0]]
