#!/usr/bin/env python

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.preprocessing import _genename


def test_convert_and_standardize_gene_identifiers_copy_and_axis_validation(
    monkeypatch,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_genename, "create_gene_synonyms", fake_gene_synonyms_cls)

    adata = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(index=["NF-kappaB", "unknown"]),
        var=pd.DataFrame(index=["Tp53", "Myc"]),
    )

    converted_var = bt.sct.pp.convert_gene_identifiers(adata, axis="var", copy=True)
    converted_obs = bt.sct.pp.standardize_gene_identifiers(
        adata,
        axis="obs",
        copy=True,
    )

    assert converted_var.var_names.tolist() == ["Trp53", "Myc"]
    assert converted_obs.obs_names.tolist() == ["Nfkb1", "unknown"]
    assert adata.var_names.tolist() == ["Tp53", "Myc"]
    assert adata.obs_names.tolist() == ["NF-kappaB", "unknown"]

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.pp.convert_gene_identifiers(adata, axis="bad")


def test_convert_gene_identifiers_accepts_explicit_gene_synonyms():
    class FakeGeneSynonyms:
        def __init__(self):
            self.calls = []

        def __call__(self, df, **kwargs):
            self.calls.append((df, kwargs))
            df.index = ["Trp53", "Myc"]

    genesyn = FakeGeneSynonyms()
    adata = ad.AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["cell"]),
        var=pd.DataFrame(index=["Tp53", "Myc"]),
    )

    result = bt.sct.pp.convert_gene_identifiers(
        adata,
        axis="var",
        genesyn=genesyn,
        copy=False,
    )

    assert result is None
    assert adata.var_names.tolist() == ["Trp53", "Myc"]
    assert len(genesyn.calls) == 1
    called_df, kwargs = genesyn.calls[0]
    assert called_df is adata.var
    assert kwargs == {
        "axis": "index",
        "input_identifier_type": "name",
        "output_identifier_type": "official_name",
        "copy": False,
    }


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
        merged = bt.sct.pp.merge_duplicate_vars(adata)

    assert not any(
        "Variable names are not unique" in str(warning.message) for warning in records
    )
    assert set(merged.var_names) == {"g1", "g2"}
    assert merged.var.loc["g1", "tag"] == "first"
    assert merged.var.loc["g2", "tag"] == "fallback"
    assert merged.obs.equals(adata.obs)
    assert adata.var_names.tolist() == ["g1", "g1", "g2", "g2"]

    ordered = merged[:, ["g1", "g2"]]
    assert ordered.X.toarray().tolist() == [[3.0, 7.0], [11.0, 15.0]]

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
        merged_without_column = bt.sct.pp.merge_duplicate_vars(
            adata_without_column,
        )

    assert not any(
        "Variable names are not unique" in str(warning.message) for warning in records
    )
    assert merged_without_column.var_names.tolist() == ["g1"]
    assert "copy_var_names" not in merged_without_column.var
    assert merged_without_column.X.toarray().tolist() == [[3.0], [7.0]]


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

    merged = bt.sct.pp.merge_duplicate_vars(adata, keep="consensus")

    assert merged.var_names.tolist() == ["Aars", "Myc"]
    assert np.isnan(merged.var.loc["Aars", "symbol"])
    assert merged.var.loc["Aars", "biotype"] == "protein_coding"
    assert merged.var.loc["Myc", "symbol"] == "Myc"
    assert merged.var.loc["Myc", "biotype"] == "protein_coding"


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

    merged = bt.sct.pp.merge_duplicate_vars(adata)

    assert merged.X.tolist() == [[3.0, 3.0], [9.0, 6.0]]
    assert merged.layers["dense"].tolist() == [[30.0, 30.0], [90.0, 60.0]]
    assert merged.layers["sparse"].toarray().tolist() == [
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
            varm={"scores": np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])},
        )

    merged = bt.sct.pp.merge_duplicate_vars(adata, varm="nan")

    assert np.isnan(merged.varm["scores"][0]).all()
    assert merged.varm["scores"][1].tolist() == [3.0, 30.0]


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
            varm={
                "scores": np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]]),
                "df_scores": pd.DataFrame(
                    {"score": [1.0, 2.0, 3.0]},
                    index=["g1", "g1", "g2"],
                ),
            },
        )

    merged = bt.sct.pp.merge_duplicate_vars(adata, varm="first")

    assert merged.varm["scores"].tolist() == [[1.0, 10.0], [3.0, 30.0]]
    assert isinstance(merged.varm["df_scores"], pd.DataFrame)
    assert merged.varm["df_scores"].index.tolist() == ["g1", "g2"]
    assert merged.varm["df_scores"]["score"].tolist() == [1.0, 3.0]


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
            varm={"scores": np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])},
        )

    merged = bt.sct.pp.merge_duplicate_vars(adata, varm="mean")

    assert merged.varm["scores"].tolist() == [[1.5, 15.0], [3.0, 30.0]]


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
            varp={"correlation": np.ones((2, 2))},
        )

    with pytest.raises(NotImplementedError, match="does not merge `.varp`"):
        bt.sct.pp.merge_duplicate_vars(adata)

    merged = bt.sct.pp.merge_duplicate_vars(adata, varp="drop")

    assert len(merged.varp) == 0
    assert merged.X.tolist() == [[2.0]]


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

    result = bt.sct.pp.merge_duplicate_vars(
        adata,
        copy=False,
    )

    assert result is None
    assert adata.var_names.tolist() == ["g1"]
    assert adata.var.loc["g1", "tag"] == "first"
    assert adata.X.toarray().tolist() == [[3.0], [7.0]]
