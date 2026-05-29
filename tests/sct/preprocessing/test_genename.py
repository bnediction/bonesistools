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
    monkeypatch.setattr(_genename, "GeneSynonyms", fake_gene_synonyms_cls)

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


def test_var_names_merge_duplicates_sums_counts_and_var_rows():
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
            obs=pd.DataFrame(index=["c1", "c2"]),
            var=pd.DataFrame(
                {
                    "symbol": ["alias", "g1", "x", "y"],
                    "tag": ["first", "preferred", "fallback", "ignored"],
                },
                index=["g1", "g1", "g2", "g2"],
            ),
        )

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        merged = bt.sct.pp.var_names_merge_duplicates(
            adata,
            var_names_column="symbol",
        )

    assert set(merged.var_names) == {"g1", "g2"}
    assert merged.var.loc["g1", "tag"] == "preferred"
    assert merged.var.loc["g2", "tag"] == "fallback"

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

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        merged_without_column = bt.sct.pp.var_names_merge_duplicates(
            adata_without_column,
        )

    assert merged_without_column.var_names.tolist() == ["g1"]
    assert "copy_var_names" not in merged_without_column.var
    assert merged_without_column.X.toarray().tolist() == [[3.0], [7.0]]
