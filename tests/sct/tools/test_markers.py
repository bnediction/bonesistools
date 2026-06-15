#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt


def _toy_logfoldchange_adata():
    return ad.AnnData(
        X=np.array(
            [
                [4.0, 2.0, 1.0],
                [4.0, 2.0, 1.0],
                [1.0, 4.0, 1.0],
                [1.0, 4.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B"])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


def _toy_rebalanced_logfoldchange_adata():
    return ad.AnnData(
        X=np.array(
            [
                [8.0],
                [8.0],
                [2.0],
                [4.0],
                [4.0],
                [4.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": ["A", "A", "B", "C", "C", "C"]},
            index=["c1", "c2", "c3", "c4", "c5", "c6"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )


def _toy_smirnov_adata():
    return ad.AnnData(
        X=np.array(
            [
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B"])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["separated", "identical"]),
    )


def _empty_group_logfoldchange_adata():
    return ad.AnnData(
        X=np.array([[1.0], [2.0]]),
        obs=pd.DataFrame({"cluster": [np.nan, np.nan]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1"]),
    )


def _zero_mean_logfoldchange_adata():
    return ad.AnnData(
        X=np.array(
            [
                [0.0],
                [0.0],
                [2.0],
                [2.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "A", "B", "B"])},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )


def _empty_var_smirnov_adata():
    return ad.AnnData(
        X=np.empty((2, 0)),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A", "B"])},
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=[]),
    )


def test_logfoldchanges_recovers_expected_deseq2_style_contrasts():
    adata = _toy_logfoldchange_adata()

    result = bt.sct.tl.logfoldchanges(
        adata,
        groupby="cluster",
    )
    observed = result.set_index(["group", "names"])["logfoldchanges"].to_dict()

    # Simple two-group DESeq2-style contrast with size factors equal to 1:
    # log2(mean normalized counts in group / mean normalized counts in reference).
    assert observed == pytest.approx(
        {
            ("A", "g1"): 2.0,
            ("A", "g2"): -1.0,
            ("A", "g3"): 0.0,
            ("B", "g1"): -2.0,
            ("B", "g2"): 1.0,
            ("B", "g3"): 0.0,
        }
    )


def test_logfoldchanges_rebalances_cluster_means():
    adata = _toy_rebalanced_logfoldchange_adata()

    result = bt.sct.tl.logfoldchanges(
        adata,
        groupby="cluster",
        cluster_rebalancing=True,
    )
    observed = result.set_index(["group", "names"])["logfoldchanges"].to_dict()

    # Cluster means are A=8, B=2, C=4. With rebalancing, the reference for each
    # group is the mean of the other cluster means, not the mean of all cells.
    assert observed == pytest.approx(
        {
            ("A", "g1"): np.log2(8 / 3),
            ("B", "g1"): np.log2(2 / 6),
            ("C", "g1"): np.log2(4 / 5),
        }
    )


def test_logfoldchanges_returns_empty_table_without_groups():
    result = bt.sct.tl.logfoldchanges(
        _empty_group_logfoldchange_adata(),
        groupby="cluster",
    )

    assert result.empty
    assert result.columns.tolist() == ["group", "names", "logfoldchanges"]


def test_logfoldchanges_handles_zero_means_without_warning():
    with warnings.catch_warnings(record=True) as warnings_record:
        warnings.simplefilter("always")
        result = bt.sct.tl.logfoldchanges(
            _zero_mean_logfoldchange_adata(),
            groupby="cluster",
        )

    messages = [str(record.message) for record in warnings_record]
    assert not any("divide by zero" in message for message in messages)

    observed = result.set_index(["group", "names"])["logfoldchanges"].to_dict()
    assert observed[("A", "g1")] == -np.inf
    assert observed[("B", "g1")] == np.inf


def test_logfoldchanges_filtering_keeps_expected_positive_ratios():
    adata = _toy_logfoldchange_adata()

    filtered = bt.sct.tl.logfoldchanges(
        adata,
        groupby="cluster",
        filter_logfoldchanges=lambda values: values > 0,
    )
    observed = filtered.set_index(["group", "names"])["logfoldchanges"].to_dict()

    assert observed == pytest.approx(
        {
            ("A", "g1"): 2.0,
            ("B", "g2"): 1.0,
        }
    )


def test_logfoldchanges_rejects_invalid_filter(mini_adata):
    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'filter_logfoldchanges'",
    ):
        bt.sct.tl.logfoldchanges(
            mini_adata,
            groupby="cluster",
            filter_logfoldchanges=cast(Any, "not callable"),
        )


def test_hypergeometric_test_returns_expected_probability():
    adata = ad.AnnData(
        X=np.ones((1, 10)),
        var=pd.DataFrame(index=[f"g{i}" for i in range(1, 11)]),
    )

    pvalue = bt.sct.tl.hypergeometric_test(
        adata,
        signature=["g1", "g2", "g3", "g4"],
        markers=["g1", "g2", "g5"],
    )

    # P(X >= 2) for X ~ Hypergeometric(M=10, n=4, N=3):
    # (C(4, 2) C(6, 1) + C(4, 3) C(6, 0)) / C(10, 3) = 1/3.
    assert pvalue == pytest.approx(1 / 3)

    with pytest.raises(TypeError, match="unsupported argument type for 'adata'"):
        bt.sct.tl.hypergeometric_test(
            cast(Any, object()),
            signature=["g1"],
            markers=["g1"],
        )


def test_smirnov_tests_recovers_expected_ks_statistics():
    adata = _toy_smirnov_adata()

    bt.sct.tl.smirnov_tests(
        adata,
        groupby="cluster",
        groups=["A"],
        reference="B",
        corr_method="bonferroni",
        key_added="ks",
    )

    result = adata.uns["ks"]["results"].set_index("names")
    assert result["group"].unique().tolist() == ["A"]
    assert result.loc["separated", "statistics"] == pytest.approx(1.0)
    assert result.loc["separated", "pvals"] == pytest.approx(1 / 3)
    assert result.loc["separated", "pvals_adj"] == pytest.approx(2 / 3)
    assert result.loc["identical", "statistics"] == pytest.approx(0.0)
    assert result.loc["identical", "pvals"] == pytest.approx(1.0)
    assert result.loc["identical", "pvals_adj"] == pytest.approx(1.0)


def test_smirnov_tests_rest_reference_sparse_layer_and_bh_cutoff():
    adata = _toy_smirnov_adata()
    adata.layers["sparse"] = csr_matrix(adata.X)

    result = bt.sct.tl.smirnov_tests(
        adata,
        groupby="cluster",
        layer="sparse",
        pval_cutoff=0.8,
        copy=True,
    )

    result = result.set_index(["group", "names"])

    assert sorted(result.index.tolist()) == [
        ("A", "separated"),
        ("B", "separated"),
    ]
    assert result["statistics"].to_dict() == pytest.approx(
        {
            ("A", "separated"): 1.0,
            ("B", "separated"): 1.0,
        }
    )
    assert result["pvals"].to_dict() == pytest.approx(
        {
            ("A", "separated"): 1 / 3,
            ("B", "separated"): 1 / 3,
        }
    )
    assert result["pvals_adj"].to_dict() == pytest.approx(
        {
            ("A", "separated"): 2 / 3,
            ("B", "separated"): 2 / 3,
        }
    )


def test_smirnov_tests_empty_gene_table_has_empty_adjusted_pvalues():
    result = bt.sct.tl.smirnov_tests(
        _empty_var_smirnov_adata(),
        groupby="cluster",
        copy=True,
    )

    assert result.empty
    assert result.columns.tolist() == [
        "group",
        "names",
        "statistics",
        "locations",
        "signs",
        "pvals",
        "pvals_adj",
    ]


def test_smirnov_tests_validates_groups_and_reference(mini_adata):
    with pytest.raises(TypeError, match="unsupported argument type for 'groups'"):
        bt.sct.tl.smirnov_tests(
            mini_adata,
            groupby="cluster",
            groups=cast(Any, 1),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'reference'"):
        bt.sct.tl.smirnov_tests(
            mini_adata,
            groupby="cluster",
            reference="missing",
        )


def test_smirnov_tests_copy_returns_dataframe(mini_adata):
    result = bt.sct.tl.smirnov_tests(
        mini_adata,
        groupby="cluster",
        groups=["A"],
        reference="B",
        corr_method="bonferroni",
        copy=True,
    )

    assert isinstance(result, pd.DataFrame)
    assert "pvals_adj" in result
    assert "smirnov_tests" not in mini_adata.uns
