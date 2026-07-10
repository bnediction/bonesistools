#!/usr/bin/env python

import warnings
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom

import bonesistools as bt
from bonesistools.omics.tools._markers import DEA_METHODS


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


def _toy_dea_adata():
    adata = ad.AnnData(
        X=np.array(
            [
                [9.0, 1.0, 2.0, 8.0, 0.0, 0.0],
                [10.0, 1.0, 2.0, 8.0, 0.0, 0.0],
                [11.0, 2.0, 2.0, 8.0, 0.0, 0.0],
                [12.0, 2.0, 2.0, 8.0, 0.0, 0.0],
                [1.0, 9.0, 2.0, 2.0, 1.0, 0.0],
                [1.0, 10.0, 2.0, 2.0, 1.0, 0.0],
                [2.0, 11.0, 2.0, 2.0, 1.0, 0.0],
                [2.0, 12.0, 2.0, 2.0, 1.0, 0.0],
                [4.0, 4.0, 2.0, 20.0, 0.0, 0.0],
                [4.0, 4.0, 2.0, 20.0, 0.0, 0.0],
                [5.0, 5.0, 2.0, 20.0, 0.0, 0.0],
                [5.0, 5.0, 2.0, 20.0, 0.0, 0.0],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4 + ["C"] * 4)},
            index=[f"c{i}" for i in range(12)],
        ),
        var=pd.DataFrame(index=["g1", "g2", "g3", "g4", "g5", "g6"]),
    )
    adata.layers["counts"] = csr_matrix(adata.X)
    return adata


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

    result = bt.omics.tl.logfoldchanges(
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

    result = bt.omics.tl.logfoldchanges(
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
    result = bt.omics.tl.logfoldchanges(
        _empty_group_logfoldchange_adata(),
        groupby="cluster",
    )

    assert result.empty
    assert result.columns.tolist() == ["group", "names", "logfoldchanges"]


def test_logfoldchanges_handles_zero_means_without_warning():
    with warnings.catch_warnings(record=True) as warnings_record:
        warnings.simplefilter("always")
        result = bt.omics.tl.logfoldchanges(
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

    filtered = bt.omics.tl.logfoldchanges(
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
        bt.omics.tl.logfoldchanges(
            mini_adata,
            groupby="cluster",
            filter_logfoldchanges=cast(Any, "not callable"),
        )


@pytest.mark.parametrize(
    "method",
    ["welch", "welch_overestimate", "wilcoxon"],
)
def test_dea_matches_underlying_statistical_tests(method):
    adata = _toy_dea_adata()

    result = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method=method,
        correction=None,
        alpha=None,
    )
    if method == "wilcoxon":
        statistics = bt.omics.tl.wilcoxon_tests(
            adata,
            groupby="cluster",
            groups=["A"],
            background=["B"],
            correction=None,
        )
    else:
        statistics = bt.omics.tl.welch_tests(
            adata,
            groupby="cluster",
            groups=["A"],
            background=["B"],
            correction=None,
            overestimate_variance=method == "welch_overestimate",
        )

    statistics = statistics.reset_index().rename(columns={"names": "feature"})
    statistics = statistics[["feature", "group", "statistics", "pvals", "pvals_adj"]]
    pd.testing.assert_frame_equal(
        result.loc[:, statistics.columns],
        statistics,
        check_dtype=False,
    )


def test_dea_returns_expected_logfoldchanges_with_explicit_background():
    adata = _toy_dea_adata()

    result = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        correction=None,
        alpha=None,
    )
    observed = result.set_index("feature")["logfoldchanges"].to_dict()

    assert observed["g1"] == pytest.approx(np.log2(10.5 / 1.5))
    assert observed["g2"] == pytest.approx(np.log2(1.5 / 10.5))
    assert observed["g3"] == pytest.approx(0.0)
    assert observed["g4"] == pytest.approx(2.0)
    assert observed["g5"] == -np.inf
    assert np.isnan(observed["g6"])


def test_dea_logfoldchanges_support_logged_expression():
    adata = _toy_dea_adata()
    expression_mtx = cast(np.ndarray, adata.X)
    adata.layers["log"] = np.log1p(expression_mtx)

    counts = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        correction=None,
        alpha=None,
    ).set_index("feature")
    logged = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        expression="log",
        is_log=True,
        correction=None,
        alpha=None,
    ).set_index("feature")

    pd.testing.assert_series_equal(
        logged["logfoldchanges"],
        counts["logfoldchanges"],
        check_names=False,
    )


def test_dea_background_rest_and_named_background_are_distinct():
    adata = _toy_dea_adata()

    explicit = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        correction=None,
        alpha=None,
    ).set_index("feature")
    rest = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background="rest",
        method="wilcoxon",
        correction=None,
        alpha=None,
    ).set_index("feature")

    assert explicit.loc["g4", "logfoldchanges"] == pytest.approx(2.0)
    assert rest.loc["g4", "logfoldchanges"] == pytest.approx(np.log2(8.0 / 11.0))


def test_dea_alpha_filters_adjusted_pvalues():
    adata = _toy_dea_adata()

    complete = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="welch",
        alpha=None,
    )
    filtered = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="welch",
        alpha=0.05,
    )
    expected = complete.loc[complete["pvals_adj"] <= 0.05].reset_index(drop=True)

    assert len(filtered) < len(complete)
    pd.testing.assert_frame_equal(filtered, expected)


def test_dea_filters_logfoldchanges_after_pvalue_filtering():
    adata = _toy_dea_adata()

    complete = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        correction=None,
        alpha=None,
    )
    filtered = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="wilcoxon",
        correction=None,
        alpha=None,
        filter_logfoldchanges=lambda values: np.abs(values) >= 1,
    )
    expected = complete.loc[
        np.abs(complete["logfoldchanges"].to_numpy(dtype=float)) >= 1
    ].reset_index(drop=True)

    pd.testing.assert_frame_equal(filtered, expected)


def test_dea_preserves_infinite_logfoldchanges():
    result = bt.omics.tl.dea(
        _zero_mean_logfoldchange_adata(),
        groupby="cluster",
        method="wilcoxon",
        correction=None,
        alpha=None,
    )
    observed = result.set_index(["group", "feature"])["logfoldchanges"].to_dict()

    assert observed[("A", "g1")] == -np.inf
    assert observed[("B", "g1")] == np.inf


def test_dea_var_subset_and_output_columns():
    adata = _toy_dea_adata()
    adata.var["selected"] = [True, False, False, True, False, False]

    result = bt.omics.tl.dea(
        adata,
        groupby="cluster",
        groups=["A"],
        background=["B"],
        method="welch",
        var_subset="selected",
        alpha=None,
    )

    assert result["feature"].tolist() == ["g4", "g1"]
    assert result.columns.tolist() == [
        "feature",
        "group",
        "statistics",
        "pvals",
        "pvals_adj",
        "logfoldchanges",
    ]


def test_dea_output_schema_is_identical_for_all_methods():
    adata = _toy_dea_adata()
    expected_columns = [
        "feature",
        "group",
        "statistics",
        "pvals",
        "pvals_adj",
        "logfoldchanges",
    ]

    for method in DEA_METHODS:
        result = bt.omics.tl.dea(
            adata,
            groupby="cluster",
            groups=["A"],
            background=["B"],
            method=method,
            alpha=None,
        )

        assert result.columns.tolist() == expected_columns


def test_dea_rejects_invalid_method(mini_adata):
    with pytest.raises(ValueError):
        bt.omics.tl.dea(
            mini_adata,
            groupby="cluster",
            method=cast(Any, "bad"),
        )


def test_ora_returns_expected_hypergeometric_probability():
    result = bt.omics.tl.ora(
        query_set=["g1", "g2", "g5"],
        signatures={"sig": ["g1", "g2", "g3", "g4"]},
        background=[f"g{i}" for i in range(1, 11)],
        correction="bonferroni",
        include_overlap=True,
    )

    expected = hypergeom.sf(2 - 1, 10, 4, 3)
    assert result.index.name == "signature"
    assert result.index.tolist() == ["sig"]
    assert result.loc["sig", "pvals"] == pytest.approx(expected)
    assert result.loc["sig", "pvals_adj"] == pytest.approx(expected)
    assert result.loc["sig", "observed_overlap"] == 2
    assert result.attrs["query_size"] == 3
    assert result.attrs["background_size"] == 10
    assert "query_size" not in result.columns
    assert "background_size" not in result.columns
    assert result.loc["sig", "expected_overlap"] == pytest.approx(1.2)
    assert result.loc["sig", "fold_enrichment"] == pytest.approx(2 / 1.2)
    assert result.loc["sig", "signature_size"] == 4
    assert result.loc["sig", "overlap"] == ("g1", "g2")
    assert result.columns.tolist() == [
        "pvals",
        "pvals_adj",
        "observed_overlap",
        "expected_overlap",
        "fold_enrichment",
        "signature_size",
        "overlap",
    ]


def test_ora_hides_overlap_by_default():
    result = bt.omics.tl.ora(
        query_set=["g1", "g2"],
        signatures={"sig": ["g1", "g3"]},
        background=["g1", "g2", "g3"],
    )

    assert "overlap" not in result.columns


def test_ora_filters_query_and_signatures_by_background():
    result = bt.omics.tl.ora(
        query_set=["g1", "g2", "missing_query"],
        signatures={
            "tested": ["g1", "missing_signature"],
            "skipped": ["missing_only"],
        },
        background=["g1", "g2", "g3"],
        correction="bonferroni",
    )

    assert result.index.tolist() == ["tested"]
    assert result.attrs["query_size"] == 2
    assert result.loc["tested", "signature_size"] == 1
    assert result.loc["tested", "observed_overlap"] == 1
    assert result.loc["tested", "pvals"] == pytest.approx(hypergeom.sf(0, 3, 1, 2))


def test_ora_rejects_empty_background_and_filtered_query():
    with pytest.raises(ValueError):
        bt.omics.tl.ora(
            query_set=["g1"],
            signatures={"sig": ["g1"]},
            background=[],
        )

    with pytest.raises(ValueError):
        bt.omics.tl.ora(
            query_set=["missing"],
            signatures={"sig": ["g1"]},
            background=["g1"],
        )


def test_ora_rejects_invalid_gene_collections_and_signatures():
    with pytest.raises(TypeError):
        bt.omics.tl.ora(
            query_set=cast(Any, "g1"),
            signatures={"sig": ["g1"]},
            background=["g1"],
        )

    with pytest.raises(TypeError):
        bt.omics.tl.ora(
            query_set=["g1"],
            signatures={"sig": ["g1"]},
            background=cast(Any, ["g1", 1]),
        )

    with pytest.raises(ValueError):
        bt.omics.tl.ora(
            query_set=["g1"],
            signatures=[("sig", ["g1"]), ("sig", ["g2"])],
            background=["g1", "g2"],
        )

    with pytest.raises(TypeError):
        bt.omics.tl.ora(
            query_set=["g1"],
            signatures=cast(Any, [["g1"]]),
            background=["g1"],
        )


def test_ora_mapping_and_sequence_signatures_are_equivalent():
    mapping = bt.omics.tl.ora(
        query_set=["g1", "g2"],
        signatures={"A": ["g1", "g3"], "B": ["g2"]},
        background=["g1", "g2", "g3", "g4"],
        correction="bonferroni",
    )
    sequence = bt.omics.tl.ora(
        query_set=["g1", "g2"],
        signatures=[("A", ["g1", "g3"]), ("B", ["g2"])],
        background=["g1", "g2", "g3", "g4"],
        correction="bonferroni",
    )

    pd.testing.assert_frame_equal(mapping, sequence)


def test_ora_correction_methods_are_explicit():
    kwargs = {
        "query_set": ["g1", "g2"],
        "signatures": {
            "A": ["g1"],
            "B": ["g1", "g2"],
            "C": ["g3"],
        },
        "background": [f"g{i}" for i in range(1, 11)],
    }

    bonferroni = bt.omics.tl.ora(**kwargs, correction="bonferroni")
    bh = bt.omics.tl.ora(**kwargs, correction="benjamini-hochberg")

    assert bonferroni["pvals_adj"].tolist() == pytest.approx(
        np.minimum(bonferroni["pvals"].to_numpy(dtype=float) * 3, 1.0)
    )
    assert bh["pvals_adj"].to_dict() == pytest.approx(
        {
            "B": 1 / 15,
            "A": 0.3,
            "C": 1.0,
        }
    )

    with pytest.raises(TypeError):
        bt.omics.tl.ora(**kwargs, correction=cast(Any, None))


def test_ora_sorting_uses_observed_overlap_to_break_pvalue_ties(monkeypatch):
    from bonesistools.omics.tools import _markers

    monkeypatch.setattr(
        _markers,
        "_hypergeometric_pvalue",
        lambda **__: 0.5,
    )

    result = bt.omics.tl.ora(
        query_set=["g1", "g2"],
        signatures={
            "small_overlap": ["g1"],
            "large_overlap": ["g1", "g2"],
        },
        background=["g1", "g2", "g3"],
        correction="bonferroni",
    )

    assert result.index.tolist() == ["large_overlap", "small_overlap"]


def test_smirnov_tests_recovers_expected_ks_statistics():
    adata = _toy_smirnov_adata()

    bt.omics.tl.smirnov_tests(
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

    result = bt.omics.tl.smirnov_tests(
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
    result = bt.omics.tl.smirnov_tests(
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
        bt.omics.tl.smirnov_tests(
            mini_adata,
            groupby="cluster",
            groups=cast(Any, 1),
        )

    with pytest.raises(ValueError, match="invalid argument value for 'reference'"):
        bt.omics.tl.smirnov_tests(
            mini_adata,
            groupby="cluster",
            reference="missing",
        )


def test_smirnov_tests_copy_returns_dataframe(mini_adata):
    result = bt.omics.tl.smirnov_tests(
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
