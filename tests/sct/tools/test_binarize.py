#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import bonesistools as bt


def _toy_binarizer_adata():
    adata = ad.AnnData(
        X=np.array(
            [
                [8.0, 1.0, 2.0],
                [9.0, 1.0, 2.0],
                [10.0, 1.0, 2.0],
                [1.0, 8.0, 2.0],
                [1.0, 9.0, 2.0],
                [1.0, 10.0, 2.0],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"macrostate": pd.Categorical(["A", "A", "A", "B", "B", "B"])},
            index=["c1", "c2", "c3", "c4", "c5", "c6"],
        ),
        var=pd.DataFrame(
            {"selected": [True, True, False]},
            index=["up_in_A", "up_in_B", "unchanged"],
        ),
    )
    adata.layers["log"] = np.log1p(adata.X)
    return adata


def test_dea_binarizer_keeps_complete_dea_before_binarization_filtering():
    adata = _toy_binarizer_adata()
    binarizer = bt.sct.tl.DEABinarizer(
        method="wilcoxon",
        correction=None,
        alpha=1.0,
        min_abs_logfoldchange=2.0,
    )

    abstraction = binarizer.fit_binarize(
        adata,
        obs="macrostate",
        expression="log",
        is_log=True,
        var_subset="selected",
    )

    assert binarizer.groups_.tolist() == ["A", "B"]
    assert binarizer.features_.tolist() == ["up_in_A", "up_in_B"]
    assert sorted(binarizer.dea_["feature"].unique()) == ["up_in_A", "up_in_B"]
    assert len(binarizer.dea_) == 4
    assert binarizer.dea_["logfoldchanges"].abs().min() >= 2.0

    expected = pd.DataFrame(
        [[1.0, 0.0], [0.0, 1.0]],
        index=pd.Index(["A", "B"]),
        columns=pd.Index(["up_in_A", "up_in_B"]),
    )
    pd.testing.assert_frame_equal(abstraction, expected)
    pd.testing.assert_frame_equal(binarizer.abstraction_, expected)


def test_dea_binarizer_applies_logfoldchange_threshold_only_in_binarize():
    adata = _toy_binarizer_adata()
    binarizer = bt.sct.tl.DEABinarizer(
        method="wilcoxon",
        correction=None,
        alpha=1.0,
        min_abs_logfoldchange=10.0,
    )

    binarizer.fit(adata, obs="macrostate")
    abstraction = binarizer.binarize()

    assert len(binarizer.dea_) == 6
    assert binarizer.dea_["logfoldchanges"].abs().max() < 10.0
    assert abstraction.isna().all(axis=None)


def test_dea_binarizer_requires_fit_before_binarize():
    with pytest.raises(RuntimeError):
        bt.sct.tl.DEABinarizer().binarize()
