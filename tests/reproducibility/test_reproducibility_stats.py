#!/usr/bin/env python

import os

import anndata as ad
import pandas as pd
import pytest

import bonesistools as bt
from tests.regression.omics.toy_data import make_nestorowa_hvg_adata

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def _load_nestorowa_hvg() -> ad.AnnData:

    return make_nestorowa_hvg_adata()


def test_wilcoxon_tests_are_reproducible_on_nestorowa():

    adata = _load_nestorowa_hvg()

    whole = bt.omics.tl.wilcoxon_tests(
        adata,
        groupby="label",
        groups="all",
        correction="benjamini-hochberg",
    )
    memory_limited = bt.omics.tl.wilcoxon_tests(
        adata,
        groupby="label",
        groups="all",
        correction="benjamini-hochberg",
        max_memory="1MB",
    )

    pd.testing.assert_frame_equal(memory_limited, whole)
