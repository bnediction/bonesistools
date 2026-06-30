#!/usr/bin/env python

import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

import bonesistools as bt

sys.path.insert(0, str(Path(__file__).parents[1] / "sct"))
from toy_data import make_nestorowa_hvg_adata  # noqa: E402

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def _load_nestorowa_hvg() -> ad.AnnData:

    return make_nestorowa_hvg_adata()


def test_wilcoxon_tests_are_reproducible_on_nestorowa():

    adata = _load_nestorowa_hvg()

    whole = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="label",
        groups="all",
        correction="benjamini-hochberg",
    )
    memory_limited = bt.sct.tl.wilcoxon_tests(
        adata,
        groupby="label",
        groups="all",
        correction="benjamini-hochberg",
        max_memory="1MB",
    )

    pd.testing.assert_frame_equal(memory_limited, whole)
