#!/usr/bin/env python

import os

import pandas as pd
import pytest

import bonesistools as bt

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def test_wilcoxon_tests_are_reproducible_on_nestorowa():

    adata = bt.sct.datasets.nestorowa()

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
