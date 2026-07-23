#!/usr/bin/env python

import os

import numpy as np
import pytest
from llvmlite import binding

from ._umap_diagnostics import load_expected, run_umap_diagnostics

_RAW_SPECTRAL_CHECKPOINT = "raw_spectral_layout"


@pytest.fixture(scope="module")
def umap_diagnostics():

    return run_umap_diagnostics()


def test_golden_umap_numerical_checkpoints(umap_diagnostics):
    expected = load_expected()
    assert tuple(umap_diagnostics) == tuple(expected)

    context = (
        f"host_cpu={binding.get_host_cpu_name()}, "
        f"host_fma={bool(binding.get_host_cpu_features().get('fma'))}, "
        f"NUMBA_CPU_NAME={os.environ.get('NUMBA_CPU_NAME')!r}, "
        f"NUMBA_CPU_FEATURES={os.environ.get('NUMBA_CPU_FEATURES')!r}"
    )
    for checkpoint, expected_values in expected.items():
        if checkpoint == _RAW_SPECTRAL_CHECKPOINT:
            np.testing.assert_allclose(
                umap_diagnostics[checkpoint],
                expected_values,
                rtol=1e-12,
                atol=1e-13,
                err_msg=(
                    f"first divergent UMAP checkpoint: {checkpoint!r} ({context})"
                ),
            )
            continue

        np.testing.assert_array_equal(
            umap_diagnostics[checkpoint],
            expected_values,
            err_msg=f"first divergent UMAP checkpoint: {checkpoint!r} ({context})",
        )
