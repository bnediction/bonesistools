#!/usr/bin/env python

import pytest

from ._omics_workflow import (
    assert_array_keys,
    assert_close_arrays,
    assert_equal_arrays,
    available_expected_names,
    load_expected,
    run_omics_workflow,
)


@pytest.fixture(scope="module")
def golden_outputs():

    return run_omics_workflow()


@pytest.mark.parametrize("name", available_expected_names())
def test_golden_omics_workflow_matches_expected_outputs(golden_outputs, name):
    result = golden_outputs[name]
    expected = load_expected(name)

    assert_array_keys(result, expected)

    if name == "qc":
        assert_equal_arrays(result, expected, ("obs_columns", "var_columns"))
        assert_close_arrays(
            result,
            expected,
            ("obs_values", "var_values"),
            rtol=1e-10,
            atol=1e-12,
        )
    elif name.startswith("hvg_"):
        assert_equal_arrays(result, expected, ("mask", "selected_names"))
        assert_close_arrays(
            result,
            expected,
            ("rank", "score"),
            rtol=1e-10,
            atol=1e-12,
            equal_nan=True,
        )
    elif name == "pca":
        assert_equal_arrays(result, expected, ("hvg_indices",))
        assert_close_arrays(
            result,
            expected,
            ("embedding", "loadings_hvg", "variance", "variance_ratio"),
            rtol=1e-6,
            atol=1e-8,
        )
    elif name == "neighbors":
        assert_equal_arrays(
            result,
            expected,
            (
                "distances_indices",
                "distances_indptr",
                "distances_shape",
                "connectivities_indices",
                "connectivities_indptr",
                "connectivities_shape",
            ),
        )
        assert_close_arrays(
            result,
            expected,
            ("distances_data", "connectivities_data"),
            rtol=1e-10,
            atol=1e-12,
        )
    elif name in {"spectral", "tsne", "umap"}:
        assert_close_arrays(
            result,
            expected,
            ("embedding",),
            rtol=1e-6,
            atol=1e-6,
        )
    elif name == "knnsc":
        assert_equal_arrays(result, expected, ("obs_names", "cluster_names"))
        assert_close_arrays(
            result,
            expected,
            ("shortest_path_lengths",),
            rtol=1e-6,
            atol=1e-8,
        )
    else:
        raise AssertionError(f"unexpected golden output: {name}")
