#!/usr/bin/env python

import numpy as np
import pytest
from sklearn.manifold import trustworthiness
from typing_extensions import Literal

from ._omics_workflow import (
    assert_array_keys,
    assert_close_arrays,
    assert_equal_arrays,
    available_expected_names,
    load_expected,
    run_omics_workflow,
)

_PORTABLE_TSNE_QUALITY_NEIGHBORS = (15, 30, 50)
_PORTABLE_TSNE_QUALITY_MARGIN = 0.01
_STRICT_HVG_LOESS_SCORE_ATOL = 2e-15


@pytest.fixture(scope="module")
def golden_outputs():

    return run_omics_workflow()


@pytest.mark.parametrize("name", available_expected_names())
def test_golden_omics_workflow_matches_expected_outputs(
    golden_outputs,
    name,
    golden_mode: Literal["strict", "portable"],
):
    result = golden_outputs[name]
    expected = load_expected(name)

    assert_array_keys(result, expected)

    if golden_mode == "strict":
        if name == "hvg_loess":
            assert_equal_arrays(
                result,
                expected,
                tuple(key for key in expected if key != "score"),
            )
            assert_close_arrays(
                result,
                expected,
                ("score",),
                rtol=0.0,
                atol=_STRICT_HVG_LOESS_SCORE_ATOL,
                equal_nan=True,
            )
        else:
            assert_equal_arrays(result, expected, expected)
        return

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
    elif name in {"spectral", "umap"}:
        assert_close_arrays(
            result,
            expected,
            ("embedding",),
            rtol=1e-6,
            atol=1e-6,
        )
    elif name == "tsne":
        _assert_portable_tsne_quality(
            result["embedding"],
            expected["embedding"],
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


def _assert_portable_tsne_quality(
    embedding: np.ndarray,
    expected_embedding: np.ndarray,
) -> None:

    assert embedding.shape == expected_embedding.shape
    assert embedding.dtype == expected_embedding.dtype
    assert np.all(np.isfinite(embedding))

    representation = load_expected("pca")["embedding"]
    for n_neighbors in _PORTABLE_TSNE_QUALITY_NEIGHBORS:
        expected_trustworthiness = trustworthiness(
            representation,
            expected_embedding,
            n_neighbors=n_neighbors,
        )
        observed_trustworthiness = trustworthiness(
            representation,
            embedding,
            n_neighbors=n_neighbors,
        )
        assert observed_trustworthiness >= (
            expected_trustworthiness - _PORTABLE_TSNE_QUALITY_MARGIN
        )

        expected_continuity = trustworthiness(
            expected_embedding,
            representation,
            n_neighbors=n_neighbors,
        )
        observed_continuity = trustworthiness(
            embedding,
            representation,
            n_neighbors=n_neighbors,
        )
        assert observed_continuity >= (
            expected_continuity - _PORTABLE_TSNE_QUALITY_MARGIN
        )
