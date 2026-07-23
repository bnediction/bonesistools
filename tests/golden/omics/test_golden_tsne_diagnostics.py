#!/usr/bin/env python

import platform

import numpy as np
import pytest
import sklearn
from sklearn.manifold import TSNE
from threadpoolctl import threadpool_limits

from ._tsne_diagnostics import (
    EARLY_EXAGGERATION,
    LEARNING_RATE,
    N_COMPONENTS,
    N_ITER,
    PERPLEXITY,
    _load_frozen_input,
    load_expected,
    run_tsne_diagnostics,
)


@pytest.fixture(scope="module")
def tsne_diagnostics():

    return run_tsne_diagnostics()


def test_golden_tsne_numerical_checkpoints(tsne_diagnostics):
    expected = load_expected()
    assert tuple(tsne_diagnostics) == tuple(expected)

    context = (
        f"platform={platform.platform()}, "
        f"machine={platform.machine()}, "
        f"scikit-learn={sklearn.__version__}"
    )
    for checkpoint, expected_values in expected.items():
        if checkpoint == "squared_knn_distances":
            np.testing.assert_allclose(
                tsne_diagnostics[checkpoint],
                expected_values,
                rtol=2e-14,
                atol=1e-12,
                err_msg=(
                    f"first divergent t-SNE checkpoint: {checkpoint!r} "
                    f"({context})"
                ),
            )
            continue
        np.testing.assert_array_equal(
            tsne_diagnostics[checkpoint],
            expected_values,
            err_msg=(
                f"first divergent t-SNE checkpoint: {checkpoint!r} "
                f"({context})"
            ),
        )


def test_golden_tsne_diagnostic_matches_public_solver(tsne_diagnostics):
    with threadpool_limits(limits=1):
        embedding = TSNE(
            n_components=N_COMPONENTS,
            perplexity=PERPLEXITY,
            early_exaggeration=EARLY_EXAGGERATION,
            learning_rate=LEARNING_RATE,
            max_iter=N_ITER,
            metric="euclidean",
            random_state=np.random.RandomState(0),
            n_jobs=1,
            method="barnes_hut",
            init="random",
        ).fit_transform(_load_frozen_input())

    np.testing.assert_array_equal(
        tsne_diagnostics[f"iteration_{N_ITER:03d}"],
        embedding,
    )
