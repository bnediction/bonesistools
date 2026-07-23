#!/usr/bin/env python

import importlib
import platform
from pathlib import Path
from typing import Dict, cast

import numpy as np
import pytest
import sklearn
from sklearn.manifold import TSNE
from threadpoolctl import threadpool_limits
from typing_extensions import Protocol


class _TSNEDiagnosticsModule(Protocol):
    EARLY_EXAGGERATION: float
    LEARNING_RATE: float
    N_COMPONENTS: int
    N_ITER: int
    PERPLEXITY: float

    def _load_frozen_input(self) -> np.ndarray: ...

    def load_expected(
        self,
        output_path: Path = ...,
    ) -> Dict[str, np.ndarray]: ...

    def run_tsne_diagnostics(self) -> Dict[str, np.ndarray]: ...


_tsne_diagnostics = cast(
    _TSNEDiagnosticsModule,
    importlib.import_module("tests.golden.omics._tsne_diagnostics"),
)

_NUMERICAL_TOLERANCES = {
    "squared_knn_distances": (2e-14, 1e-12),
    "conditional_probabilities": (2e-15, 1e-18),
    "joint_probabilities": (2e-14, 1e-18),
}


@pytest.fixture(scope="module")
def tsne_diagnostics():

    return _tsne_diagnostics.run_tsne_diagnostics()


def test_golden_tsne_numerical_checkpoints(tsne_diagnostics):
    expected = _tsne_diagnostics.load_expected()
    assert tuple(tsne_diagnostics) == tuple(expected)

    context = (
        f"platform={platform.platform()}, "
        f"machine={platform.machine()}, "
        f"scikit-learn={sklearn.__version__}"
    )
    for checkpoint, expected_values in expected.items():
        if checkpoint in _NUMERICAL_TOLERANCES:
            rtol, atol = _NUMERICAL_TOLERANCES[checkpoint]
            np.testing.assert_allclose(
                tsne_diagnostics[checkpoint],
                expected_values,
                rtol=rtol,
                atol=atol,
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
            n_components=_tsne_diagnostics.N_COMPONENTS,
            perplexity=_tsne_diagnostics.PERPLEXITY,
            early_exaggeration=_tsne_diagnostics.EARLY_EXAGGERATION,
            learning_rate=_tsne_diagnostics.LEARNING_RATE,
            max_iter=_tsne_diagnostics.N_ITER,
            metric="euclidean",
            random_state=np.random.RandomState(0),
            n_jobs=1,
            method="barnes_hut",
            init="random",
        ).fit_transform(_tsne_diagnostics._load_frozen_input())

    np.testing.assert_array_equal(
        tsne_diagnostics[
            f"iteration_{_tsne_diagnostics.N_ITER:03d}"
        ],
        embedding,
    )
