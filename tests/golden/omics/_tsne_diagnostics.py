#!/usr/bin/env python

from __future__ import annotations

import hashlib
import importlib
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Optional, Tuple, cast

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
from typing_extensions import Protocol


class _SparseProbabilityMatrix(Protocol):
    data: np.ndarray
    indices: np.ndarray
    indptr: np.ndarray

    def copy(self) -> "_SparseProbabilityMatrix": ...

    def __imul__(self, value: float) -> "_SparseProbabilityMatrix": ...

    def __itruediv__(self, value: float) -> "_SparseProbabilityMatrix": ...


_BinarySearchPerplexity = Callable[[np.ndarray, float, int], np.ndarray]
_GradientDescent = Callable[..., Tuple[np.ndarray, float, int]]
_JointProbabilitiesNN = Callable[
    [csr_matrix, float, int],
    _SparseProbabilityMatrix,
]
_KLDivergenceBH = Callable[..., Tuple[float, np.ndarray]]

_utils_module = importlib.import_module("sklearn.manifold._utils")
_tsne_module = importlib.import_module("sklearn.manifold._t_sne")
_binary_search_perplexity = cast(
    _BinarySearchPerplexity,
    getattr(_utils_module, "_binary_search_perplexity"),
)
_gradient_descent = cast(
    _GradientDescent,
    getattr(_tsne_module, "_gradient_descent"),
)
_joint_probabilities_nn = cast(
    _JointProbabilitiesNN,
    getattr(_tsne_module, "_joint_probabilities_nn"),
)
_kl_divergence_bh = cast(
    _KLDivergenceBH,
    getattr(_tsne_module, "_kl_divergence_bh"),
)

GOLDEN_DIR = Path(__file__).parent
EXPECTED_DIR = GOLDEN_DIR / "expected"
DIAGNOSTIC_DIR = EXPECTED_DIR / "tsne_diagnostics"
EXPECTED_PATH = DIAGNOSTIC_DIR / "checkpoints.npz"

N_COMPONENTS = 2
N_ITER = 300
PERPLEXITY = 30.0
EARLY_EXAGGERATION = 12.0
LEARNING_RATE = 1000.0
EXPLORATION_N_ITER = 250
CHECKPOINT_ITERATIONS = (1, 2, 5, 10, 25, 50, 100, 250, N_ITER)


def run_tsne_diagnostics() -> Dict[str, np.ndarray]:
    """Return ordered numerical checkpoints from the golden t-SNE calculation."""

    data = _load_frozen_input()
    squared_distances = _squared_neighbor_distances(data)
    distance_rows = squared_distances.data.reshape(data.shape[0], -1).astype(
        np.float32,
        copy=False,
    )
    conditional_probabilities = _binary_search_perplexity(
        distance_rows,
        PERPLEXITY,
        0,
    )
    joint_probabilities = _joint_probabilities_nn(
        squared_distances.copy(),
        PERPLEXITY,
        0,
    )

    random_state = np.random.RandomState(0)
    initialization = 1e-4 * random_state.standard_normal(
        size=(data.shape[0], N_COMPONENTS)
    ).astype(np.float32)
    optimizer = _optimizer_checkpoints(joint_probabilities, initialization)

    checkpoints = {
        "knn_structure_digest": _digest_arrays(
            squared_distances.indices,
            squared_distances.indptr,
        ),
        "squared_knn_distances": squared_distances.data.copy(),
        "conditional_probabilities_sample": _sample(
            conditional_probabilities
        ),
        "conditional_probabilities_digest": _digest_arrays(
            conditional_probabilities
        ),
        "joint_probability_structure_digest": _digest_arrays(
            joint_probabilities.indices,
            joint_probabilities.indptr,
        ),
        "joint_probabilities_sample": _sample(joint_probabilities.data),
        "joint_probabilities_digest": _digest_arrays(
            joint_probabilities.data
        ),
        "random_initialization": initialization,
    }
    checkpoints.update(optimizer)
    return checkpoints


def save_expected(
    checkpoints: Dict[str, np.ndarray],
    output_path: Path = EXPECTED_PATH,
) -> None:

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cast(Any, np.savez_compressed)(output_path, **checkpoints)


def load_expected(output_path: Path = EXPECTED_PATH) -> Dict[str, np.ndarray]:

    with np.load(output_path) as expected:
        return {key: expected[key] for key in expected.files}


def _load_frozen_input() -> np.ndarray:

    with np.load(EXPECTED_DIR / "pca.npz") as pca:
        return np.asarray(pca["embedding"])


def _squared_neighbor_distances(data: np.ndarray) -> csr_matrix:

    n_neighbors = min(data.shape[0] - 1, int(3.0 * PERPLEXITY + 1))
    nearest_neighbors = NearestNeighbors(
        algorithm="auto",
        n_jobs=1,
        n_neighbors=n_neighbors,
        metric="euclidean",
    )
    nearest_neighbors.fit(data)
    distances = cast(
        csr_matrix,
        nearest_neighbors.kneighbors_graph(mode="distance"),
    )
    distances.data **= 2
    distances.sort_indices()
    return distances


def _optimizer_checkpoints(
    probabilities: Any,
    initialization: np.ndarray,
) -> Dict[str, np.ndarray]:

    embeddings: Dict[int, np.ndarray] = {}
    initial_gradient: Optional[np.ndarray] = None
    iteration = 0

    def objective(
        parameters: np.ndarray,
        *args: Any,
        **kwargs: Any,
    ) -> Tuple[float, np.ndarray]:
        nonlocal initial_gradient, iteration

        if iteration in CHECKPOINT_ITERATIONS:
            embeddings[iteration] = parameters.reshape(
                initialization.shape
            ).copy()
        error, gradient = _kl_divergence_bh(
            parameters,
            *args,
            **kwargs,
        )
        if iteration == 0:
            initial_gradient = gradient.copy()
        iteration += 1
        return error, gradient

    optimized_probabilities = probabilities.copy()
    optimized_probabilities *= EARLY_EXAGGERATION
    parameters, error, last_iteration = _gradient_descent(
        objective,
        initialization.ravel(),
        it=0,
        max_iter=EXPLORATION_N_ITER,
        n_iter_check=50,
        n_iter_without_progress=EXPLORATION_N_ITER,
        momentum=0.5,
        learning_rate=LEARNING_RATE,
        min_gain=0.01,
        min_grad_norm=1e-7,
        args=[
            optimized_probabilities,
            max(N_COMPONENTS - 1, 1),
            initialization.shape[0],
            N_COMPONENTS,
        ],
        kwargs={
            "angle": 0.5,
            "skip_num_points": 0,
            "verbose": False,
            "num_threads": 1,
        },
    )

    optimized_probabilities /= EARLY_EXAGGERATION
    parameters, error, last_iteration = _gradient_descent(
        objective,
        parameters,
        it=last_iteration + 1,
        max_iter=N_ITER,
        n_iter_check=50,
        n_iter_without_progress=300,
        momentum=0.8,
        learning_rate=LEARNING_RATE,
        min_gain=0.01,
        min_grad_norm=1e-7,
        args=[
            optimized_probabilities,
            max(N_COMPONENTS - 1, 1),
            initialization.shape[0],
            N_COMPONENTS,
        ],
        kwargs={
            "angle": 0.5,
            "skip_num_points": 0,
            "verbose": False,
            "num_threads": 1,
        },
    )
    embeddings[N_ITER] = parameters.reshape(initialization.shape).copy()

    if initial_gradient is None or last_iteration + 1 != N_ITER:
        raise RuntimeError("t-SNE did not reach every diagnostic checkpoint")
    missing = set(CHECKPOINT_ITERATIONS) - set(embeddings)
    if missing:
        raise RuntimeError(
            f"t-SNE did not return diagnostic checkpoints: {sorted(missing)!r}"
        )

    checkpoints = {
        "initial_gradient": initial_gradient,
    }
    for checkpoint in CHECKPOINT_ITERATIONS:
        checkpoints[f"iteration_{checkpoint:03d}"] = embeddings[checkpoint]
    checkpoints["final_kl_divergence"] = np.asarray([error], dtype=np.float64)
    return checkpoints


def _sample(values: np.ndarray, size: int = 256) -> np.ndarray:

    flattened = np.asarray(values).ravel()
    if flattened.size <= size:
        return flattened.copy()
    indices = np.linspace(0, flattened.size - 1, num=size, dtype=np.int64)
    return flattened[indices]


def _digest_arrays(*arrays: Iterable[Any]) -> np.ndarray:

    digest = hashlib.sha256()
    for values in arrays:
        array = np.ascontiguousarray(values)
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return np.frombuffer(digest.digest(), dtype=np.uint8).copy()
