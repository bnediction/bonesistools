#!/usr/bin/env python

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, Dict, Tuple, cast

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from bonesistools.omics.tools._embedding import (
    _canonicalize_umap_curve_parameters,
    _orient_umap_spectral_layout,
    _prepare_umap_graph_for_embedding,
    _reproducible_umap_embedding,
)

GOLDEN_DIR = Path(__file__).parent
EXPECTED_DIR = GOLDEN_DIR / "expected"
DIAGNOSTIC_DIR = EXPECTED_DIR / "umap_diagnostics"
EXPECTED_PATH = DIAGNOSTIC_DIR / "checkpoints.npz"

N_COMPONENTS = 2
N_EPOCHS = 500
CHECKPOINT_EPOCHS = (0, 1, 2, 5, 10, 25, 50, 100, 250, N_EPOCHS)


def run_umap_diagnostics() -> Dict[str, np.ndarray]:
    """Return ordered numerical checkpoints from the golden UMAP calculation."""

    data, graph = _load_frozen_inputs()
    umap_module = importlib.import_module("umap.umap_")
    find_ab_params = cast(Any, getattr(umap_module, "find_ab_params"))
    spectral_layout = cast(Any, getattr(umap_module, "spectral_layout"))
    noisy_scale_coords = cast(Any, getattr(umap_module, "noisy_scale_coords"))

    a, b = _canonicalize_umap_curve_parameters(
        *cast(Tuple[float, float], find_ab_params(1.0, 0.5))
    )
    prepared_graph = cast(
        coo_matrix,
        _prepare_umap_graph_for_embedding(graph, N_EPOCHS),
    )
    random_state = np.random.RandomState(0)
    raw_layout = np.asarray(
        spectral_layout(
            data,
            prepared_graph,
            N_COMPONENTS,
            random_state,
            metric="euclidean",
            metric_kwds={},
        )
    )
    canonical_layout = _orient_umap_spectral_layout(raw_layout)
    noisy_initialization = np.asarray(
        noisy_scale_coords(
            canonical_layout,
            random_state,
            max_coord=10,
            noise=0.0001,
        )
    )
    normalized_initialization = _normalize_initialization(noisy_initialization)
    optimizer_rng_state = _optimizer_rng_state(random_state, umap_module)

    embedding_function = _reproducible_umap_embedding(umap_module)
    final_embedding, auxiliary = cast(
        Tuple[np.ndarray, Dict[str, Any]],
        embedding_function(
            data=data,
            graph=prepared_graph,
            n_components=N_COMPONENTS,
            initial_alpha=1.0,
            a=a,
            b=b,
            gamma=1.0,
            negative_sample_rate=5,
            n_epochs=list(CHECKPOINT_EPOCHS),
            init=noisy_initialization,
            random_state=random_state,
            metric="euclidean",
            metric_kwds={},
            densmap=False,
            densmap_kwds={},
            output_dens=False,
            verbose=False,
        ),
    )
    embeddings = cast(Any, auxiliary)["embedding_list"]
    if len(embeddings) != len(CHECKPOINT_EPOCHS):
        raise RuntimeError("UMAP did not return every requested diagnostic checkpoint")

    checkpoints = {
        "curve_parameters": np.asarray([a, b], dtype=np.float64),
        "prepared_graph_data": np.asarray(prepared_graph.data),
        "prepared_graph_row": np.asarray(prepared_graph.row),
        "prepared_graph_col": np.asarray(prepared_graph.col),
        "raw_spectral_layout": raw_layout,
        "canonical_spectral_layout": canonical_layout,
        "noisy_initialization": noisy_initialization,
        "normalized_initialization": normalized_initialization,
        "optimizer_rng_state": optimizer_rng_state,
    }
    for epoch, embedding in zip(CHECKPOINT_EPOCHS, embeddings):
        checkpoints[f"epoch_{epoch:03d}"] = np.asarray(embedding)
    checkpoints["final_embedding"] = np.asarray(final_embedding)
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


def _load_frozen_inputs() -> Tuple[np.ndarray, coo_matrix]:

    with np.load(EXPECTED_DIR / "pca.npz") as pca:
        data = np.asarray(pca["embedding"])
    with np.load(EXPECTED_DIR / "neighbors.npz") as neighbors:
        graph = cast(
            coo_matrix,
            csr_matrix(
                (
                    neighbors["connectivities_data"],
                    neighbors["connectivities_indices"],
                    neighbors["connectivities_indptr"],
                ),
                shape=tuple(neighbors["connectivities_shape"]),
            ).tocoo(),
        )
    return data, graph


def _normalize_initialization(initialization: np.ndarray) -> np.ndarray:

    minimum = np.min(initialization, axis=0)
    maximum = np.max(initialization, axis=0)
    return np.asarray(
        10.0 * (initialization - minimum) / (maximum - minimum),
        dtype=np.float32,
        order="C",
    )


def _optimizer_rng_state(
    random_state: np.random.RandomState,
    umap_module: Any,
) -> np.ndarray:

    cloned_random_state = np.random.RandomState()
    cloned_random_state.set_state(random_state.get_state())
    return cloned_random_state.randint(
        getattr(umap_module, "INT32_MIN"),
        getattr(umap_module, "INT32_MAX"),
        3,
    ).astype(np.int64)
