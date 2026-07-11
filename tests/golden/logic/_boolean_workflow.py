#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Tuple, cast

import numpy as np

import bonesistools as bt

GOLDEN_DIR = Path(__file__).parent
DATA_DIR = GOLDEN_DIR / "data"
EXPECTED_DIR = GOLDEN_DIR / "expected"
DEATH_RECEPTOR_CELL_FATE_PATH = DATA_DIR / "death_receptor_cell_fate.zginml"
SEA_URCHIN_DORSAL_VENTRAL_AXIS_PATH = DATA_DIR / "sea_urchin_dorsal_ventral_axis.sbml"
SYNTHETIC_30_NODE_DYNAMICS_PATH = DATA_DIR / "synthetic_30_node_dynamics.zginml"


def run_boolean_workflow() -> Dict[str, Dict[str, Any]]:

    death_receptor_model = bt.logic.io.read_zginml(DEATH_RECEPTOR_CELL_FATE_PATH)
    synthetic_model = bt.logic.io.read_zginml(SYNTHETIC_30_NODE_DYNAMICS_PATH)

    return {
        "death_receptor_cell_fate_mp": _collect_reachable_attractors(
            death_receptor_model,
            _death_receptor_cell_fate_initial_states(death_receptor_model),
            updates=("synchronous", "most-permissive"),
        ),
        "synthetic_30_node_dynamics": _collect_reachable_attractors(
            synthetic_model,
            _synthetic_30_node_initial_states(synthetic_model),
            updates=(
                "synchronous",
                "asynchronous",
                "general",
                "most-permissive",
            ),
        ),
    }


def save_expected(
    outputs: Mapping[str, Mapping[str, Any]],
    output_dir: Path = EXPECTED_DIR,
) -> None:

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, values in outputs.items():
        np.savez_compressed(output_dir / f"{name}.npz", **values)


def load_expected(
    name: str,
    expected_dir: Path = EXPECTED_DIR,
) -> Dict[str, np.ndarray]:

    with np.load(expected_dir / f"{name}.npz") as expected:
        return {key: expected[key] for key in expected.files}


def available_expected_names(expected_dir: Path = EXPECTED_DIR) -> Tuple[str, ...]:

    return tuple(sorted(path.stem for path in expected_dir.glob("*.npz")))


def _collect_reachable_attractors(
    model: Any,
    states: Mapping[str, Mapping[str, int]],
    *,
    updates: Tuple[
        str,
        ...,
    ],
) -> Dict[str, Any]:

    network = model.get("boolean_network")
    components = tuple(network.keys())

    case_states = []
    case_updates = []
    case_offsets = [0]
    attractor_sizes = []
    config_offsets = [0]
    configs = []

    for state_name, state in states.items():
        for update in updates:
            case_states.append(state_name)
            case_updates.append(update)
            attractors = network.reachable_attractors(
                state,
                update=cast(Any, update),
            )

            for attractor in _canonical_attractors(attractors, components):
                attractor_sizes.append(len(attractor))
                configs.extend(attractor)
                config_offsets.append(len(configs))

            case_offsets.append(len(attractor_sizes))

    return {
        "components": np.asarray(components, dtype=str),
        "rules": np.asarray([network.rule(component) for component in components]),
        "state_names": np.asarray(tuple(states), dtype=str),
        "state_values": _state_values_matrix(states, components),
        "case_states": np.asarray(case_states, dtype=str),
        "case_updates": np.asarray(case_updates, dtype=str),
        "case_offsets": np.asarray(case_offsets, dtype=np.int64),
        "attractor_sizes": np.asarray(attractor_sizes, dtype=np.int64),
        "config_offsets": np.asarray(config_offsets, dtype=np.int64),
        "configs": np.asarray(configs, dtype=np.int8),
        "graph_shape": np.asarray(
            [
                model.get("influence_graph").number_of_nodes(),
                model.get("influence_graph").number_of_edges(),
            ],
            dtype=np.int64,
        ),
    }


def _synthetic_30_node_initial_states(
    model: Any,
) -> Dict[str, Dict[str, int]]:

    return {
        "toggle_pair": _plain_state(model.initial_states["toggle_pair"]),
    }


def _death_receptor_cell_fate_initial_states(
    model: Any,
) -> Dict[str, Dict[str, int]]:

    archive_initial = _plain_state(model.initial_states["FAS1"])

    return {
        "archive_initial": dict(archive_initial),
        "tnf_fadd": {
            **archive_initial,
            "TNF": 1,
            "FASL": 0,
            "FADD": 1,
        },
        "fas_fadd": {
            **archive_initial,
            "TNF": 0,
            "FASL": 1,
            "FADD": 1,
        },
        "tnf_fas_fadd": {
            **archive_initial,
            "TNF": 1,
            "FASL": 1,
            "FADD": 1,
        },
        "no_ligand_fadd": {
            **archive_initial,
            "TNF": 0,
            "FASL": 0,
            "FADD": 1,
        },
        "tnf_no_fadd": {
            **archive_initial,
            "TNF": 1,
            "FASL": 0,
            "FADD": 0,
        },
    }


def _canonical_attractors(
    attractors: Tuple[bt.logic.ba.ConfigurationSet, ...],
    components: Tuple[str, ...],
) -> Tuple[Tuple[Tuple[int, ...], ...], ...]:

    canonical = []
    for attractor in attractors:
        configs = tuple(
            sorted(
                tuple(int(configuration[component]) for component in components)
                for configuration in attractor.enumerate()
            )
        )
        canonical.append(configs)

    return tuple(sorted(canonical))


def _state_values_matrix(
    states: Mapping[str, Mapping[str, int]],
    components: Tuple[str, ...],
) -> np.ndarray:

    values = np.full((len(states), len(components)), -1, dtype=np.int8)

    for row, state in enumerate(states.values()):
        for column, component in enumerate(components):
            if component in state:
                values[row, column] = state[component]

    return values


def _plain_state(state: Mapping[str, Any]) -> Dict[str, int]:

    return {
        key: int(cast(Any, value).value) if hasattr(value, "value") else int(value)
        for key, value in state.items()
    }
