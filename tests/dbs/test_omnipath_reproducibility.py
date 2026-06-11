#!/usr/bin/env python

import hashlib
import json
import os
from pathlib import Path

import pytest

from bonesistools.databases.omnipath import dorothea

REFERENCE_DIR = Path(__file__).with_name("data")
MODERN_REFERENCE = REFERENCE_DIR / "dorothea_current_mouse_A.sha256"
LEGACY_REFERENCE = REFERENCE_DIR / "dorothea_legacy_mouse.sha256"


def _signature_hash(graph):
    signature = {
        "nodes": sorted(str(node) for node in graph.nodes()),
        "edges": sorted(
            (str(source), str(target), int(data["sign"]))
            for source, target, data in graph.edges(data=True)
        ),
    }
    payload = json.dumps(signature, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _expected_hash(path):
    return path.read_text(encoding="utf-8").strip().split()[0]


@pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires live OmniPath access",
)
def test_dorothea_current_mouse_level_a_signature_is_reproducible():
    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version="latest",
        flavor="modern",
        hcop_version="bundled",
    )

    assert _signature_hash(graph) == _expected_hash(MODERN_REFERENCE)


@pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires live OmniPath archive access",
)
def test_dorothea_legacy_mouse_level_a_signature_is_reproducible():
    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version="2025-08-13",
        flavor="legacy",
        hcop_version="bundled",
    )

    assert _signature_hash(graph) == _expected_hash(LEGACY_REFERENCE)
