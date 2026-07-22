#!/usr/bin/env python

import hashlib
import json
import os
from pathlib import Path

import pytest

from bonesistools.resources.omnipath import dorothea

REFERENCE_DIR = Path(__file__).parent
ARCHIVE_VERSION = "2025-08-13"
MODERN_REFERENCE = REFERENCE_DIR / "dorothea_modern_mouse_A_2025-08-13.sha256"
LEGACY_REFERENCE = REFERENCE_DIR / "dorothea_legacy_mouse.sha256"

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


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


def test_dorothea_modern_mouse_level_a_signature_is_reproducible():
    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version=ARCHIVE_VERSION,
        flavor="modern",
        hcop_version="bundled",
        compatibility=True,
    )

    assert _signature_hash(graph) == _expected_hash(MODERN_REFERENCE)


def test_dorothea_legacy_mouse_level_a_signature_is_reproducible():
    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version=ARCHIVE_VERSION,
        flavor="legacy",
        hcop_version="bundled",
        compatibility=True,
    )

    assert _signature_hash(graph) == _expected_hash(LEGACY_REFERENCE)
