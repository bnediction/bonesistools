#!/usr/bin/env python

import json
import os
from pathlib import Path

import pytest

from bonesistools.boolpy.influence_graph import InfluenceGraph
from bonesistools.databases.omnipath import dorothea

REFERENCE_DIR = Path(__file__).with_name("data")
MODERN_REFERENCE = REFERENCE_DIR / "dorothea_mouse_level_a_modern_decoupler.json"
LEGACY_REFERENCE = REFERENCE_DIR / "dorothea_mouse_legacy_decoupler.json"


def _edge_signatures(graph):
    return {
        (source, target, data["sign"])
        for source, target, data in graph.edges(data=True)
    }


def _reference_graph(reference):
    graph = InfluenceGraph()
    graph.add_nodes_from(reference["nodes"])
    for edge in reference["edges"]:
        graph.add_edge(edge["source"], edge["target"], sign=edge["sign"])
    return graph


@pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_OMNIPATH_INTEGRATION") != "1",
    reason="requires live OmniPath archive and decoupler HCOP translation",
)
def test_dorothea_mouse_modern_matches_decoupler2_reference():
    reference = json.loads(MODERN_REFERENCE.read_text())
    expected = _reference_graph(reference)

    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version="latest",
        flavor="modern",
    )

    assert isinstance(graph, InfluenceGraph)
    assert sorted(graph.nodes()) == sorted(expected.nodes())
    assert _edge_signatures(graph) == _edge_signatures(expected)
    assert graph.number_of_nodes() == expected.number_of_nodes()
    assert graph.number_of_edges() == expected.number_of_edges()


@pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_OMNIPATH_INTEGRATION") != "1",
    reason="requires live OmniPath archive",
)
def test_dorothea_mouse_legacy_matches_decoupler1_reference():
    reference = json.loads(LEGACY_REFERENCE.read_text())
    expected = _reference_graph(reference)

    graph = dorothea(
        organism="mouse",
        levels=["A"],
        version="2025-08-13",
        flavor="legacy",
    )

    assert isinstance(graph, InfluenceGraph)
    assert sorted(graph.nodes()) == sorted(expected.nodes())
    assert _edge_signatures(graph) == _edge_signatures(expected)
    assert graph.number_of_nodes() == expected.number_of_nodes()
    assert graph.number_of_edges() == expected.number_of_edges()
