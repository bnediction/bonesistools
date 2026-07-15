#!/usr/bin/env python

from hashlib import sha256
from pathlib import Path

import pytest

import bonesistools as bt

from ._boolean_workflow import (
    EXPECTED_DIR,
    SEA_URCHIN_DORSAL_VENTRAL_AXIS_PATH,
)

_EXPECTED_NETWORK_PATH = EXPECTED_DIR / "sea_urchin_dorsal_ventral_axis.bnet"
_EXPECTED_SOURCE_SHA256 = (
    "699c6747774278e1966cb442808b716bc8a3b4cee87979f1fac45ccac465ae06"
)
_EXPECTED_NETWORK_SHA256 = (
    "e9f00617af7c5d97b771584b628a24d8d04b4365a5423fb96949a66ddccc6f5f"
)
_EXPECTED_INFLUENCE_SHA256 = (
    "5a3279f2da1fcdaab36b6db3bd98055d3472d0454e81c086461be9b41ed10d47"
)
_EXPECTED_MULTIVALUED_COMPONENTS = {
    "Alk1_2_3_6": 2,
    "Alk4_5_7": 2,
    "Bmp2_4": 2,
    "Bmp2_4_In": 2,
    "Chordin": 2,
    "Chordin_In": 2,
    "Nodal": 2,
    "Nodal_In": 3,
    "Smad1_4_5_8": 2,
    "Smad2_3_4": 2,
}


@pytest.fixture(scope="module")
def sea_urchin_model():
    return bt.logic.io.read_sbml(SEA_URCHIN_DORSAL_VENTRAL_AXIS_PATH)


def test_golden_sea_urchin_reference_files_are_unchanged():
    assert _file_sha256(SEA_URCHIN_DORSAL_VENTRAL_AXIS_PATH) == (
        _EXPECTED_SOURCE_SHA256
    )
    assert _file_sha256(_EXPECTED_NETWORK_PATH) == _EXPECTED_NETWORK_SHA256


def test_golden_sea_urchin_metadata_is_complete(sea_urchin_model):
    metadata = sea_urchin_model.metadata
    species = metadata["qualitative_species"]

    assert metadata["format"] == "sbml"
    assert metadata["sbml_level"] == 3
    assert metadata["sbml_version"] == 1
    assert metadata["model"]["name"] == ("Floc&#39;hlay2020 - SeaUrchin_model_ginsim")
    assert len(species) == 31
    assert {
        component: attributes["max_level"]
        for component, attributes in species.items()
        if attributes["max_level"] > 1
    } == _EXPECTED_MULTIVALUED_COMPONENTS
    assert len(metadata["transitions"]) == 24
    assert len(metadata["layout"]) == 31
    assert sea_urchin_model.initial_states == {}


def test_golden_sea_urchin_rules_match_biolqm_reference(sea_urchin_model):
    network = sea_urchin_model.get("boolean_network")
    expected = bt.logic.io.read_bnet(_EXPECTED_NETWORK_PATH, check=False)
    expected.pop("targets")
    expected.validate()

    assert network.components == expected.components

    for component in sorted(expected.components):
        expected_rule = network.ba.parse(expected.rule(component))
        assert bt.logic.ba.equivalence(
            network[component],
            expected_rule,
            ba=network.ba,
        ), component


def test_golden_sea_urchin_influence_graph_is_stable(sea_urchin_model):
    network = sea_urchin_model.get("boolean_network")
    graph = sea_urchin_model.get("influence_graph")
    expected_graph = network.to_influence_graph()
    influences = _influences(graph)

    assert graph.number_of_nodes() == 42
    assert graph.number_of_edges() == 113
    assert set(graph) == network.components
    assert influences == _influences(expected_graph)
    assert _influence_sha256(influences) == _EXPECTED_INFLUENCE_SHA256
    assert _edge_signs(graph, "Nodal_In_b1", "Alk4_5_7_b1") == (1,)
    assert _edge_signs(graph, "Nodal_In_b3", "Alk4_5_7_b2") == (1,)
    assert _edge_signs(graph, "Bmp2_4_In_b1", "Alk1_2_3_6_b1") == (-1, 1)
    assert _edge_signs(graph, "Chordin_In_b1", "Alk1_2_3_6_b1") == (-1,)
    assert _edge_signs(graph, "Smad2_3_4_b1", "Smad1_4_5_8_b1") == (-1,)


def test_golden_sea_urchin_threshold_nodes_preserve_source_metadata(
    sea_urchin_model,
):
    graph = sea_urchin_model.get("influence_graph")

    for threshold in (1, 2, 3):
        attributes = graph.nodes[f"Nodal_In_b{threshold}"]

        assert attributes["sbml_component"] == "Nodal_In"
        assert attributes["sbml_max_level"] == 3
        assert attributes["sbml_threshold"] == threshold
        assert attributes["sbml_x"] == "419"
        assert attributes["sbml_y"] == "52"
        assert attributes["sbml_width"] == "85"
        assert attributes["sbml_height"] == "35"


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _influences(graph):
    return tuple(
        sorted(
            (str(source), str(target), int(attributes["sign"]))
            for source, target, attributes in graph.edges(data=True)
        )
    )


def _influence_sha256(influences) -> str:
    content = "\n".join(
        f"{source}\t{target}\t{sign}" for source, target, sign in influences
    )
    return sha256(content.encode()).hexdigest()


def _edge_signs(graph, source: str, target: str):
    edge_data = graph.get_edge_data(source, target)

    if edge_data is None:
        return ()

    return tuple(sorted(int(attributes["sign"]) for attributes in edge_data.values()))
