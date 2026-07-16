#!/usr/bin/env python

from pathlib import Path
from textwrap import dedent

import pytest

import bonesistools as bt

_UNSUPPORTED_SBML_DOCUMENT = (
    "unsupported SBML document "
    "(expected an SBML Level 3 document using the Qual package)"
)


def test_read_sbml_rejects_documents_without_level_3_qual(tmp_path):
    paths = [
        _write_text(
            tmp_path / "level2.sbml",
            """
            <sbml xmlns="http://www.sbml.org/sbml/level2/version4"
              level="2" version="4">
              <model id="metabolic"/>
            </sbml>
            """,
        ),
        _write_text(
            tmp_path / "core.sbml",
            """
            <sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
              level="3" version="1">
              <model id="metabolic"/>
            </sbml>
            """,
        ),
    ]

    for path in paths:
        with pytest.raises(ValueError) as error:
            bt.logic.io.read_sbml(path)

        assert str(error.value) == _UNSUPPORTED_SBML_DOCUMENT


def test_read_sbml_converts_boolean_terms_and_default(tmp_path):
    path = _write_text(
        tmp_path / "boolean.sbml",
        _sbml_document(
            species="""
              <qual:qualitativeSpecies qual:id="A" qual:maxLevel="1"
                qual:constant="true" qual:initialLevel="1"/>
              <qual:qualitativeSpecies qual:id="B" qual:maxLevel="1"
                qual:constant="false"/>
              <qual:qualitativeSpecies qual:id="C" qual:maxLevel="1"
                qual:constant="false"/>
            """,
            transitions="""
              <qual:transition qual:id="tr_B">
                <qual:listOfInputs>
                  <qual:input qual:id="tr_B_in" qual:qualitativeSpecies="A"
                    qual:sign="positive" qual:transitionEffect="none"/>
                </qual:listOfInputs>
                <qual:listOfOutputs>
                  <qual:output qual:id="tr_B_out" qual:qualitativeSpecies="B"
                    qual:transitionEffect="assignmentLevel"/>
                </qual:listOfOutputs>
                <qual:listOfFunctionTerms>
                  <qual:defaultTerm qual:resultLevel="0"/>
                  <qual:functionTerm qual:resultLevel="1">
                    <math xmlns="http://www.w3.org/1998/Math/MathML">
                      <apply><eq/><ci>A</ci><cn type="integer">1</cn></apply>
                    </math>
                  </qual:functionTerm>
                </qual:listOfFunctionTerms>
              </qual:transition>
              <qual:transition qual:id="tr_C">
                <qual:listOfInputs>
                  <qual:input qual:id="tr_C_in" qual:qualitativeSpecies="A"
                    qual:sign="negative" qual:transitionEffect="none"/>
                </qual:listOfInputs>
                <qual:listOfOutputs>
                  <qual:output qual:id="tr_C_out" qual:qualitativeSpecies="C"
                    qual:transitionEffect="assignmentLevel"/>
                </qual:listOfOutputs>
                <qual:listOfFunctionTerms>
                  <qual:defaultTerm qual:resultLevel="1"/>
                  <qual:functionTerm qual:resultLevel="0">
                    <math xmlns="http://www.w3.org/1998/Math/MathML">
                      <apply><eq/><ci>A</ci><cn type="integer">1</cn></apply>
                    </math>
                  </qual:functionTerm>
                </qual:listOfFunctionTerms>
              </qual:transition>
            """,
        ),
    )

    model = bt.logic.io.read_sbml(file=path)

    assert model.boolean_network.rules == {
        "A": "A",
        "B": "A",
        "C": "~A",
    }
    assert model.initial_conditions() == {"default": {"A": 1}}
    _assert_influence_graph_matches_boolean_network(model)


def test_read_sbml_booleanizes_multivalued_species_and_initial_levels(tmp_path):
    path = _write_text(
        tmp_path / "multivalued.sbml",
        _sbml_document(
            species="""
              <qual:qualitativeSpecies qual:id="X-signal" qual:name="Signal"
                qual:maxLevel="3" qual:constant="true" qual:initialLevel="2"/>
              <qual:qualitativeSpecies qual:id="Y" qual:maxLevel="2"
                qual:constant="false"/>
            """,
            transitions="""
              <qual:transition qual:id="tr_Y">
                <qual:listOfInputs>
                  <qual:input qual:id="tr_Y_in" qual:qualitativeSpecies="X-signal"
                    qual:sign="dual" qual:transitionEffect="none"/>
                </qual:listOfInputs>
                <qual:listOfOutputs>
                  <qual:output qual:id="tr_Y_out" qual:qualitativeSpecies="Y"
                    qual:transitionEffect="assignmentLevel"/>
                </qual:listOfOutputs>
                <qual:listOfFunctionTerms>
                  <qual:defaultTerm qual:resultLevel="0"/>
                  <qual:functionTerm qual:resultLevel="1">
                    <math xmlns="http://www.w3.org/1998/Math/MathML">
                      <apply>
                        <leq/><ci>X-signal</ci><cn type="integer">1</cn>
                      </apply>
                    </math>
                  </qual:functionTerm>
                  <qual:functionTerm qual:resultLevel="2">
                    <math xmlns="http://www.w3.org/1998/Math/MathML">
                      <apply>
                        <geq/><ci>X-signal</ci><cn type="integer">2</cn>
                      </apply>
                    </math>
                  </qual:functionTerm>
                </qual:listOfFunctionTerms>
              </qual:transition>
            """,
            layout="""
              <layout:listOfLayouts
                xmlns:layout="http://www.sbml.org/sbml/level3/version1/layout/version1">
                <layout:layout layout:id="layout">
                  <layout:listOfSpeciesGlyphs>
                    <layout:speciesGlyph layout:id="glyph_X"
                      layout:species="X-signal">
                      <layout:boundingBox>
                        <layout:position layout:x="10" layout:y="20"/>
                        <layout:dimensions layout:width="80" layout:height="40"/>
                      </layout:boundingBox>
                    </layout:speciesGlyph>
                  </layout:listOfSpeciesGlyphs>
                </layout:layout>
              </layout:listOfLayouts>
            """,
        ),
    )

    model = bt.logic.io.read_sbml(path)
    network = model.boolean_network
    graph = model.influence_graph

    assert network.components == {
        "X_signal_b1",
        "X_signal_b2",
        "X_signal_b3",
        "Y_b1",
        "Y_b2",
    }
    assert model.initial_conditions()["default"] == {
        "X_signal_b1": 1,
        "X_signal_b2": 1,
        "X_signal_b3": 0,
    }
    assert model.metadata()["boolean_component_names"] == {"X-signal": "X_signal"}
    assert network.rules["X_signal_b3"] == ("X_signal_b1 & X_signal_b2 & X_signal_b3")
    assert graph.nodes["X_signal_b1"]["sbml_component"] == "X-signal"
    assert graph.nodes["X_signal_b1"]["sbml_threshold"] == 1
    assert graph.nodes["X_signal_b3"]["sbml_threshold"] == 3
    assert graph.nodes["X_signal_b1"]["sbml_x"] == "10"
    assert graph.nodes["X_signal_b2"]["sbml_width"] == "80"

    low = network.next_configuration(
        {
            "X_signal_b1": 1,
            "X_signal_b2": 0,
            "X_signal_b3": 0,
            "Y_b1": 0,
            "Y_b2": 0,
        }
    )
    high = network.next_configuration(
        {
            "X_signal_b1": 1,
            "X_signal_b2": 1,
            "X_signal_b3": 0,
            "Y_b1": 1,
            "Y_b2": 0,
        }
    )
    assert (low["Y_b1"], low["Y_b2"]) == (1, 0)
    assert (high["Y_b1"], high["Y_b2"]) == (1, 1)
    _assert_influence_graph_matches_boolean_network(model)


def _assert_influence_graph_matches_boolean_network(model):
    imported = model.influence_graph
    expected = model.boolean_network.to_influence_graph()

    assert set(imported) == set(expected)
    assert {
        (source, target, attributes["sign"])
        for source, target, attributes in imported.edges(data=True)
    } == {
        (source, target, attributes["sign"])
        for source, target, attributes in expected.edges(data=True)
    }


def _sbml_document(species: str, transitions: str, layout: str = "") -> str:
    return f"""
    <?xml version="1.0" encoding="UTF-8"?>
    <sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1"
      xmlns:layout="http://www.sbml.org/sbml/level3/version1/layout/version1"
      level="3" version="1" qual:required="true" layout:required="false">
      <model id="logical_model">
        <qual:listOfQualitativeSpecies>
          {species}
        </qual:listOfQualitativeSpecies>
        <qual:listOfTransitions>
          {transitions}
        </qual:listOfTransitions>
        {layout}
      </model>
    </sbml>
    """


def _write_text(path: Path, content: str) -> Path:
    path.write_text(dedent(content).strip() + "\n")
    return path
