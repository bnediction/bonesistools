#!/usr/bin/env python

from pathlib import Path
from textwrap import dedent
from typing import Union
from zipfile import ZipFile

import pytest

import bonesistools as bt


def _write_text(path: Union[str, Path], content: str) -> Path:
    path = Path(path)
    path.write_text(dedent(content).strip() + "\n")
    return path


def _toy_ginml() -> str:
    return """
    <?xml version="1.0" encoding="UTF-8"?>
    <gxl>
      <graph class="regulatory" id="toy" nodeorder="A B C D E">
        <node id="A" maxvalue="1" input="true"/>
        <node id="B" maxvalue="1">
          <parameter val="1"/>
        </node>
        <node id="C" maxvalue="1">
          <value val="1">
            <exp str="A &amp; B"/>
          </value>
        </node>
        <node id="D" maxvalue="1">
          <value val="1">
            <exp str="!C"/>
          </value>
        </node>
        <node id="E" maxvalue="1">
          <value val="1">
            <exp str="A | D"/>
          </value>
        </node>
        <edge id="A:C" from="A" to="C" minvalue="1" sign="positive"/>
        <edge id="B:C" from="B" to="C" minvalue="1" sign="positive"/>
        <edge id="C:D" from="C" to="D" minvalue="1" sign="negative"/>
        <edge id="D:E" from="D" to="E" minvalue="1" sign="positive"/>
      </graph>
    </gxl>
    """


def test_read_ginml_returns_executable_model_with_boolean_network(tmp_path):
    path = _write_text(tmp_path / "toy.ginml", _toy_ginml())

    model = bt.bpy.io.read_ginml(path)

    assert isinstance(model, bt.bpy.io.ExecutableModel)
    assert model.metadata["graph"]["id"] == "toy"
    assert model.metadata["node_order"] == ["A", "B", "C", "D", "E"]

    graph = model.get("influence_graph")
    assert graph.number_of_nodes() == 5
    assert graph.number_of_edges() == 4
    assert graph.edge_sign("C", "D") == -1

    bn = model.get("boolean_network")
    assert bn.rules == {
        "A": "A",
        "B": "1",
        "C": "A & B",
        "D": "~C",
        "E": "A | D",
    }


def test_read_ginml_threshold_encodes_simple_multivalued_inputs(tmp_path):
    path = _write_text(
        tmp_path / "multi.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="multi" nodeorder="RA Target">
            <node id="RA" maxvalue="2" input="true"/>
            <node id="Target" maxvalue="1">
              <value val="1">
                <exp str="RA:2"/>
              </value>
            </node>
            <edge id="RA:Target" from="RA" to="Target"
              effects="1:positive 2:positive"/>
          </graph>
        </gxl>
        """,
    )

    model = bt.bpy.io.read_ginml(path)

    assert model.metadata["model_type"] == "multi-valued"
    assert model.metadata["multi_valued_components"] == {"RA": 2}
    assert model.get("boolean_network").rules == {
        "RA_b1": "RA_b1",
        "RA_b2": "RA_b1 & RA_b2",
        "Target": "RA_b1 & RA_b2",
    }
    assert model.get("influence_graph").edge_sign("RA", "Target") == 1


def test_read_ginml_keeps_unsupported_boolean_network_unavailable(tmp_path):
    path = _write_text(
        tmp_path / "unsupported.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="unsupported">
            <node id="RA" maxvalue="2" input="true"/>
            <node id="Target" maxvalue="1">
              <value val="1">
                <exp str="RA:3"/>
              </value>
            </node>
            <edge id="RA:Target" from="RA" to="Target" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    model = bt.bpy.io.read_ginml(path)

    assert model.boolean_network is None
    assert model.get("influence_graph").number_of_edges() == 1
    with pytest.raises(ValueError, match="Boolean network unavailable"):
        model.get("boolean_network")


def test_read_zginml_preserves_archive_metadata_and_companion_data(tmp_path):
    archive = tmp_path / "toy.zginml"

    with ZipFile(archive, "w") as zf:
        zf.writestr("GINsim-data/regulatoryGraph.ginml", _toy_ginml().strip())
        zf.writestr(
            "GINsim-data/initialState",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <initialStates>
                  <initialState name="state_A" value="A;1 B;0"/>
                  <initialState name="" value="A;0 B;1"/>
                  <initialState name="" value="A;1 B;1"/>
                  <input name="input_A" value="A;1"/>
                </initialStates>
                """).strip(),
        )
        zf.writestr(
            "GINsim-data/reg2dyn_parameters",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <simulationParameters nodeOrder="A B C D E">
                  <parameter name="run_A" updating="Asynchronous">
                    <initstates>
                      <row name="state_A"/>
                    </initstates>
                    <inputs>
                      <row name="input_A"/>
                    </inputs>
                  </parameter>
                </simulationParameters>
                """).strip(),
        )
        zf.writestr(
            "GINsim-data/modelSimplifier",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <modelModifierConfig>
                  <modelSimplifications>
                    <simplificationConfig name="Reduction"
                      strict="true" removeList="D E"/>
                  </modelSimplifications>
                </modelModifierConfig>
                """).strip(),
        )
        zf.writestr(
            "GINsim-data/mutant",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <perturbationConfig>
                  <listOfPerturbations>
                    <mutant name="A KO">
                      <change target="A" min="0" max="0"/>
                    </mutant>
                  </listOfPerturbations>
                  <listOfUsers>
                    <user key="simulation::A" value="A KO"/>
                  </listOfUsers>
                </perturbationConfig>
                """).strip(),
        )
        zf.writestr("GINsim-data/notes.txt", "preserve me")

    model = bt.bpy.io.read_zginml(archive)

    assert model.metadata["format"] == "zginml"
    assert (
        model.metadata["archive"]["main_ginml"] == "GINsim-data/regulatoryGraph.ginml"
    )
    assert "GINsim-data/notes.txt" in model.metadata["archive"]["files"]
    assert model.metadata["unknown_companion_files"] == ["GINsim-data/notes.txt"]

    assert model.initial_states["state_A"] == {"A": 1, "B": 0}
    assert model.initial_states["initialState_2"] == {"A": 0, "B": 1}
    assert model.initial_states["initialState_3"] == {"A": 1, "B": 1}
    assert model.initial_states["input_A"] == {"A": 1}
    assert model.metadata["initial_state_metadata"][1] == {
        "name": "initialState_2",
        "original_name": "",
        "section": "initialState",
    }
    assert model.metadata["simulation_parameters"]["parameters"]["run_A"] == {
        "attributes": {"name": "run_A", "updating": "Asynchronous"},
        "initstates": ["state_A"],
        "inputs": ["input_A"],
    }
    assert model.metadata["model_simplifier"]["simplifications"] == [
        {"name": "Reduction", "strict": "true", "removeList": ["D", "E"]}
    ]
    assert model.perturbations["A KO"] == [
        {
            "type": "change",
            "attributes": {"target": "A", "min": "0", "max": "0"},
        }
    ]
    assert model.metadata["perturbation_users"] == [
        {"key": "simulation::A", "value": "A KO"}
    ]


def test_executable_model_get_returns_available_containers_and_validates_attribute():
    model = bt.bpy.io.ExecutableModel()

    assert model.get("initial_states") == {}
    assert model.get("perturbations") == {}
    assert model.get("metadata") == {}

    with pytest.raises(ValueError, match="unsupported executable model attribute"):
        model.get("unknown")


def test_read_zginml_requires_main_ginml_file(tmp_path):
    archive = tmp_path / "empty.zginml"

    with ZipFile(archive, "w") as zf:
        zf.writestr("GINsim-data/notes.txt", "no model")

    with pytest.raises(ValueError, match="no '.ginml' file"):
        bt.bpy.io.read_zginml(archive)


def test_read_zginml_threshold_encodes_multivalued_initial_states(tmp_path):
    archive = tmp_path / "multi.zginml"

    with ZipFile(archive, "w") as zf:
        zf.writestr(
            "model.ginml",
            _multivalued_ginml().strip(),
        )
        zf.writestr(
            "initialState",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <initialStates>
                  <input name="RA_pharmacological" value="RA;2"/>
                  <input name="RA_physiological" value="RA;1"/>
                </initialStates>
                """).strip(),
        )

    model = bt.bpy.io.read_zginml(archive)

    assert model.initial_states["RA_pharmacological"] == {
        "RA_b1": 1,
        "RA_b2": 1,
    }
    assert model.initial_states["RA_physiological"] == {
        "RA_b1": 1,
        "RA_b2": 0,
    }


def _multivalued_ginml() -> str:
    return dedent("""
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="multi" nodeorder="RA Target">
            <node id="RA" maxvalue="2" input="true"/>
            <node id="Target" maxvalue="1">
              <value val="1">
                <exp str="RA:2"/>
              </value>
            </node>
            <edge id="RA:Target" from="RA" to="Target"
              effects="1:positive 2:positive"/>
          </graph>
        </gxl>
        """)
