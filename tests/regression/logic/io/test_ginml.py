#!/usr/bin/env python

from itertools import product
from pathlib import Path
from textwrap import dedent
from typing import Union
from zipfile import ZipFile

import pytest

import bonesistools as bt
from bonesistools.logic.input_output._executable_model import ExecutableModel


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
        <edge id="A:E" from="A" to="E" minvalue="1" sign="positive"/>
        <edge id="D:E" from="D" to="E" minvalue="1" sign="positive"/>
      </graph>
    </gxl>
    """


def test_read_ginml_returns_executable_model_with_boolean_network(tmp_path):
    path = _write_text(tmp_path / "toy.ginml", _toy_ginml())

    model = bt.logic.io.read_ginml(file=path)

    assert isinstance(model, ExecutableModel)
    assert repr(model).startswith(
        "ExecutableModel(boolean_network(components=5), "
        "influence_graph(nodes=5, edges=6), "
    )
    assert model.metadata["graph"]["id"] == "toy"
    assert model.metadata["node_order"] == ["A", "B", "C", "D", "E"]

    graph = model.get("influence_graph")
    assert graph.number_of_nodes() == 5
    assert graph.number_of_edges() == 6
    assert graph.edge_sign("A", "A") == 1
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
            <node id="RA" name="Retinoic acid" maxvalue="2" input="true"
              custom="shared"/>
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

    model = bt.logic.io.read_ginml(path)

    assert model.metadata["model_type"] == "multi-valued"
    assert model.metadata["multi_valued_components"] == {"RA": 2}
    assert model.get("boolean_network").rules == {
        "RA_b1": "RA_b1",
        "RA_b2": "RA_b1 & RA_b2",
        "Target": "RA_b1 & RA_b2",
    }
    graph = model.get("influence_graph")
    assert set(graph) == {"RA_b1", "RA_b2", "Target"}
    assert graph.nodes["RA_b1"]["ginml_threshold"] == 1
    assert graph.nodes["RA_b2"]["ginml_threshold"] == 2
    assert graph.nodes["RA_b1"]["ginml_component"] == "RA"
    assert graph.nodes["RA_b2"]["ginml_component"] == "RA"
    assert graph.nodes["RA_b1"]["ginml_name"] == "Retinoic acid"
    assert graph.nodes["RA_b2"]["ginml_name"] == "Retinoic acid"
    assert graph.nodes["RA_b1"]["ginml_custom"] == "shared"
    assert graph.nodes["RA_b2"]["ginml_custom"] == "shared"
    assert graph.edge_sign("RA_b1", "Target") == 1
    assert graph.edge_sign("RA_b2", "Target") == 1
    assert graph.edge_sign("RA_b1", "RA_b2") == 1
    assert graph["RA_b1"]["Target"][0]["id"] == "RA:Target"
    assert graph["RA_b2"]["Target"][0]["id"] == "RA:Target"
    expected_graph = model.get("boolean_network").to_influence_graph()
    assert set(graph) == set(expected_graph)
    assert {
        (source, target, attributes["sign"])
        for source, target, attributes in graph.edges(data=True)
    } == {
        (source, target, attributes["sign"])
        for source, target, attributes in expected_graph.edges(data=True)
    }


def test_read_ginml_regularizes_multivalued_target_thresholds(tmp_path):
    path = _write_text(
        tmp_path / "regularized.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="regularized" nodeorder="A B X">
            <node id="A" maxvalue="1" input="true"/>
            <node id="B" maxvalue="1" input="true"/>
            <node id="X" maxvalue="2">
              <value val="1"><exp str="A"/></value>
              <value val="2"><exp str="B"/></value>
            </node>
            <edge id="A:X" from="A" to="X" minvalue="1" sign="positive"/>
            <edge id="B:X" from="B" to="X" minvalue="1" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    network = bt.logic.io.read_ginml(path).get("boolean_network")

    assert network.components == {"A", "B", "X_b1", "X_b2"}

    for a, b, x_b1, x_b2 in product((0, 1), repeat=4):
        configuration = {"A": a, "B": b, "X_b1": x_b1, "X_b2": x_b2}
        observed = network.next_configuration(configuration)
        expected = {
            "A": a,
            "B": b,
            "X_b1": int(a or b or (x_b1 and x_b2)),
            "X_b2": int(x_b1 and b),
        }

        assert observed == expected
        assert observed["X_b2"] <= observed["X_b1"]


def test_read_ginml_expands_multivalued_intervals_by_exact_level(tmp_path):
    path = _write_text(
        tmp_path / "interval.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="interval" nodeorder="X Target">
            <node id="X" maxvalue="3" input="true"/>
            <node id="Target" maxvalue="1">
              <value val="1"><exp str="X:1"/></value>
            </node>
            <edge id="X:Target:1" from="X" to="Target"
              minvalue="1" maxvalue="2" sign="positive"/>
            <edge id="X:Target:2" from="X" to="Target"
              minvalue="3" sign="negative"/>
          </graph>
        </gxl>
        """,
    )

    network = bt.logic.io.read_ginml(path).get("boolean_network")

    assert network.rules["Target"] == ("(X_b1 & ~X_b2) | (X_b1 & X_b2 & ~X_b3)")
    assert (
        network.next_state("Target", {"X_b1": 1, "X_b2": 0, "X_b3": 1, "Target": 0})
        == 1
    )
    assert (
        network.next_state("Target", {"X_b1": 1, "X_b2": 1, "X_b3": 1, "Target": 0})
        == 0
    )


def test_read_ginml_negates_complete_multivalued_regulator_condition(tmp_path):
    path = _write_text(
        tmp_path / "negated-multivalued.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="negated" nodeorder="X Target">
            <node id="X" maxvalue="2" input="true"/>
            <node id="Target" maxvalue="1">
              <value val="1"><exp str="!X"/></value>
            </node>
            <edge id="X:Target" from="X" to="Target"
              minvalue="1" maxvalue="2" sign="negative"/>
          </graph>
        </gxl>
        """,
    )

    network = bt.logic.io.read_ginml(path).get("boolean_network")

    assert network.rules["Target"] == "~((X_b1 & ~X_b2) | (X_b1 & X_b2))"
    assert network.next_state("Target", {"X_b1": 0, "X_b2": 0, "Target": 0}) == 1
    assert network.next_state("Target", {"X_b1": 1, "X_b2": 0, "Target": 0}) == 0
    assert network.next_state("Target", {"X_b1": 1, "X_b2": 1, "Target": 0}) == 0


def test_read_ginml_explicit_context_overrides_multivalued_function(tmp_path):
    path = _write_text(
        tmp_path / "override.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="override" nodeorder="S V X">
            <node id="S" maxvalue="2" input="true"/>
            <node id="V" maxvalue="1" input="true"/>
            <node id="X" maxvalue="2">
              <parameter idActiveInteractions="S:X:1 V:X" val="2"/>
              <value val="2"><exp str="S:2"/></value>
              <value val="1">
                <exp str="(V &amp; !S) | (S:1 &amp; !V)"/>
              </value>
            </node>
            <edge id="S:X" from="S" to="X"
              effects="1:positive 2:positive"/>
            <edge id="V:X" from="V" to="X"
              minvalue="1" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    network = bt.logic.io.read_ginml(path).get("boolean_network")
    observed = network.next_configuration(
        {"S_b1": 1, "S_b2": 0, "V": 1, "X_b1": 1, "X_b2": 0}
    )

    assert observed["X_b1"] == 1
    assert observed["X_b2"] == 1


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

    with pytest.warns(UserWarning, match="Boolean-network conversion failed"):
        model = bt.logic.io.read_ginml(path)

    assert model.boolean_network is None
    assert model.get("influence_graph").number_of_edges() == 1
    with pytest.raises(ValueError, match="Boolean network unavailable"):
        model.get("boolean_network")


def test_read_ginml_converts_boolean_regulatory_contexts(tmp_path):
    path = _write_text(
        tmp_path / "contexts.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="contexts">
            <node id="A" maxvalue="1" input="true"/>
            <node id="B" maxvalue="1">
              <parameter idActiveInteractions="A:B" val="1"/>
            </node>
            <node id="C" maxvalue="1">
              <parameter val="1"/>
            </node>
            <node id="D" maxvalue="1" basevalue="1"/>
            <node id="G" maxvalue="1" basevalue="1"/>
            <node id="E" maxvalue="1">
              <parameter idActiveInteractions="A:E" val="1"/>
              <value val="1"><exp str="A"/></value>
            </node>
            <node id="X-Y" maxvalue="1" input="true"/>
            <node id="F" maxvalue="1">
              <value val="1"><exp str="X-Y &amp; A"/></value>
            </node>
            <edge id="A:B" from="A" to="B" minvalue="1" sign="positive"/>
            <edge id="B:C" from="B" to="C" minvalue="1" sign="negative"/>
            <edge id="A:E" from="A" to="E" minvalue="1" sign="positive"/>
            <edge id="A:G" from="A" to="G" minvalue="1" sign="negative"/>
            <edge id="X-Y:F" from="X-Y" to="F" minvalue="1" sign="positive"/>
            <edge id="A:F" from="A" to="F" minvalue="1" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    model = bt.logic.io.read_ginml(path)

    assert model.get("boolean_network").rules == {
        "A": "A",
        "B": "A",
        "C": "~B",
        "D": "1",
        "G": "~A",
        "E": "A",
        "X_Y": "X_Y",
        "F": "X_Y & A",
    }
    assert model.metadata["boolean_component_names"] == {"X-Y": "X_Y"}


def test_read_ginml_rejects_colliding_normalized_component_names(tmp_path):
    path = _write_text(
        tmp_path / "collision.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="collision">
            <node id="A-B" maxvalue="1" input="true"/>
            <node id="A_B" maxvalue="1" input="true"/>
            <edge id="A-B:A_B" from="A-B" to="A_B"
              minvalue="1" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    with pytest.warns(
        UserWarning,
        match="Boolean-network conversion failed.*duplicate names: A_B",
    ):
        model = bt.logic.io.read_ginml(path)

    assert model.boolean_network is None
    assert "duplicate names: A_B" in model.metadata["boolean_network_reason"]


def test_read_ginml_resolves_named_node_styles(tmp_path):
    path = _write_text(
        tmp_path / "styles.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="styles">
            <nodestyle background="#ffffff" foreground="#000000"
              text="#000000" shape="RECTANGLE" width="120" height="45"/>
            <nodestyle name="Arrest" background="#ff9999"/>
            <node id="CDKN1A" name="Cyclin-dependent kinase inhibitor 1A"
              maxvalue="1">
              <parameter val="1"/>
              <nodevisualsetting x="10" y="20" style="Arrest"/>
            </node>
            <node id="RB1" maxvalue="1">
              <parameter val="1"/>
              <nodevisualsetting x="30" y="40" style=""/>
            </node>
            <edge id="CDKN1A:RB1" from="CDKN1A" to="RB1" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    model = bt.logic.io.read_ginml(path)
    graph = model.get("influence_graph")

    assert model.metadata["node_styles"] == {
        "default": {
            "background": "#ffffff",
            "foreground": "#000000",
            "text": "#000000",
            "shape": "RECTANGLE",
            "width": "120",
            "height": "45",
        },
        "Arrest": {"background": "#ff9999"},
    }
    assert model.metadata["node_style_groups"] == {"Arrest": ["CDKN1A"]}
    assert model.metadata["node_visual_settings"]["CDKN1A"] == {
        "background": "#ff9999",
        "foreground": "#000000",
        "text": "#000000",
        "shape": "RECTANGLE",
        "width": "120",
        "height": "45",
        "x": "10",
        "y": "20",
        "style": "Arrest",
    }
    assert model.metadata["node_visual_settings"]["RB1"] == {
        "background": "#ffffff",
        "foreground": "#000000",
        "text": "#000000",
        "shape": "RECTANGLE",
        "width": "120",
        "height": "45",
        "x": "30",
        "y": "40",
        "style": "",
    }
    assert graph.nodes["CDKN1A"]["name"] == "Cyclin-dependent kinase inhibitor 1A"
    assert graph.nodes["CDKN1A"]["ginml_style"] == "Arrest"
    assert graph.nodes["CDKN1A"]["ginml_background"] == "#ff9999"
    assert graph.nodes["CDKN1A"]["ginml_x"] == "10"
    assert "ginml_style" not in graph.nodes["RB1"]
    assert graph.nodes["RB1"]["ginml_background"] == "#ffffff"


def test_read_ginml_preserves_direct_node_visual_settings(tmp_path):
    path = _write_text(
        tmp_path / "direct-visual.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl>
          <graph class="regulatory" id="direct">
            <node id="MEK" maxvalue="1">
              <parameter val="1"/>
              <nodevisualsetting>
                <ellipse x="746" y="459" width="80" height="30"
                  backgroundColor="#ffcc66" foregroundColor="#000000"/>
              </nodevisualsetting>
            </node>
            <node id="ERK" maxvalue="1">
              <parameter val="1"/>
              <nodevisualsetting>
                <rect x="953" y="421" width="80" height="30"
                  backgroundColor="#cccccc" foregroundColor="#111111"/>
              </nodevisualsetting>
            </node>
            <edge id="MEK:ERK" from="MEK" to="ERK" sign="positive"/>
          </graph>
        </gxl>
        """,
    )

    model = bt.logic.io.read_ginml(path)
    graph = model.get("influence_graph")

    assert model.metadata["node_styles"] == {}
    assert model.metadata["node_style_groups"] == {}
    assert model.metadata["node_visual_settings"]["MEK"] == {
        "shape": "ellipse",
        "x": "746",
        "y": "459",
        "width": "80",
        "height": "30",
        "background": "#ffcc66",
        "foreground": "#000000",
    }
    assert graph.nodes["MEK"]["ginml_shape"] == "ellipse"
    assert graph.nodes["MEK"]["ginml_background"] == "#ffcc66"
    assert graph.nodes["MEK"]["ginml_foreground"] == "#000000"
    assert graph.nodes["MEK"]["ginml_width"] == "80"
    assert graph.nodes["ERK"]["ginml_shape"] == "rect"
    assert graph.nodes["ERK"]["ginml_background"] == "#cccccc"


def test_read_ginml_preserves_graph_edge_styles_and_annotations(tmp_path):
    path = _write_text(
        tmp_path / "metadata.ginml",
        """
        <?xml version="1.0" encoding="UTF-8"?>
        <gxl xmlns:xlink="http://www.w3.org/1999/xlink">
          <graph class="regulatory" id="metadata" nodeorder="A B">
            <attr name="display.node" value="name"/>
            <annotation>
              <comment>Graph note</comment>
              <linklist>
                <link xlink:href="https://example.org/model"/>
              </linklist>
            </annotation>
            <nodestyle background="#ffffff" text="#000000"/>
            <edgestyle color="#000000" line_width="1"
              properties="positive:#00c800 negative:#c80000 dual:#0000c8"/>
            <edgestyle name="Tick" line_width="3"/>
            <node id="A" maxvalue="1" input="true" basevalue="1">
              <nodevisualsetting x="10" y="20" style=""/>
            </node>
            <node id="B" maxvalue="1">
              <value val="1"><exp str="A"/></value>
              <annotation><comment>Node note</comment></annotation>
              <nodevisualsetting x="30" y="40" style=""/>
            </node>
            <edge id="A:B" from="A" to="B" minvalue="1" sign="positive">
              <annotation><comment>Edge note</comment></annotation>
              <edgevisualsetting points="10,20 30,40" anchor="NE" style="Tick"/>
            </edge>
          </graph>
        </gxl>
        """,
    )

    model = bt.logic.io.read_ginml(path)
    graph = model.get("influence_graph")
    edge = next(iter(graph.get_edge_data("A", "B").values()))

    assert model.metadata["graph_attributes"] == {"display.node": "name"}
    assert model.metadata["annotations"][0]["children"][0] == {
        "tag": "comment",
        "attributes": {},
        "text": "Graph note",
    }
    assert model.metadata["annotations"][0]["children"][1]["children"][0] == {
        "tag": "link",
        "attributes": {"href": "https://example.org/model"},
    }
    assert (
        model.metadata["nodes"]["B"]["annotations"][0]["children"][0]["text"]
        == "Node note"
    )
    assert model.metadata["edge_styles"]["Tick"] == {"line_width": "3"}
    assert model.metadata["edge_style_groups"] == {"Tick": ["A:B"]}
    assert model.metadata["edge_visual_settings"]["A:B"] == {
        "color": "#00c800",
        "line_width": "3",
        "properties": "positive:#00c800 negative:#c80000 dual:#0000c8",
        "points": "10,20 30,40",
        "anchor": "NE",
        "style": "Tick",
    }
    assert (
        model.metadata["edge_annotations"]["A:B"][0]["children"][0]["text"]
        == "Edge note"
    )
    assert graph.graph["ginml_id"] == "metadata"
    assert graph.graph["ginml_display_node"] == "name"
    assert graph.nodes["A"]["ginml_basevalue"] == "1"
    assert edge["ginml_style"] == "Tick"
    assert edge["ginml_color"] == "#00c800"
    assert edge["ginml_line_width"] == "3"
    assert edge["ginml_points"] == "10,20 30,40"


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
                  <priorityClassList id="Priorities">
                    <class name="Fast" mode="1" rank="1" content="A B "/>
                  </priorityClassList>
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
                  <stripOutput key="E"/>
                </modelModifierConfig>
                """).strip(),
        )
        zf.writestr(
            "GINsim-data/avatar_parameters",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <avatarParameters nodeOrder="A B C D E">
                  <parameter name="AVATAR" avatarparameters="algorithm=0">
                    <stateList states="{A=[1]}" namestates="state_A"/>
                  </parameter>
                </avatarParameters>
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

    model = bt.logic.io.read_zginml(file=archive)

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
    assert model.metadata["simulation_parameters"]["priority_class_lists"] == [
        {
            "attributes": {"id": "Priorities"},
            "classes": [
                {
                    "name": "Fast",
                    "mode": "1",
                    "rank": "1",
                    "content": ["A", "B"],
                }
            ],
        }
    ]
    assert model.metadata["model_simplifier"]["simplifications"] == [
        {"name": "Reduction", "strict": "true", "removeList": ["D", "E"]}
    ]
    assert model.metadata["model_simplifier"]["strip_outputs"] == [{"key": "E"}]
    assert model.metadata["avatar_parameters"]["parameters"] == [
        {
            "tag": "parameter",
            "attributes": {"name": "AVATAR", "avatarparameters": "algorithm=0"},
            "children": [
                {
                    "tag": "stateList",
                    "attributes": {"states": "{A=[1]}", "namestates": "state_A"},
                }
            ],
        }
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


def test_executable_model_is_not_promoted_from_io_namespace():
    assert "ExecutableModel" not in dir(bt.logic.io)
    assert not hasattr(bt.logic.io, "ExecutableModel")


def test_imported_model_get_returns_available_containers_and_validates_attribute():
    model = ExecutableModel()

    assert model.get("initial_states") == {}
    assert model.get("perturbations") == {}
    assert model.get("metadata") == {}

    with pytest.raises(ValueError, match="unsupported executable model attribute"):
        model.get("unknown")


def test_executable_model_repr_summarizes_large_mappings_by_counts():
    model = ExecutableModel(
        initial_states={"state_A": bt.logic.ba.Hypercube({"A": 1})},
        perturbations={"A KO": [{"target": "A", "min": 0, "max": 0}]},
        metadata={"archive": {"files": ["model.ginml", "notes.txt"]}},
    )

    representation = repr(model)

    assert representation == (
        "ExecutableModel(boolean_network=None, influence_graph=None, "
        "initial_states(hypercubes=1), perturbations(n=1), metadata(entries=1))"
    )
    assert "Hypercube" not in representation
    assert "state_A" not in representation
    assert "model.ginml" not in representation


def test_read_zginml_requires_main_ginml_file(tmp_path):
    archive = tmp_path / "empty.zginml"

    with ZipFile(archive, "w") as zf:
        zf.writestr("GINsim-data/notes.txt", "no model")

    with pytest.raises(ValueError, match="no '.ginml' file"):
        bt.logic.io.read_zginml(archive)


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

    model = bt.logic.io.read_zginml(archive)

    assert model.initial_states["RA_pharmacological"] == {
        "RA_b1": 1,
        "RA_b2": 1,
    }
    assert model.initial_states["RA_physiological"] == {
        "RA_b1": 1,
        "RA_b2": 0,
    }


def test_read_zginml_normalizes_initial_state_component_names(tmp_path):
    archive = tmp_path / "normalized.zginml"

    with ZipFile(archive, "w") as zf:
        zf.writestr(
            "model.ginml",
            _multivalued_ginml().replace("RA", "RA-signal").strip(),
        )
        zf.writestr(
            "initialState",
            dedent("""
                <?xml version="1.0" encoding="UTF-8"?>
                <initialStates>
                  <input name="high" value="RA-signal;2"/>
                </initialStates>
                """).strip(),
        )

    model = bt.logic.io.read_zginml(archive)

    assert model.initial_states["high"] == {
        "RA_signal_b1": 1,
        "RA_signal_b2": 1,
    }
    assert set(model.initial_states["high"]) <= model.get("boolean_network").components


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
