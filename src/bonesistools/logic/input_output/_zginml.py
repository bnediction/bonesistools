#!/usr/bin/env python

"""
ZGINML writer for executable logical models.
"""

from __future__ import annotations

import warnings
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from xml.etree import ElementTree

from ..boolean_algebra import Hypercube
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph
from ._booleanization import _syntactic_influences

_DEFAULT_MAIN_GINML = "GINsim-data/regulatoryGraph.ginml"
_BOOLEANIZED_NODE_GAP = 15.0
_XLINK_NAMESPACE = "http://www.w3.org/1999/xlink"


def _write_zginml(
    file: Path,
    *,
    boolean_network: Optional[BooleanNetwork],
    influence_graph: Optional[InfluenceGraph],
    initial_conditions: Mapping[str, Hypercube],
    parameters: Mapping[str, Any],
    perturbations: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> None:
    """Write protected executable-model data to a ZGINML archive."""

    rebuild_source = _can_rebuild_source_ginml(boolean_network, metadata)
    source_content = _source_ginml_bytes(metadata) if rebuild_source else None

    if (
        source_content is not None
        and boolean_network is not None
        and not _source_matches_network(source_content, boolean_network)
    ):
        rebuild_source = False
        source_content = None

    if source_content is not None:
        main_content = source_content
    elif boolean_network is not None:
        main_content = _boolean_ginml_bytes(
            boolean_network,
            influence_graph,
            metadata,
        )
    else:
        raise ValueError(
            "cannot save executable model: no Boolean network or reusable "
            "GINML source metadata is available"
        )

    main_name = _main_archive_name(metadata)
    booleanized_source = bool(
        boolean_network is not None
        and metadata.get("model_type") == "multi-valued"
        and not rebuild_source
    )
    files = {main_name: main_content}
    files.update(
        _companion_files(
            main_name,
            initial_conditions=initial_conditions,
            parameters=parameters,
            perturbations=perturbations,
            metadata=metadata,
            booleanized_source=booleanized_source,
        )
    )

    unknown_files = metadata.get("unknown_companion_files", [])
    if unknown_files:
        warnings.warn(
            "ZGINML export omitted companion files whose contents were not "
            "retained during import: " + ", ".join(str(name) for name in unknown_files),
            UserWarning,
            stacklevel=3,
        )

    _write_archive(file, files)


def _can_rebuild_source_ginml(
    boolean_network: Optional[BooleanNetwork],
    metadata: Mapping[str, Any],
) -> bool:
    """Test whether parsed source metadata can safely rebuild the main GINML."""

    nodes = metadata.get("nodes")
    logical_parameters = metadata.get("raw_logical_parameters")

    if not isinstance(nodes, Mapping) or not isinstance(
        logical_parameters,
        Mapping,
    ):
        return False

    if boolean_network is None:
        return True

    if metadata.get("model_type") != "boolean":
        return False

    renamed = metadata.get("boolean_component_names", {})
    if not isinstance(renamed, Mapping):
        renamed = {}

    source_components = {str(renamed.get(component, component)) for component in nodes}
    return source_components == boolean_network.components


def _source_matches_network(content: bytes, network: BooleanNetwork) -> bool:
    """Test the retained source representation against the protected network."""

    from ._ginml import _parse_ginml_bytes

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        source_network, _, _ = _parse_ginml_bytes(
            content,
            source="<retained GINML metadata>",
        )

    return source_network is not None and source_network.rules == network.rules


def _source_ginml_bytes(metadata: Mapping[str, Any]) -> bytes:
    """Rebuild a parsed Boolean GINML graph from retained metadata."""

    root = ElementTree.Element("gxl")
    graph_attributes = _string_attributes(metadata.get("graph", {}))
    graph_attributes.setdefault("class", "regulatory")
    graph = ElementTree.SubElement(root, "graph", graph_attributes)

    _append_graph_attributes(graph, metadata)
    _append_style_definitions(graph, metadata)

    nodes = metadata.get("nodes", {})
    logical_parameters = metadata.get("raw_logical_parameters", {})

    if isinstance(nodes, Mapping):
        for node_id, node_metadata in nodes.items():
            if not isinstance(node_metadata, Mapping):
                continue

            node = ElementTree.SubElement(
                graph,
                "node",
                _source_node_attributes(str(node_id), node_metadata),
            )
            raw = logical_parameters.get(node_id, {})
            if not isinstance(raw, Mapping):
                raw = {}

            _append_logical_parameters(node, raw)
            _append_metadata_elements(node, node_metadata.get("annotations", []))
            _append_source_visual_setting(
                node,
                "nodevisualsetting",
                node_metadata.get("source_visual"),
                node_metadata.get("visual", {}),
                node=True,
            )
            _append_xml_fragments(node, raw.get("unsupported_children", []))

    edge_annotations = metadata.get("edge_annotations", {})
    edge_visuals = metadata.get("edge_visual_settings", {})
    edge_source_visuals = metadata.get("edge_source_visual_settings", {})
    edges = metadata.get("edges", [])

    if isinstance(edges, Sequence):
        for edge_metadata in edges:
            if not isinstance(edge_metadata, Mapping):
                continue

            edge = ElementTree.SubElement(
                graph,
                "edge",
                _string_attributes(edge_metadata),
            )
            edge_id = edge_metadata.get("id")

            if isinstance(edge_annotations, Mapping) and edge_id is not None:
                _append_metadata_elements(edge, edge_annotations.get(edge_id, []))

            if isinstance(edge_visuals, Mapping) and edge_id is not None:
                source_visual = (
                    edge_source_visuals.get(edge_id)
                    if isinstance(edge_source_visuals, Mapping)
                    else None
                )
                _append_source_visual_setting(
                    edge,
                    "edgevisualsetting",
                    source_visual,
                    edge_visuals.get(edge_id, {}),
                    node=False,
                )

    _append_metadata_elements(graph, metadata.get("annotations", []))
    _append_unsupported_sections(graph, metadata)
    return _ginml_document_bytes(root)


def _boolean_ginml_bytes(
    boolean_network: BooleanNetwork,
    influence_graph: Optional[InfluenceGraph],
    metadata: Mapping[str, Any],
) -> bytes:
    """Serialize the executable Boolean network as a Boolean GINML graph."""

    root = ElementTree.Element("gxl")
    graph_attributes = _string_attributes(metadata.get("graph", {}))
    graph_attributes["class"] = "regulatory"
    graph_attributes.setdefault("id", "bonesistools_model")
    graph_attributes["nodeorder"] = " ".join(boolean_network)
    graph = ElementTree.SubElement(root, "graph", graph_attributes)

    _append_graph_attributes(graph, metadata)
    _append_style_definitions(graph, metadata)

    source_nodes = metadata.get("nodes", {})
    raw_parameters = metadata.get("raw_logical_parameters", {})

    for component, rule in boolean_network.items():
        graph_data = _node_data(influence_graph, component)
        source_id = str(graph_data.get("ginml_component", component))
        source_metadata = (
            source_nodes.get(source_id, {}) if isinstance(source_nodes, Mapping) else {}
        )
        if not isinstance(source_metadata, Mapping):
            source_metadata = {}

        attributes = _boolean_node_attributes(
            component,
            boolean_network.rule(component),
            graph_data,
            source_metadata,
        )
        node = ElementTree.SubElement(graph, "node", attributes)

        if attributes.get("input") != "true":
            value = ElementTree.SubElement(node, "value", {"val": "1"})
            ElementTree.SubElement(
                value,
                "exp",
                {"str": boolean_network.rule(component).replace("~", "!")},
            )

        _append_metadata_elements(node, source_metadata.get("annotations", []))
        _append_visual_setting(
            node,
            "nodevisualsetting",
            _booleanized_node_visual(graph_data, metadata),
            node=True,
        )

        source_raw = (
            raw_parameters.get(source_id, {})
            if isinstance(raw_parameters, Mapping)
            else {}
        )
        if isinstance(source_raw, Mapping):
            _append_xml_fragments(node, source_raw.get("unsupported_children", []))

    _append_boolean_edges(
        graph,
        boolean_network,
        influence_graph,
        metadata,
    )
    _append_metadata_elements(graph, metadata.get("annotations", []))
    _append_unsupported_sections(graph, metadata)
    return _ginml_document_bytes(root)


def _append_boolean_edges(
    graph_element: ElementTree.Element,
    boolean_network: BooleanNetwork,
    influence_graph: Optional[InfluenceGraph],
    metadata: Mapping[str, Any],
) -> None:
    """Append all interactions needed by Boolean expressions and graph data."""

    signs_by_pair: Dict[Tuple[str, str], set] = {}
    data_by_pair: Dict[Tuple[str, str], List[Mapping[str, Any]]] = {}

    for source, target, sign in _syntactic_influences(boolean_network):
        signs_by_pair.setdefault((source, target), set()).add(sign)

    if influence_graph is not None:
        for source, target, data in influence_graph.edges(data=True):
            pair = (str(source), str(target))
            signs_by_pair.setdefault(pair, set()).add(int(data["sign"]))
            data_by_pair.setdefault(pair, []).append(data)

    edge_annotations = metadata.get("edge_annotations", {})

    for source, target in sorted(signs_by_pair):
        signs = signs_by_pair[(source, target)]
        data = data_by_pair.get((source, target), [])
        raw_id = next(
            (str(item["id"]) for item in data if item.get("id") not in [None, ""]),
            None,
        )
        edge_id = f"{source}:{target}"
        attributes = {
            "id": edge_id,
            "from": source,
            "to": target,
            "minvalue": "1",
            "sign": _ginml_sign(signs),
        }
        attributes.update(_generated_edge_attributes(data))
        edge = ElementTree.SubElement(graph_element, "edge", attributes)

        if isinstance(edge_annotations, Mapping) and raw_id is not None:
            _append_metadata_elements(edge, edge_annotations.get(raw_id, []))

        visual = _visual_from_graph_data(data[0]) if data else {}
        _append_visual_setting(
            edge,
            "edgevisualsetting",
            visual,
            node=False,
        )


def _companion_files(
    main_name: str,
    *,
    initial_conditions: Mapping[str, Hypercube],
    parameters: Mapping[str, Any],
    perturbations: Mapping[str, Any],
    metadata: Mapping[str, Any],
    booleanized_source: bool,
) -> Dict[str, bytes]:
    """Serialize recognized GINsim companion files."""

    files = {}
    initial = _initial_conditions_bytes(
        initial_conditions,
        metadata,
        booleanized_source=booleanized_source,
    )
    if initial is not None or _had_companion(metadata, "initialstate"):
        files[_companion_name(main_name, metadata, "initialState")] = (
            initial if initial is not None else _empty_xml_bytes("initialStates")
        )

    perturbation_content = _perturbations_bytes(
        perturbations,
        metadata,
        booleanized_source=booleanized_source,
    )
    if perturbation_content is not None or _had_companion(metadata, "mutant"):
        files[_companion_name(main_name, metadata, "mutant")] = (
            perturbation_content
            if perturbation_content is not None
            else _empty_perturbations_bytes()
        )

    avatar_data = parameters.get(
        "avatar_parameters",
        metadata.get("avatar_parameters"),
    )
    avatar = _avatar_parameters_bytes(avatar_data, metadata, booleanized_source)
    if avatar is not None:
        files[_companion_name(main_name, metadata, "avatar_parameters")] = avatar

    simplifier = _model_simplifier_bytes(
        metadata.get("model_simplifier"),
        metadata,
        booleanized_source,
    )
    if simplifier is not None:
        files[_companion_name(main_name, metadata, "modelSimplifier")] = simplifier

    simulation_data = parameters.get(
        "simulation_parameters",
        metadata.get("simulation_parameters"),
    )
    simulation = _simulation_parameters_bytes(
        simulation_data,
        metadata,
        booleanized_source,
    )
    if simulation is not None:
        files[_companion_name(main_name, metadata, "reg2dyn_parameters")] = simulation

    represented_parameters = {"avatar_parameters", "simulation_parameters"}
    unsupported_parameters = sorted(set(parameters) - represented_parameters)
    if unsupported_parameters:
        warnings.warn(
            "ZGINML export omitted parameters without a recognized GINsim "
            "representation: " + ", ".join(unsupported_parameters),
            UserWarning,
            stacklevel=4,
        )

    return files


def _initial_conditions_bytes(
    initial_conditions: Mapping[str, Hypercube],
    metadata: Mapping[str, Any],
    *,
    booleanized_source: bool,
) -> Optional[bytes]:
    """Serialize named Boolean hypercubes as GINsim initial conditions."""

    unparsed = metadata.get("unparsed_initial_states", {})
    if not initial_conditions and not unparsed:
        return None

    root = ElementTree.Element("initialStates")
    state_metadata = metadata.get("initial_state_metadata", [])
    details = {
        item.get("name"): item
        for item in state_metadata
        if isinstance(item, Mapping) and item.get("name") is not None
    }
    component_names = _serialized_component_names(
        metadata,
        booleanized_source=booleanized_source,
    )

    for name, condition in initial_conditions.items():
        detail = details.get(name, {})
        section = str(detail.get("section", "initialState"))
        if section not in ["initialState", "input"]:
            section = "initialState"

        original_name = detail.get("original_name", name)
        attributes = {
            "name": name if original_name is None else str(original_name),
            "value": _state_assignment(condition, component_names),
        }
        ElementTree.SubElement(root, section, attributes)

    if isinstance(unparsed, Mapping) and not booleanized_source:
        for name, data in unparsed.items():
            if not isinstance(data, Mapping) or "value" not in data:
                continue
            section = str(data.get("section", "initialState"))
            attributes = {
                "name": str(data.get("original_name", name) or ""),
                "value": str(data["value"]),
            }
            ElementTree.SubElement(root, section, attributes)

    return _xml_document_bytes(root)


def _perturbations_bytes(
    perturbations: Mapping[str, Any],
    metadata: Mapping[str, Any],
    *,
    booleanized_source: bool,
) -> Optional[bytes]:
    """Serialize recognized GINsim perturbation definitions."""

    users = metadata.get("perturbation_users", [])
    if not perturbations and not users:
        return None

    root = ElementTree.Element("perturbationConfig")
    perturbation_list = ElementTree.SubElement(root, "listOfPerturbations")

    for name, changes in perturbations.items():
        mutant = ElementTree.SubElement(
            perturbation_list,
            "mutant",
            {"name": str(name)},
        )
        for change_type, attributes in _serialized_perturbation_changes(
            changes,
            metadata,
            booleanized_source=booleanized_source,
        ):
            ElementTree.SubElement(mutant, change_type, attributes)

    user_list = ElementTree.SubElement(root, "listOfUsers")
    if isinstance(users, Iterable):
        for user in users:
            if isinstance(user, Mapping):
                ElementTree.SubElement(user_list, "user", _string_attributes(user))

    return _xml_document_bytes(root)


def _avatar_parameters_bytes(
    data: Any,
    metadata: Mapping[str, Any],
    booleanized_source: bool,
) -> Optional[bytes]:
    """Serialize parsed AVATAR parameters."""

    if not isinstance(data, Mapping):
        return None

    attributes = _top_level_attributes(data, excluded=["parameters"])
    _transform_node_order_attribute(
        attributes,
        metadata,
        booleanized_source=booleanized_source,
    )
    root = ElementTree.Element("avatarParameters", attributes)
    _append_metadata_elements(root, data.get("parameters", []))
    return _xml_document_bytes(root).replace(b"&#10;", b"\n")


def _model_simplifier_bytes(
    data: Any,
    metadata: Mapping[str, Any],
    booleanized_source: bool,
) -> Optional[bytes]:
    """Serialize parsed model-simplifier settings."""

    if not isinstance(data, Mapping):
        return None

    root = ElementTree.Element(
        "modelModifierConfig",
        _string_attributes(data.get("attributes", {})),
    )
    simplifications = ElementTree.SubElement(root, "modelSimplifications")

    for item in data.get("simplifications", []):
        if not isinstance(item, Mapping):
            continue
        attributes = _string_attributes(item)
        if booleanized_source and "removeList" in item:
            attributes["removeList"] = " ".join(
                _transform_components(item["removeList"], metadata)
            )
        ElementTree.SubElement(
            simplifications,
            "simplificationConfig",
            attributes,
        )

    for item in data.get("strip_outputs", []):
        if isinstance(item, Mapping):
            ElementTree.SubElement(root, "stripOutput", _string_attributes(item))

    users = data.get("users", [])
    if users:
        user_list = ElementTree.SubElement(root, "listOfUsers")
        for user in users:
            if isinstance(user, Mapping):
                ElementTree.SubElement(user_list, "user", _string_attributes(user))

    return _xml_document_bytes(root)


def _simulation_parameters_bytes(
    data: Any,
    metadata: Mapping[str, Any],
    booleanized_source: bool,
) -> Optional[bytes]:
    """Serialize parsed simulation parameters."""

    if not isinstance(data, Mapping):
        return None

    excluded = [
        "parameters",
        "parameter_sections",
        "priority_class_lists",
    ]
    attributes = _top_level_attributes(data, excluded=excluded)
    _transform_node_order_attribute(
        attributes,
        metadata,
        booleanized_source=booleanized_source,
    )
    root = ElementTree.Element("simulationParameters", attributes)

    for priority_list in data.get("priority_class_lists", []):
        if not isinstance(priority_list, Mapping):
            continue
        element = ElementTree.SubElement(
            root,
            "priorityClassList",
            _string_attributes(priority_list.get("attributes", {})),
        )
        for priority_class in priority_list.get("classes", []):
            if not isinstance(priority_class, Mapping):
                continue
            class_attributes = _string_attributes(priority_class)
            if booleanized_source and "content" in priority_class:
                class_attributes["content"] = " ".join(
                    _transform_priority_content(
                        priority_class["content"],
                        metadata,
                    )
                )
            ElementTree.SubElement(element, "class", class_attributes)

    parameter_sections = data.get("parameter_sections", {})
    for name, parameter in data.get("parameters", {}).items():
        if not isinstance(parameter, Mapping):
            continue
        parameter_attributes = _string_attributes(parameter.get("attributes", {}))
        parameter_attributes.setdefault("name", str(name))
        element = ElementTree.SubElement(root, "parameter", parameter_attributes)
        _append_named_rows(element, "initstates", parameter.get("initstates", []))
        _append_named_rows(element, "inputs", parameter.get("inputs", []))
        if isinstance(parameter_sections, Mapping):
            _append_metadata_elements(element, parameter_sections.get(name, []))

    return _xml_document_bytes(root)


def _serialized_perturbation_changes(
    changes: Any,
    metadata: Mapping[str, Any],
    *,
    booleanized_source: bool,
) -> List[Tuple[str, Dict[str, str]]]:
    """Normalize public and parsed perturbation schemas for XML output."""

    if not isinstance(changes, Iterable) or isinstance(changes, (str, bytes)):
        warnings.warn(
            "ZGINML export ignored a perturbation with an unsupported definition",
            UserWarning,
            stacklevel=5,
        )
        return []

    serialized = []
    for change in changes:
        if not isinstance(change, Mapping):
            continue

        change_type = str(change.get("type", "change"))
        raw_attributes = change.get("attributes", change)
        if not isinstance(raw_attributes, Mapping):
            continue
        attributes = dict(raw_attributes)
        attributes.pop("type", None)

        if "value" in attributes:
            attributes.setdefault("min", attributes["value"])
            attributes.setdefault("max", attributes["value"])
            attributes.pop("value", None)

        if booleanized_source and change_type == "change":
            converted = _booleanized_perturbation(attributes, metadata)
            if converted is not None:
                serialized.extend((change_type, item) for item in converted)
                continue

        serialized.append((change_type, _string_attributes(attributes)))

    return serialized


def _booleanized_perturbation(
    attributes: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> Optional[List[Dict[str, str]]]:
    """Translate a fixed multi-valued perturbation to threshold variables."""

    target = attributes.get("target")
    if not isinstance(target, str):
        return None

    multi_valued = metadata.get("multi_valued_components", {})
    renamed = metadata.get("boolean_component_names", {})
    if not isinstance(multi_valued, Mapping):
        multi_valued = {}
    if not isinstance(renamed, Mapping):
        renamed = {}

    normalized_target = str(renamed.get(target, target))
    if target not in multi_valued:
        converted = _string_attributes(attributes)
        converted["target"] = normalized_target
        return [converted]

    minimum_value = attributes.get("min")
    maximum_value = attributes.get("max")
    if minimum_value is None or maximum_value is None:
        return None

    minimum = int(minimum_value)
    maximum = int(maximum_value)
    max_level = int(multi_valued[target])

    if minimum == -1 and maximum == -1:
        return []

    if minimum < 0 or maximum > max_level or minimum > maximum:
        warnings.warn(
            f"ZGINML export could not Booleanize invalid perturbation range for "
            f"{target!r}; its original definition was retained",
            UserWarning,
            stacklevel=6,
        )
        return None

    converted = []
    for threshold in range(1, max_level + 1):
        if threshold <= minimum:
            value = "1"
        elif threshold > maximum:
            value = "0"
        else:
            continue

        converted.append(
            {
                "target": f"{normalized_target}_b{threshold}",
                "min": value,
                "max": value,
            }
        )

    return converted


def _source_node_attributes(
    node_id: str,
    metadata: Mapping[str, Any],
) -> Dict[str, str]:
    """Recover source node attributes with conservative fallbacks."""

    attributes = _string_attributes(metadata.get("attributes", {}))
    attributes.setdefault("id", node_id)
    attributes.setdefault("maxvalue", str(metadata.get("maxvalue", 1)))
    if metadata.get("name") is not None:
        attributes.setdefault("name", str(metadata["name"]))
    if bool(metadata.get("input")):
        attributes.setdefault("input", "true")
    return attributes


def _boolean_node_attributes(
    component: str,
    rule: str,
    graph_data: Mapping[str, Any],
    source_metadata: Mapping[str, Any],
) -> Dict[str, str]:
    """Build Boolean node attributes without multi-valued provenance fields."""

    attributes = {}
    source_attributes = source_metadata.get("attributes", {})
    if isinstance(source_attributes, Mapping):
        attributes.update(
            _string_attributes(
                {
                    key: value
                    for key, value in source_attributes.items()
                    if key not in ["basevalue", "id", "input", "maxvalue", "name"]
                }
            )
        )

    attributes.update(_generated_node_attributes(graph_data))
    attributes["id"] = component
    attributes["maxvalue"] = "1"
    name = graph_data.get("name")
    if name not in [None, "", component]:
        attributes["name"] = str(name)
    if rule == component:
        attributes["input"] = "true"
    return attributes


def _append_logical_parameters(
    node: ElementTree.Element,
    metadata: Mapping[str, Any],
) -> None:
    """Append retained value expressions and logical parameters."""

    for value_data in metadata.get("values", []):
        if not isinstance(value_data, Mapping):
            continue
        attributes = _string_attributes(value_data.get("attributes", {}))
        attributes.setdefault("val", str(value_data.get("val", 1)))
        value = ElementTree.SubElement(node, "value", attributes)
        for expression in value_data.get("expressions", []):
            ElementTree.SubElement(value, "exp", {"str": str(expression)})

    for parameter in metadata.get("parameters", []):
        if isinstance(parameter, Mapping):
            ElementTree.SubElement(node, "parameter", _string_attributes(parameter))


def _append_graph_attributes(
    graph: ElementTree.Element,
    metadata: Mapping[str, Any],
) -> None:
    """Append named graph attributes."""

    attributes = metadata.get("graph_attributes", {})
    if not isinstance(attributes, Mapping):
        return

    for name, value in attributes.items():
        ElementTree.SubElement(
            graph,
            "attr",
            {"name": str(name), "value": str(value)},
        )


def _append_style_definitions(
    graph: ElementTree.Element,
    metadata: Mapping[str, Any],
) -> None:
    """Append retained node and edge style definitions."""

    for item in metadata.get("node_style_definitions", []):
        if isinstance(item, Mapping):
            ElementTree.SubElement(graph, "nodestyle", _string_attributes(item))

    for item in metadata.get("edge_style_definitions", []):
        if isinstance(item, Mapping):
            ElementTree.SubElement(graph, "edgestyle", _string_attributes(item))


def _append_visual_setting(
    parent: ElementTree.Element,
    tag: str,
    visual: Any,
    *,
    node: bool,
) -> None:
    """Append a direct visual setting reconstructed from normalized metadata."""

    if not isinstance(visual, Mapping) or not visual:
        return

    shape_key = "shape" if node else "path"
    shape = visual.get(shape_key)
    parent_attributes = {}
    if "style" in visual and visual.get("style") is not None:
        parent_attributes["style"] = str(visual["style"])

    if shape is None:
        parent_attributes.update(
            _visual_attributes(
                visual,
                excluded=["style", shape_key],
                node=node,
            )
        )
        ElementTree.SubElement(parent, tag, parent_attributes)
        return

    setting = ElementTree.SubElement(parent, tag, parent_attributes)
    ElementTree.SubElement(
        setting,
        str(shape),
        _visual_attributes(
            visual,
            excluded=["style", shape_key],
            node=node,
        ),
    )


def _append_source_visual_setting(
    parent: ElementTree.Element,
    tag: str,
    source_visual: Any,
    visual: Any,
    *,
    node: bool,
) -> None:
    """Append retained source visual XML or rebuild normalized metadata."""

    if isinstance(source_visual, Mapping) and source_visual.get("tag") == tag:
        parent.append(_element_from_metadata(source_visual))
        return

    _append_visual_setting(parent, tag, visual, node=node)


def _append_metadata_elements(parent: ElementTree.Element, values: Any) -> None:
    """Append elements represented by generic parsed XML metadata."""

    if not isinstance(values, Iterable) or isinstance(values, (str, bytes)):
        return

    for value in values:
        if isinstance(value, Mapping) and "tag" in value:
            parent.append(_element_from_metadata(value))


def _append_unsupported_sections(
    graph: ElementTree.Element,
    metadata: Mapping[str, Any],
) -> None:
    """Preserve retained unsupported graph sections as XML fragments."""

    sections = metadata.get("unsupported_sections", [])
    if not isinstance(sections, Iterable):
        return

    for section in sections:
        if isinstance(section, Mapping) and "xml" in section:
            _append_xml_fragments(graph, [section["xml"]])


def _append_xml_fragments(parent: ElementTree.Element, fragments: Any) -> None:
    """Append independently retained XML fragments when they remain valid."""

    if not isinstance(fragments, Iterable) or isinstance(fragments, (str, bytes)):
        return

    for fragment in fragments:
        try:
            parent.append(ElementTree.fromstring(str(fragment)))
        except ElementTree.ParseError:
            warnings.warn(
                "ZGINML export omitted malformed retained XML metadata",
                UserWarning,
                stacklevel=5,
            )


def _element_from_metadata(metadata: Mapping[str, Any]) -> ElementTree.Element:
    """Recreate one XML element from generic parsed metadata."""

    tag = str(metadata["tag"])
    attributes = {}
    raw_attributes = metadata.get("attributes", {})
    if isinstance(raw_attributes, Mapping):
        for name, value in raw_attributes.items():
            xml_name = str(name)
            if tag == "link" and xml_name == "href":
                xml_name = f"{{{_XLINK_NAMESPACE}}}href"
            attributes[xml_name] = str(value)

    element = ElementTree.Element(tag, attributes)
    if metadata.get("text") is not None:
        element.text = str(metadata["text"])
    _append_metadata_elements(element, metadata.get("children", []))
    return element


def _serialized_component_names(
    metadata: Mapping[str, Any],
    *,
    booleanized_source: bool,
) -> Dict[str, str]:
    """Map Boolean component names to names used by the emitted main graph."""

    if booleanized_source:
        return {}

    renamed = metadata.get("boolean_component_names", {})
    if not isinstance(renamed, Mapping):
        return {}
    return {str(target): str(source) for source, target in renamed.items()}


def _state_assignment(
    condition: Hypercube,
    component_names: Mapping[str, str],
) -> str:
    """Serialize fixed coordinates of a Boolean hypercube."""

    return " ".join(
        f"{component_names.get(component, component)};{condition[component].value}"
        for component in condition
        if condition[component].is_fixed
    )


def _transform_node_order_attribute(
    attributes: Dict[str, str],
    metadata: Mapping[str, Any],
    *,
    booleanized_source: bool,
) -> None:
    """Expand multi-valued node orders to emitted threshold components."""

    if not booleanized_source or "nodeOrder" not in attributes:
        return
    attributes["nodeOrder"] = " ".join(
        _transform_components(attributes["nodeOrder"].split(), metadata)
    )


def _transform_components(values: Any, metadata: Mapping[str, Any]) -> List[str]:
    """Expand source components to normalized Boolean threshold components."""

    if isinstance(values, str):
        values = values.split()

    multi_valued = metadata.get("multi_valued_components", {})
    renamed = metadata.get("boolean_component_names", {})
    if not isinstance(multi_valued, Mapping):
        multi_valued = {}
    if not isinstance(renamed, Mapping):
        renamed = {}

    transformed = []
    for value in values:
        source = str(value)
        component = str(renamed.get(source, source))
        maximum = int(multi_valued.get(source, 1))
        if maximum == 1:
            transformed.append(component)
        else:
            transformed.extend(
                f"{component}_b{threshold}" for threshold in range(1, maximum + 1)
            )
    return transformed


def _transform_priority_content(
    values: Any,
    metadata: Mapping[str, Any],
) -> List[str]:
    """Expand component names while retaining priority direction suffixes."""

    if isinstance(values, str):
        values = values.split()

    transformed = []
    for value in values:
        text = str(value)
        component, separator, suffix = text.partition(",")
        for expanded in _transform_components([component], metadata):
            transformed.append(expanded if not separator else f"{expanded},{suffix}")
    return transformed


def _append_named_rows(
    parent: ElementTree.Element,
    section: str,
    names: Any,
) -> None:
    """Append a GINsim row-list section."""

    if not names:
        return
    section_element = ElementTree.SubElement(parent, section)
    for name in names:
        ElementTree.SubElement(section_element, "row", {"name": str(name)})


def _top_level_attributes(
    metadata: Mapping[str, Any],
    *,
    excluded: Iterable[str],
) -> Dict[str, str]:
    """Serialize scalar root attributes while excluding parsed child fields."""

    excluded = set(excluded)
    return _string_attributes(
        {key: value for key, value in metadata.items() if key not in excluded}
    )


def _generated_node_attributes(data: Mapping[str, Any]) -> Dict[str, str]:
    """Recover custom source attributes stored with a Boolean graph node."""

    reserved = {
        "basevalue",
        "component",
        "maxvalue",
        "name",
        "threshold",
        "input",
        "style",
        "shape",
        "x",
        "y",
        "width",
        "height",
        "background",
        "foreground",
        "text",
    }
    return {
        key[len("ginml_") :]: str(value)
        for key, value in data.items()
        if key.startswith("ginml_")
        and key[len("ginml_") :] not in reserved
        and not isinstance(value, (Mapping, list, tuple, set))
    }


def _generated_edge_attributes(
    values: Sequence[Mapping[str, Any]],
) -> Dict[str, str]:
    """Recover non-structural GINML attributes from one emitted interaction."""

    if not values:
        return {}

    reserved = {
        "style",
        "path",
        "points",
        "anchor",
        "color",
        "line_width",
        "line_style",
        "routage",
        "generated",
        "interactions",
    }
    return {
        key[len("ginml_") :]: str(value)
        for key, value in values[0].items()
        if key.startswith("ginml_")
        and key[len("ginml_") :] not in reserved
        and not isinstance(value, (Mapping, list, tuple, set))
    }


def _booleanized_node_visual(
    data: Mapping[str, Any],
    metadata: Mapping[str, Any],
) -> Dict[str, Any]:
    """Adapt source visuals for threshold-expanded Boolean components."""

    visual = _visual_from_graph_data(data)
    if metadata.get("model_type") != "multi-valued":
        return visual

    if visual.get("shape") == "point":
        visual["shape"] = "ellipse"
        visual.setdefault("width", "45")
        visual.setdefault("height", "25")
        visual.setdefault("background", "#ffffff")
        visual.setdefault("foreground", "#000000")

    source = data.get("ginml_component")
    threshold = data.get("ginml_threshold")
    multi_valued = metadata.get("multi_valued_components", {})
    if source is None or threshold is None or not isinstance(multi_valued, Mapping):
        return visual

    try:
        threshold = int(threshold)
        maximum = int(multi_valued[source])
        x = float(visual["x"])
        width = float(visual.get("width", 45.0))
    except (KeyError, TypeError, ValueError):
        return visual

    if maximum <= 1 or threshold <= 1:
        return visual

    x += (threshold - 1) * (width + _BOOLEANIZED_NODE_GAP)
    visual["x"] = str(int(x)) if x.is_integer() else str(x)
    return visual


def _visual_from_graph_data(data: Mapping[str, Any]) -> Dict[str, Any]:
    """Extract normalized GINML visual attributes from graph data."""

    visual_keys = {
        "style",
        "shape",
        "path",
        "x",
        "y",
        "width",
        "height",
        "background",
        "foreground",
        "text",
        "color",
        "line_width",
        "line_style",
        "routage",
        "points",
        "anchor",
        "properties",
    }
    return {
        key[len("ginml_") :]: value
        for key, value in data.items()
        if key.startswith("ginml_") and key[len("ginml_") :] in visual_keys
    }


def _node_data(
    influence_graph: Optional[InfluenceGraph],
    component: str,
) -> Mapping[str, Any]:
    """Return graph-node data when available."""

    if influence_graph is None or component not in influence_graph:
        return {}
    return influence_graph.nodes[component]


def _ginml_sign(signs: Iterable[int]) -> str:
    """Convert one or two signed influences to a GINML sign."""

    signs = set(signs)
    if signs == {1}:
        return "positive"
    if signs == {-1}:
        return "negative"
    return "dual"


def _visual_attributes(
    visual: Mapping[str, Any],
    *,
    excluded: Iterable[str],
    node: bool,
) -> Dict[str, str]:
    """Convert normalized visual keys back to GINML XML attributes."""

    excluded = set(excluded)
    renamed = {
        "background": "backgroundColor",
        "foreground": "foregroundColor",
        "text": "textColor",
        "color": "foregroundColor" if node else "line_color",
    }
    return {
        renamed.get(key, key): str(value)
        for key, value in visual.items()
        if key not in excluded and value is not None
    }


def _string_attributes(values: Any) -> Dict[str, str]:
    """Convert a mapping to XML-compatible string attributes."""

    if not isinstance(values, Mapping):
        return {}

    attributes = {}
    for key, value in values.items():
        if value is None:
            continue
        if isinstance(value, bool):
            attributes[str(key)] = "true" if value else "false"
        elif isinstance(value, (list, tuple)):
            attributes[str(key)] = " ".join(str(item) for item in value)
        elif isinstance(value, (Mapping, set)):
            continue
        else:
            attributes[str(key)] = str(value)
    return attributes


def _main_archive_name(metadata: Mapping[str, Any]) -> str:
    """Reuse a safe source archive path for the main GINML when available."""

    archive = metadata.get("archive", {})
    if isinstance(archive, Mapping):
        candidate = archive.get("main_ginml")
        if isinstance(candidate, str) and candidate.lower().endswith(".ginml"):
            safe = _safe_archive_name(candidate)
            if safe is not None:
                return safe
    return _DEFAULT_MAIN_GINML


def _companion_name(
    main_name: str,
    metadata: Mapping[str, Any],
    basename: str,
) -> str:
    """Reuse a source companion path or place it beside the main GINML."""

    archive = metadata.get("archive", {})
    files = archive.get("files", []) if isinstance(archive, Mapping) else []
    for candidate in files:
        if Path(str(candidate)).name.lower() == basename.lower():
            safe = _safe_archive_name(str(candidate))
            if safe is not None:
                return safe

    parent = PurePosixPath(main_name).parent
    return str(parent / basename) if str(parent) != "." else basename


def _had_companion(metadata: Mapping[str, Any], basename: str) -> bool:
    """Test whether the source archive contained a named companion file."""

    archive = metadata.get("archive", {})
    files = archive.get("files", []) if isinstance(archive, Mapping) else []
    return any(Path(str(name)).name.lower() == basename.lower() for name in files)


def _safe_archive_name(name: str) -> Optional[str]:
    """Reject absolute and parent-traversing archive member names."""

    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        return None
    return str(path)


def _write_archive(file: Path, files: Mapping[str, bytes]) -> None:
    """Write deterministic compressed archive entries."""

    with zipfile.ZipFile(file, "w") as archive:
        for name in sorted(files):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, files[name])


def _ginml_document_bytes(root: ElementTree.Element) -> bytes:
    """Serialize GINML XML with its standard document declaration."""

    ElementTree.register_namespace("xlink", _XLINK_NAMESPACE)
    _indent_xml(root)
    body = ElementTree.tostring(root, encoding="utf-8")
    return (
        b'<?xml version="1.0" encoding="UTF-8"?>\n'
        b'<!DOCTYPE gxl SYSTEM "http://ginsim.org/GINML_2_2.dtd">\n' + body + b"\n"
    )


def _xml_document_bytes(root: ElementTree.Element) -> bytes:
    """Serialize a companion XML document."""

    _indent_xml(root)
    body = ElementTree.tostring(root, encoding="utf-8")
    return b'<?xml version="1.0" encoding="UTF-8"?>\n' + body + b"\n"


def _empty_xml_bytes(tag: str) -> bytes:
    """Serialize an empty companion root."""

    return _xml_document_bytes(ElementTree.Element(tag))


def _empty_perturbations_bytes() -> bytes:
    """Serialize an empty GINsim perturbation configuration."""

    root = ElementTree.Element("perturbationConfig")
    ElementTree.SubElement(root, "listOfPerturbations")
    ElementTree.SubElement(root, "listOfUsers")
    return _xml_document_bytes(root)


def _indent_xml(element: ElementTree.Element, level: int = 0) -> None:
    """Indent an XML tree without requiring Python 3.9 ElementTree APIs."""

    indentation = "\n" + "  " * level
    child_indentation = "\n" + "  " * (level + 1)

    if len(element):
        if not element.text or not element.text.strip():
            element.text = child_indentation
        for child in element:
            _indent_xml(child, level + 1)
        if not element[-1].tail or not element[-1].tail.strip():
            element[-1].tail = indentation

    if level and (not element.tail or not element.tail.strip()):
        element.tail = indentation
