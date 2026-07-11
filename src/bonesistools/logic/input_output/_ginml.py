#!/usr/bin/env python

"""
GINML and ZGINML readers.
"""

from __future__ import annotations

import re
import warnings
import zipfile
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)
from xml.etree import ElementTree

from ..boolean_algebra import Hypercube
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph
from ._executable_model import ExecutableModel


def read_ginml(path: Union[str, Path]) -> ExecutableModel:
    """
    Read a GINML logical model.

    GINML import is conservative. Boolean rules are converted to a
    `BooleanNetwork` only when the logical parameters are safely representable.
    When available, its Boolean components and rules define the returned
    `InfluenceGraph`; otherwise, signed GINML interactions are used directly.
    Component characters reserved by the Boolean-rule syntax are replaced by
    underscores; changed names are recorded in
    `model.metadata["boolean_component_names"]`.
    """

    path = Path(path)
    content = path.read_bytes()
    return _read_ginml_bytes(
        content,
        source=str(path),
    )


def read_zginml(path: Union[str, Path]) -> ExecutableModel:
    """
    Read a ZGINML archive.

    The main `.ginml` file is parsed with `read_ginml` logic. Companion files
    are parsed conservatively when recognized. Unknown companion file names are
    preserved in metadata.
    """

    path = Path(path)

    with zipfile.ZipFile(path) as archive:
        archive_files = {
            name: archive.read(name)
            for name in sorted(archive.namelist())
            if not name.endswith("/")
        }

    ginml_files = [name for name in archive_files if name.lower().endswith(".ginml")]

    if not ginml_files:
        raise ValueError(f"invalid ZGINML archive: no '.ginml' file found in {path}")

    main_ginml = _select_main_ginml_file(ginml_files)
    model = _read_ginml_bytes(
        archive_files[main_ginml],
        source=str(path),
    )

    model.metadata["format"] = "zginml"
    model.metadata["archive"] = {
        "path": str(path),
        "main_ginml": main_ginml,
        "files": sorted(archive_files),
    }

    _parse_zginml_companion_files(model, archive_files, main_ginml)

    return model


def _read_ginml_bytes(
    content: bytes,
    *,
    source: str,
) -> ExecutableModel:
    root = ElementTree.fromstring(content.strip())
    graph_element = _find_graph_element(root)

    if graph_element is None:
        raise ValueError("invalid GINML file: missing graph element")

    node_styles = _parse_ginml_node_styles(graph_element)
    edge_styles = _parse_ginml_edge_styles(graph_element)
    nodes = _parse_ginml_nodes(graph_element, node_styles=node_styles)
    edges = _parse_ginml_edges(graph_element, edge_styles=edge_styles)
    metadata = _ginml_metadata(
        source=source,
        graph_element=graph_element,
        node_styles=node_styles,
        edge_styles=edge_styles,
        nodes=nodes,
        edges=edges,
    )

    boolean_network = _build_boolean_network(nodes, edges, metadata)
    influence_graph = _build_influence_graph(
        nodes,
        edges,
        metadata,
        boolean_network=boolean_network,
    )

    return ExecutableModel(
        boolean_network=boolean_network,
        influence_graph=influence_graph,
        metadata=metadata,
    )


def _find_graph_element(root: ElementTree.Element) -> Optional[ElementTree.Element]:
    if _local_name(root.tag) == "graph":
        return root

    for element in root.iter():
        if _local_name(element.tag) == "graph":
            return element

    return None


def _parse_ginml_nodes(
    graph_element: ElementTree.Element,
    *,
    node_styles: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    nodes = {}

    for element in _children(graph_element, "node"):
        node_id = element.get("id")

        if node_id is None:
            raise ValueError("invalid GINML file: node without 'id' attribute")

        maxvalue = _parse_int_attribute(element, "maxvalue", default=1)
        input_node = _parse_bool_attribute(element, "input", default=False)
        values = []
        parameters = []
        unsupported_children = []
        visual_settings: Dict[str, Any] = {}
        annotations = []

        for child in list(element):
            child_name = _local_name(child.tag)

            if child_name == "value":
                values.append(
                    {
                        "val": _parse_int_attribute(child, "val", default=1),
                        "expressions": [
                            str(exp.get("str"))
                            for exp in _children(child, "exp")
                            if exp.get("str") is not None
                        ],
                        "attributes": dict(child.attrib),
                    }
                )

            elif child_name == "parameter":
                parameters.append(dict(child.attrib))

            elif child_name == "nodevisualsetting":
                visual_settings = _parse_ginml_node_visual_settings(
                    child,
                    node_styles=node_styles,
                )

            elif child_name == "annotation":
                annotations.append(_ginml_xml_metadata(child))

            else:
                unsupported_children.append(
                    ElementTree.tostring(child, encoding="unicode")
                )

        nodes[node_id] = {
            "id": node_id,
            "name": element.get("name"),
            "maxvalue": maxvalue,
            "input": input_node,
            "attributes": dict(element.attrib),
            "values": values,
            "parameters": parameters,
            "visual": visual_settings,
            "annotations": annotations,
            "unsupported_children": unsupported_children,
        }

    return nodes


def _parse_ginml_node_styles(
    graph_element: ElementTree.Element,
) -> Dict[str, Dict[str, str]]:
    node_styles = {}

    for element in _children(graph_element, "nodestyle"):
        attributes = dict(element.attrib)
        name = attributes.pop("name", None)
        node_styles["default" if name is None else name] = attributes

    return node_styles


def _parse_ginml_edge_styles(
    graph_element: ElementTree.Element,
) -> Dict[str, Dict[str, str]]:
    edge_styles = {}

    for element in _children(graph_element, "edgestyle"):
        attributes = dict(element.attrib)
        name = attributes.pop("name", None)
        edge_styles["default" if name is None else name] = attributes

    return edge_styles


def _parse_ginml_node_visual_settings(
    element: ElementTree.Element,
    *,
    node_styles: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    attributes = dict(element.attrib)
    style = attributes.get("style")
    resolved: Dict[str, Any] = {}

    if "default" in node_styles:
        resolved.update(node_styles["default"])

    if style is not None and style != "" and style in node_styles:
        resolved.update(node_styles[style])

    resolved.update(_normalize_ginml_visual_attributes(attributes))

    for child in list(element):
        shape = _local_name(child.tag)
        resolved["shape"] = shape
        resolved.update(_normalize_ginml_visual_attributes(dict(child.attrib)))

    return resolved


def _normalize_ginml_visual_attributes(attributes: Mapping[str, Any]) -> Dict[str, Any]:
    normalized = {}

    for key, value in attributes.items():
        if key == "backgroundColor":
            normalized["background"] = value
        elif key == "foregroundColor":
            normalized["foreground"] = value
        elif key == "textColor":
            normalized["text"] = value
        elif key == "line_color":
            normalized["color"] = value
        else:
            normalized[key] = value

    return normalized


def _parse_ginml_edges(
    graph_element: ElementTree.Element,
    *,
    edge_styles: Mapping[str, Mapping[str, Any]],
) -> List[Dict[str, Any]]:
    edges = []

    for element in _children(graph_element, "edge"):
        source = element.get("from")
        target = element.get("to")

        if source is None or target is None:
            raise ValueError("invalid GINML file: edge without 'from' or 'to'")

        visual_settings: Dict[str, Any] = {}
        annotations = []
        unsupported_children = []

        for child in list(element):
            child_name = _local_name(child.tag)

            if child_name == "edgevisualsetting":
                visual_settings = _parse_ginml_edge_visual_settings(
                    child,
                    edge=element.attrib,
                    edge_styles=edge_styles,
                )
            elif child_name == "annotation":
                annotations.append(_ginml_xml_metadata(child))
            else:
                unsupported_children.append(
                    ElementTree.tostring(child, encoding="unicode")
                )

        edges.append(
            {
                "id": element.get("id"),
                "source": source,
                "target": target,
                "sign": element.get("sign"),
                "attributes": dict(element.attrib),
                "visual": visual_settings,
                "annotations": annotations,
                "unsupported_children": unsupported_children,
            }
        )

    return edges


def _parse_ginml_edge_visual_settings(
    element: ElementTree.Element,
    *,
    edge: Mapping[str, Any],
    edge_styles: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    attributes = dict(element.attrib)
    style = attributes.get("style")
    resolved: Dict[str, Any] = {}

    if "default" in edge_styles:
        resolved.update(edge_styles["default"])
        _apply_ginml_edge_sign_color(resolved, edge)

    if style is not None and style != "" and style in edge_styles:
        resolved.update(edge_styles[style])
        _apply_ginml_edge_sign_color(resolved, edge)

    resolved.update(_normalize_ginml_visual_attributes(attributes))

    for child in list(element):
        resolved["path"] = _local_name(child.tag)
        resolved.update(_normalize_ginml_visual_attributes(dict(child.attrib)))

    return resolved


def _apply_ginml_edge_sign_color(
    visual: Dict[str, Any],
    edge: Mapping[str, Any],
) -> None:
    properties = visual.get("properties")

    if not isinstance(properties, str):
        return

    colors = {}

    for item in properties.split():
        if ":" not in item:
            continue

        sign, color = item.split(":", maxsplit=1)
        colors[sign] = color

    signs = _ginml_edge_signs({"sign": edge.get("sign"), "attributes": edge})

    if signs == [1]:
        sign_name = "positive"
    elif signs == [-1]:
        sign_name = "negative"
    elif set(signs) == {-1, 1}:
        sign_name = "dual"
    else:
        return

    if sign_name in colors:
        visual["color"] = colors[sign_name]


def _ginml_metadata(
    *,
    source: str,
    graph_element: ElementTree.Element,
    node_styles: Mapping[str, Mapping[str, Any]],
    edge_styles: Mapping[str, Mapping[str, Any]],
    nodes: Mapping[str, Mapping[str, Any]],
    edges: List[Dict[str, Any]],
) -> Dict[str, Any]:
    multi_valued_components = {
        node_id: int(node["maxvalue"])
        for node_id, node in nodes.items()
        if int(node["maxvalue"]) > 1
    }

    return {
        "format": "ginml",
        "source": source,
        "graph": dict(graph_element.attrib),
        "graph_attributes": _ginml_graph_attributes(graph_element),
        "annotations": [
            _ginml_xml_metadata(element)
            for element in _children(graph_element, "annotation")
        ],
        "node_order": _parse_node_order(graph_element.get("nodeorder")),
        "node_styles": dict(node_styles),
        "node_style_definitions": [
            dict(element.attrib) for element in _children(graph_element, "nodestyle")
        ],
        "node_style_groups": _ginml_node_style_groups(nodes),
        "node_visual_settings": {
            node_id: node["visual"] for node_id, node in nodes.items() if node["visual"]
        },
        "edge_styles": dict(edge_styles),
        "edge_style_definitions": [
            dict(element.attrib) for element in _children(graph_element, "edgestyle")
        ],
        "edge_style_groups": _ginml_edge_style_groups(edges),
        "edge_visual_settings": {
            str(edge["id"]): edge["visual"]
            for edge in edges
            if edge["id"] is not None and edge["visual"]
        },
        "edge_annotations": {
            str(edge["id"]): edge["annotations"]
            for edge in edges
            if edge["id"] is not None and edge["annotations"]
        },
        "nodes": {
            node_id: {
                "name": node["name"],
                "maxvalue": node["maxvalue"],
                "input": node["input"],
                "attributes": node["attributes"],
                "visual": node["visual"],
                "annotations": node["annotations"],
            }
            for node_id, node in nodes.items()
        },
        "edges": [edge["attributes"] for edge in edges],
        "model_type": "multi-valued" if multi_valued_components else "boolean",
        "multi_valued_components": multi_valued_components,
        "raw_logical_parameters": {
            node_id: {
                "values": node["values"],
                "parameters": node["parameters"],
                "unsupported_children": node["unsupported_children"],
            }
            for node_id, node in nodes.items()
        },
        "unsupported_sections": _unsupported_ginml_graph_sections(graph_element),
    }


def _ginml_graph_attributes(
    graph_element: ElementTree.Element,
) -> Dict[str, str]:
    attributes = {}

    for element in _children(graph_element, "attr"):
        name = element.get("name")
        value = element.get("value")

        if name is not None and value is not None:
            attributes[name] = value

    return attributes


def _unsupported_ginml_graph_sections(
    graph_element: ElementTree.Element,
) -> List[Dict[str, Any]]:
    supported = {"annotation", "attr", "edge", "edgestyle", "node", "nodestyle"}
    unsupported = []

    for element in list(graph_element):
        if _local_name(element.tag) in supported:
            continue

        unsupported.append(
            {
                "section": _local_name(element.tag),
                "reason": "unsupported GINML graph section",
                "xml": ElementTree.tostring(element, encoding="unicode"),
            }
        )

    return unsupported


def _ginml_node_style_groups(
    nodes: Mapping[str, Mapping[str, Any]],
) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}

    for node_id, node in nodes.items():
        style = node["visual"].get("style")

        if style in [None, ""]:
            continue

        groups.setdefault(str(style), []).append(node_id)

    return groups


def _ginml_edge_style_groups(
    edges: Iterable[Mapping[str, Any]],
) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}

    for edge in edges:
        style = edge["visual"].get("style")
        edge_id = edge.get("id")

        if style in [None, ""] or edge_id is None:
            continue

        groups.setdefault(str(style), []).append(str(edge_id))

    return groups


def _ginml_visual_attributes(visual: Mapping[str, Any]) -> Dict[str, Any]:
    attributes = {}

    for key, value in visual.items():
        if value is None:
            continue

        if key == "style" and value == "":
            continue

        attributes[f"ginml_{key}"] = value

    return attributes


def _ginml_source_attributes(
    attributes: Mapping[str, Any],
    *,
    excluded: Iterable[str],
) -> Dict[str, Any]:
    excluded = set(excluded)

    return {
        f"ginml_{_ginml_attribute_name(key)}": value
        for key, value in attributes.items()
        if key not in excluded and value is not None
    }


def _ginml_attribute_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


def _build_influence_graph(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: List[Dict[str, Any]],
    metadata: Dict[str, Any],
    *,
    boolean_network: Optional[BooleanNetwork],
) -> Optional[InfluenceGraph]:
    if boolean_network is not None:
        return _build_booleanized_influence_graph(
            nodes,
            edges,
            metadata,
            boolean_network,
        )

    graph = InfluenceGraph()

    graph.graph.update(
        _ginml_source_attributes(metadata["graph"], excluded=[]),
    )
    graph.graph.update(
        _ginml_source_attributes(metadata["graph_attributes"], excluded=[]),
    )

    for node_id, node in nodes.items():
        graph.add_node(
            node_id,
            name=node["name"],
            maxvalue=node["maxvalue"],
            input=node["input"],
            **_ginml_source_attributes(
                node["attributes"],
                excluded=["id", "input", "maxvalue", "name"],
            ),
            **_ginml_visual_attributes(node["visual"]),
        )

    skipped_edges = []

    for edge in edges:
        signs = _ginml_edge_signs(edge)

        if not signs:
            skipped_edges.append(edge["attributes"])
            continue

        for sign in signs:
            edge_attr = {
                key: value
                for key, value in edge["attributes"].items()
                if key not in ["from", "to", "sign"]
            }
            edge_attr["raw_sign"] = edge.get("sign")
            edge_attr.update(_ginml_visual_attributes(edge["visual"]))

            try:
                graph.add_edge(
                    edge["source"],
                    edge["target"],
                    sign=cast(Any, sign),
                    **edge_attr,
                )
            except ValueError as error:
                skipped_edges.append(
                    {
                        "edge": edge["attributes"],
                        "reason": str(error),
                    }
                )

    if skipped_edges:
        metadata["unsupported_sections"].append(
            {
                "section": "edges",
                "reason": "some interactions are unsigned, unsupported or duplicated",
                "items": skipped_edges,
            }
        )

    if graph.number_of_edges() == 0:
        metadata["influence_graph_reason"] = "no signed interactions found"
        return None

    return graph


def _build_booleanized_influence_graph(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: List[Dict[str, Any]],
    metadata: Mapping[str, Any],
    boolean_network: BooleanNetwork,
) -> InfluenceGraph:
    graph = InfluenceGraph()
    component_names = _boolean_component_names(nodes)
    component_origins = {}

    graph.graph.update(
        _ginml_source_attributes(metadata["graph"], excluded=[]),
    )
    graph.graph.update(
        _ginml_source_attributes(metadata["graph_attributes"], excluded=[]),
    )

    for node_id, node in nodes.items():
        component = component_names[node_id]
        maxvalue = int(node["maxvalue"])

        if maxvalue == 1:
            boolean_components = [(component, None)]
        else:
            boolean_components = [
                (_threshold_component(component, threshold), threshold)
                for threshold in range(1, maxvalue + 1)
            ]

        for boolean_component, threshold in boolean_components:
            component_origins[boolean_component] = node_id
            attributes = _booleanized_node_attributes(
                node_id,
                node,
                boolean_component=boolean_component,
                threshold=threshold,
            )
            graph.add_node(boolean_component, **attributes)

    edges_by_components: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}

    for edge in edges:
        key = (str(edge["source"]), str(edge["target"]))
        edges_by_components.setdefault(key, []).append(edge)

    for source, target, sign in _syntactic_boolean_influences(boolean_network):
        original_source = component_origins[source]
        original_target = component_origins[target]
        candidates = edges_by_components.get(
            (original_source, original_target),
            [],
        )
        matching_candidates = [
            edge for edge in candidates if sign in _ginml_edge_signs(edge)
        ]

        if matching_candidates:
            candidates = matching_candidates

        attributes = _booleanized_edge_attributes(candidates)

        if not candidates:
            attributes["ginml_generated"] = True

        graph.add_edge(source, target, sign=cast(Any, sign), **attributes)

    return graph


def _booleanized_node_attributes(
    node_id: str,
    node: Mapping[str, Any],
    *,
    boolean_component: str,
    threshold: Optional[int],
) -> Dict[str, Any]:
    attributes = {
        "name": node["name"] if threshold is None else boolean_component,
        "maxvalue": 1,
        "input": node["input"],
        **_ginml_source_attributes(
            node["attributes"],
            excluded=["id", "input", "maxvalue", "name"],
        ),
        **_ginml_visual_attributes(node["visual"]),
    }

    if boolean_component != node_id or threshold is not None:
        attributes["ginml_component"] = node_id

    if threshold is not None:
        attributes["ginml_name"] = node["name"]
        attributes["ginml_maxvalue"] = node["maxvalue"]
        attributes["ginml_threshold"] = threshold

    return attributes


def _booleanized_edge_attributes(
    candidates: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    if len(candidates) == 1:
        edge = candidates[0]
        attributes = {
            key: value
            for key, value in edge["attributes"].items()
            if key not in ["from", "to", "sign"]
        }
        attributes["raw_sign"] = edge.get("sign")
        attributes.update(_ginml_visual_attributes(edge["visual"]))
        return attributes

    if not candidates:
        return {}

    return {
        "ginml_interactions": [
            {
                **dict(edge["attributes"]),
                "visual": dict(edge["visual"]),
            }
            for edge in candidates
        ]
    }


def _syntactic_boolean_influences(
    boolean_network: BooleanNetwork,
) -> List[Tuple[str, str, int]]:
    influences = set()

    def collect(expression: Any, target: str, sign: int = 1) -> None:
        if isinstance(expression, boolean_network.ba.Symbol):
            influences.add((str(expression.obj), target, sign))
            return

        if isinstance(expression, boolean_network.ba.NOT):
            collect(expression.args[0], target, -sign)
            return

        for argument in expression.args:
            collect(argument, target, sign)

    for target, rule in boolean_network.items():
        collect(rule, target)

    return sorted(influences)


def _build_boolean_network(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: List[Dict[str, Any]],
    metadata: Dict[str, Any],
) -> Optional[BooleanNetwork]:
    try:
        component_names = _boolean_component_names(nodes)
        rules = _boolean_rules_from_ginml_nodes(nodes, edges, component_names)
        renamed_components = {
            source: target
            for source, target in component_names.items()
            if source != target
        }

        if renamed_components:
            metadata["boolean_component_names"] = renamed_components

        return BooleanNetwork(rules)

    except ValueError as error:
        metadata["boolean_network_reason"] = str(error)
        warnings.warn(
            "GINML Boolean-network conversion failed: "
            f"{error}. Returning the original signed GINML influence graph "
            "when available.",
            UserWarning,
            stacklevel=4,
        )
        return None


def _boolean_rules_from_ginml_nodes(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: Iterable[Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> Dict[str, str]:
    rules = {}
    incoming_edges: Dict[str, List[Mapping[str, Any]]] = {
        node_id: [] for node_id in nodes
    }

    for edge in edges:
        incoming_edges.setdefault(str(edge["target"]), []).append(edge)

    for node_id, node in nodes.items():
        maxvalue = int(node["maxvalue"])
        component = component_names[node_id]

        if maxvalue == 1:
            rules[component] = _boolean_rule_for_native_node(
                node_id,
                node,
                nodes,
                incoming_edges[node_id],
                component_names,
            )
            continue

        for threshold in range(1, maxvalue + 1):
            rules[_threshold_component(component, threshold)] = (
                _boolean_rule_for_threshold_node(
                    node_id,
                    threshold,
                    node,
                    nodes,
                    incoming_edges[node_id],
                    component_names,
                )
            )

    return rules


def _boolean_component_names(
    nodes: Mapping[str, Mapping[str, Any]],
) -> Dict[str, str]:
    component_names = {
        node_id: re.sub(r"[^A-Za-z0-9_]", "_", node_id) for node_id in nodes
    }

    for node_id, component in component_names.items():
        if not component or component[0].isdigit():
            component_names[node_id] = f"_{component}"

    duplicates = [
        component
        for component in set(component_names.values())
        if list(component_names.values()).count(component) > 1
    ]

    if duplicates:
        raise ValueError(
            "Boolean component-name normalization creates duplicate names: "
            + ", ".join(sorted(duplicates))
        )

    return component_names


def _boolean_rule_for_native_node(
    node_id: str,
    node: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    incoming_edges: List[Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    if bool(node["input"]):
        return component_names[node_id]

    expressions = _expressions_for_value(node, 1)
    expression_rule = "0"

    if expressions:
        converted = [
            _convert_ginml_expression(
                expression,
                nodes,
                component_names,
                target=node_id,
                incoming_edges=incoming_edges,
            )
            for expression in expressions
        ]
        expression_rule = _join_or(converted)

    return _logical_rule_from_regulatory_contexts(
        node_id,
        node,
        nodes,
        incoming_edges,
        component_names,
        threshold=1,
        fallback_rule=expression_rule,
    )


def _logical_rule_from_regulatory_contexts(
    node_id: str,
    node: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    incoming_edges: List[Mapping[str, Any]],
    component_names: Mapping[str, str],
    *,
    threshold: int,
    fallback_rule: str = "0",
) -> str:
    interactions = _ginml_interaction_conditions(
        node_id,
        incoming_edges,
        nodes,
        component_names,
    )
    interaction_ids = {interaction_id for interaction_id, _ in interactions}
    parameter_contexts = []

    for parameter in node["parameters"]:
        unsupported = set(parameter) - {"idActiveInteractions", "val"}

        if unsupported or "val" not in parameter:
            raise ValueError(
                "unsupported logical parameters for Boolean component "
                f"{node_id!r}: expected 'val' and optional "
                "'idActiveInteractions' attributes"
            )

        value = int(parameter["val"])

        if value < 0 or value > int(node["maxvalue"]):
            raise ValueError(
                f"unsupported logical parameter value {value!r} for "
                f"component {node_id!r}"
            )

        active_edge_ids = set(_parse_node_order(parameter.get("idActiveInteractions")))
        unknown_edge_ids = active_edge_ids - interaction_ids

        if unknown_edge_ids:
            unknown = ", ".join(sorted(unknown_edge_ids))
            raise ValueError(
                f"unsupported logical parameters for {node_id!r}: unknown "
                f"active interaction(s) {unknown}"
            )

        terms = []

        for interaction_id, activity in interactions:
            terms.append(
                activity
                if interaction_id in active_edge_ids
                else _negate_boolean_condition(activity)
            )

        context = _join_and(terms)

        parameter_contexts.append((value, context))

    basevalue = int(node["attributes"].get("basevalue", 0))
    maxvalue = int(node["maxvalue"])

    if basevalue < 0 or basevalue > maxvalue:
        raise ValueError(
            f"unsupported base value {basevalue!r}: expected a value between "
            f"0 and {maxvalue}"
        )

    rule = fallback_rule

    if interactions and fallback_rule == "0":
        inactive_context = _join_and(
            _negate_boolean_condition(activity) for _, activity in interactions
        )
        rule = _override_rule_in_contexts(
            rule,
            [(basevalue, inactive_context)],
            threshold=threshold,
        )
    elif not parameter_contexts and fallback_rule == "0":
        return "1" if basevalue >= threshold else "0"

    return _override_rule_in_contexts(
        rule,
        parameter_contexts,
        threshold=threshold,
    )


def _override_rule_in_contexts(
    fallback_rule: str,
    contexts: List[Tuple[int, str]],
    *,
    threshold: int,
) -> str:
    if not contexts:
        return fallback_rule

    positive_rule = _join_or(
        context for value, context in contexts if value >= threshold
    )
    negative_rule = _join_or(
        context for value, context in contexts if value < threshold
    )
    fallback_outside_negative_contexts = fallback_rule

    if fallback_rule != "0" and negative_rule != "0":
        fallback_outside_negative_contexts = (
            f"({fallback_rule}) & {_negate_boolean_condition(negative_rule)}"
        )

    if positive_rule == fallback_outside_negative_contexts:
        return positive_rule

    return _join_or(
        rule
        for rule in [fallback_outside_negative_contexts, positive_rule]
        if rule != "0"
    )


def _ginml_interaction_conditions(
    target: str,
    incoming_edges: List[Mapping[str, Any]],
    nodes: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> List[Tuple[str, str]]:
    interactions = []

    for edge in incoming_edges:
        edge_id = edge.get("id")

        if edge_id is None:
            raise ValueError(
                f"unsupported logical parameters for {target!r}: incoming "
                "interaction without an identifier"
            )

        source = str(edge["source"])
        component = component_names[source]
        maxvalue = int(nodes[source]["maxvalue"])
        attributes = edge["attributes"]
        effect_levels = _ginml_effect_levels(attributes.get("effects"))

        if effect_levels:
            interactions.extend(
                (
                    f"{edge_id}:{level}",
                    _multivalued_range_condition(
                        component,
                        minimum=level,
                        maximum=level,
                        maxvalue=maxvalue,
                    ),
                )
                for level in effect_levels
            )
            continue

        interactions.append(
            (
                str(edge_id),
                _multivalued_range_condition(
                    component,
                    minimum=int(attributes.get("minvalue", 1)),
                    maximum=int(attributes.get("maxvalue", maxvalue)),
                    maxvalue=maxvalue,
                ),
            )
        )

    return interactions


def _negate_boolean_condition(condition: str) -> str:
    if condition == "0":
        return "1"

    if condition == "1":
        return "0"

    if condition.startswith("~") and " " not in condition:
        return condition[1:]

    return f"~({condition})"


def _boolean_rule_for_threshold_node(
    node_id: str,
    threshold: int,
    node: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
    incoming_edges: List[Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    component = component_names[node_id]
    threshold_component = _threshold_component(component, threshold)

    if bool(node["input"]):
        if threshold == 1:
            return threshold_component

        return (
            f"{_threshold_component(component, threshold - 1)} "
            f"& {threshold_component}"
        )

    expressions = []

    for value in node["values"]:
        if int(value["val"]) >= threshold:
            expressions.extend(value["expressions"])

    expression_rule = "0"

    if expressions:
        converted = [
            _convert_ginml_expression(
                expression,
                nodes,
                component_names,
                target=node_id,
                incoming_edges=incoming_edges,
            )
            for expression in expressions
        ]
        expression_rule = _join_or(converted)

    desired_rule = _logical_rule_from_regulatory_contexts(
        node_id,
        node,
        nodes,
        incoming_edges,
        component_names,
        threshold=threshold,
        fallback_rule=expression_rule,
    )

    return _regularized_threshold_rule(
        desired_rule,
        component=component,
        threshold=threshold,
        maxvalue=int(node["maxvalue"]),
    )


def _regularized_threshold_rule(
    desired_rule: str,
    *,
    component: str,
    threshold: int,
    maxvalue: int,
) -> str:
    terms = []

    if desired_rule != "0":
        desired_terms = [
            _threshold_component(component, level) for level in range(1, threshold)
        ]

        if desired_rule != "1":
            desired_terms.append(f"({desired_rule})")

        terms.append(_join_and(desired_terms))

    if threshold < maxvalue:
        terms.append(
            _join_and(
                _threshold_component(component, level)
                for level in range(1, threshold + 2)
            )
        )

    return _join_or(terms)


def _expressions_for_value(node: Mapping[str, Any], value: int) -> List[str]:
    expressions = []

    for value_data in node["values"]:
        if int(value_data["val"]) == value:
            expressions.extend(value_data["expressions"])

    return expressions


def _convert_ginml_expression(
    expression: str,
    nodes: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
    *,
    target: str,
    incoming_edges: List[Mapping[str, Any]],
) -> str:
    normalized_expression = expression

    for node_id in sorted(nodes, key=len, reverse=True):
        component = component_names[node_id]

        if component == node_id:
            continue

        normalized_expression = normalized_expression.replace(node_id, component)

    nodes_by_component = {
        component: node_id for node_id, component in component_names.items()
    }

    def replace_token(match: re.Match[str]) -> str:
        token = match.group(0)

        if token in ["true", "True"]:
            return "1"

        if token in ["false", "False"]:
            return "0"

        if ":" in token:
            component, threshold_text = token.split(":", maxsplit=1)
            level = int(threshold_text)
            node_id = nodes_by_component.get(component, component)

            if node_id not in nodes:
                raise ValueError(
                    f"unsupported logical parameters: unknown component "
                    f"{component!r}"
                )

            maxvalue = int(nodes[node_id]["maxvalue"])

            if level < 1 or level > maxvalue:
                raise ValueError(
                    f"unsupported logical parameters: threshold {token!r} is "
                    f"outside component range 1..{maxvalue}"
                )

            return _ginml_regulator_condition(
                node_id,
                level=level,
                target=target,
                incoming_edges=incoming_edges,
                nodes=nodes,
                component_names=component_names,
            )

        node_id = nodes_by_component.get(token, token)

        if node_id not in nodes:
            raise ValueError(
                f"unsupported logical parameters: unknown component {token!r}"
            )

        return _ginml_regulator_condition(
            node_id,
            level=None,
            target=target,
            incoming_edges=incoming_edges,
            nodes=nodes,
            component_names=component_names,
        )

    converted = _TOKEN_PATTERN.sub(replace_token, normalized_expression)
    return converted.replace("!", "~")


def _ginml_regulator_condition(
    source: str,
    *,
    level: Optional[int],
    target: str,
    incoming_edges: List[Mapping[str, Any]],
    nodes: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    source_edges = [edge for edge in incoming_edges if edge["source"] == source]

    if not source_edges:
        raise ValueError(
            f"unsupported logical expression for {target!r}: component "
            f"{source!r} has no incoming interaction"
        )

    component = component_names[source]
    maxvalue = int(nodes[source]["maxvalue"])
    conditions = []

    for edge in source_edges:
        attributes = edge["attributes"]
        effects = _ginml_effect_levels(attributes.get("effects"))

        if effects:
            selected_levels = effects if level is None else [level]

            for selected_level in selected_levels:
                if selected_level not in effects:
                    continue

                conditions.append(
                    _multivalued_range_condition(
                        component,
                        minimum=selected_level,
                        maximum=selected_level,
                        maxvalue=maxvalue,
                    )
                )

            continue

        minimum = int(attributes.get("minvalue", 1))
        maximum = int(attributes.get("maxvalue", maxvalue))

        if level is not None and level != minimum:
            continue

        conditions.append(
            _multivalued_range_condition(
                component,
                minimum=minimum,
                maximum=maximum,
                maxvalue=maxvalue,
            )
        )

    if not conditions:
        raise ValueError(
            f"unsupported logical expression for {target!r}: no interaction "
            f"condition found for {source!r}"
            + ("" if level is None else f" at level {level}")
        )

    condition = _join_or(conditions)
    return condition if len(conditions) == 1 else f"({condition})"


def _ginml_effect_levels(effects: Any) -> List[int]:
    if not isinstance(effects, str):
        return []

    levels = []

    for effect in effects.split():
        if ":" not in effect:
            continue

        level, _ = effect.split(":", maxsplit=1)
        levels.append(int(level))

    return levels


def _multivalued_range_condition(
    component: str,
    *,
    minimum: int,
    maximum: int,
    maxvalue: int,
) -> str:
    if minimum == 1 and maximum == 1 and maxvalue == 1:
        return component

    level_conditions = []

    for level in range(minimum, maximum + 1):
        if level == 0:
            level_conditions.append(f"~{_threshold_component(component, 1)}")
            continue

        terms = [
            _threshold_component(component, threshold)
            for threshold in range(1, level + 1)
        ]

        if level < maxvalue:
            terms.append(f"~{_threshold_component(component, level + 1)}")

        level_conditions.append(_join_and(terms))

    return f"({_join_or(level_conditions)})"


def _threshold_component(component: str, threshold: int) -> str:
    return f"{component}_b{threshold}"


def _join_or(expressions: Iterable[str]) -> str:
    expressions = [expression for expression in expressions if expression]

    if not expressions:
        return "0"

    if len(expressions) == 1:
        return expressions[0]

    return " | ".join(f"({expression})" for expression in expressions)


def _join_and(expressions: Iterable[str]) -> str:
    expressions = [expression for expression in expressions if expression != "1"]

    if "0" in expressions:
        return "0"

    if not expressions:
        return "1"

    return " & ".join(expressions)


def _select_main_ginml_file(ginml_files: List[str]) -> str:
    for name in sorted(ginml_files):
        if Path(name).name == "regulatoryGraph.ginml":
            return name

    return sorted(ginml_files)[0]


def _parse_zginml_companion_files(
    model: ExecutableModel,
    archive_files: Mapping[str, bytes],
    main_ginml: str,
) -> None:
    unknown_companion_files = []
    unparsed_initial_states = {}
    initial_state_metadata = []

    for name, content in archive_files.items():
        if name == main_ginml:
            continue

        basename = Path(name).name.lower()

        if basename == "initialstate":
            renamed_components = model.metadata.get("boolean_component_names", {})
            component_names = {
                node_id: renamed_components.get(node_id, node_id)
                for node_id in model.metadata.get("nodes", {})
            }
            parsed, unparsed, state_metadata = _parse_initial_states(
                content,
                model.metadata.get("nodes", {}),
                component_names,
            )
            model.initial_states.update(parsed)
            unparsed_initial_states.update(unparsed)
            initial_state_metadata.extend(state_metadata)
            continue

        if basename == "mutant":
            perturbations, perturbation_users = _parse_perturbations(content)
            model.perturbations.update(perturbations)

            if perturbation_users:
                model.metadata["perturbation_users"] = perturbation_users

            continue

        if basename == "avatar_parameters":
            model.metadata["avatar_parameters"] = _parse_avatar_parameters(content)
            continue

        if basename == "modelsimplifier":
            model.metadata["model_simplifier"] = _parse_model_simplifier(content)
            continue

        if basename == "reg2dyn_parameters":
            model.metadata["simulation_parameters"] = _parse_simulation_parameters(
                content
            )
            continue

        unknown_companion_files.append(name)

    if unparsed_initial_states:
        model.metadata["unparsed_initial_states"] = unparsed_initial_states

    if initial_state_metadata:
        model.metadata["initial_state_metadata"] = initial_state_metadata

    if unknown_companion_files:
        model.metadata["unknown_companion_files"] = sorted(unknown_companion_files)


def _parse_initial_states(
    content: bytes,
    nodes: Mapping[str, Any],
    component_names: Mapping[str, str],
) -> Tuple[
    Dict[str, Hypercube],
    Dict[str, Dict[str, Any]],
    List[Dict[str, Any]],
]:
    root = ElementTree.fromstring(content.strip())
    initial_states = {}
    unparsed = {}
    metadata = []
    section_counts: Dict[str, int] = {}

    for element in list(root):
        element_name = _local_name(element.tag)

        if element_name not in ["initialState", "input"]:
            continue

        section_counts[element_name] = section_counts.get(element_name, 0) + 1

        original_name = element.get("name")
        value = element.get("value")

        if value is None:
            continue

        name = _initial_state_name(
            original_name,
            section=element_name,
            index=section_counts[element_name],
            existing_names=initial_states,
        )

        try:
            initial_states[name] = Hypercube(
                _parse_state_assignment(value, nodes, component_names)
            )
            metadata.append(
                {
                    "name": name,
                    "original_name": original_name,
                    "section": element_name,
                }
            )
        except ValueError as error:
            unparsed[name] = {
                "section": element_name,
                "original_name": original_name,
                "value": value,
                "reason": str(error),
            }

    return initial_states, unparsed, metadata


def _initial_state_name(
    name: Optional[str],
    *,
    section: str,
    index: int,
    existing_names: Mapping[str, Hypercube],
) -> str:
    if name is None or name == "":
        name = f"{section}_{index}"

    if name not in existing_names:
        return name

    suffix = 2

    while f"{name}_{suffix}" in existing_names:
        suffix += 1

    return f"{name}_{suffix}"


def _parse_avatar_parameters(content: bytes) -> Dict[str, Any]:
    root = ElementTree.fromstring(content.strip())
    metadata: Dict[str, Any] = dict(root.attrib)

    if "nodeOrder" in metadata:
        metadata["nodeOrder"] = _parse_node_order(metadata["nodeOrder"])

    metadata["parameters"] = [
        _ginml_xml_metadata(element) for element in _children(root, "parameter")
    ]

    return metadata


def _parse_model_simplifier(content: bytes) -> Dict[str, Any]:
    root = ElementTree.fromstring(content.strip())
    simplifications: List[Dict[str, Any]] = []
    strip_outputs = []
    users = []

    for element in root.iter():
        element_name = _local_name(element.tag)

        if element_name == "stripOutput":
            strip_outputs.append(dict(element.attrib))
            continue

        if element_name == "user":
            users.append(dict(element.attrib))
            continue

        if element_name != "simplificationConfig":
            continue

        attributes: Dict[str, Any] = dict(element.attrib)

        if "removeList" in attributes:
            attributes["removeList"] = _parse_node_order(attributes["removeList"])

        simplifications.append(attributes)

    return {
        "attributes": dict(root.attrib),
        "simplifications": simplifications,
        "strip_outputs": strip_outputs,
        "users": users,
    }


def _parse_simulation_parameters(content: bytes) -> Dict[str, Any]:
    root = ElementTree.fromstring(content.strip())
    metadata: Dict[str, Any] = dict(root.attrib)

    if "nodeOrder" in metadata:
        metadata["nodeOrder"] = _parse_node_order(metadata["nodeOrder"])

    parameters: Dict[str, Dict[str, Any]] = {}
    parameter_sections: Dict[str, List[Dict[str, Any]]] = {}

    for element in _children(root, "parameter"):
        name = element.get("name")

        if name is None:
            continue

        parameters[name] = {
            "attributes": dict(element.attrib),
            "initstates": _row_names(element, "initstates"),
            "inputs": _row_names(element, "inputs"),
        }

        sections = [
            _ginml_xml_metadata(child)
            for child in list(element)
            if _local_name(child.tag) not in ["initstates", "inputs"]
        ]

        if sections:
            parameter_sections[name] = sections

    metadata["parameters"] = parameters
    metadata["parameter_sections"] = parameter_sections
    metadata["priority_class_lists"] = [
        {
            "attributes": dict(element.attrib),
            "classes": [
                {
                    **dict(class_element.attrib),
                    "content": _parse_node_order(class_element.get("content")),
                }
                for class_element in _children(element, "class")
            ],
        }
        for element in _children(root, "priorityClassList")
    ]

    return metadata


def _row_names(
    element: ElementTree.Element,
    section: str,
) -> List[str]:
    rows = []

    for section_element in _children(element, section):
        for row in _children(section_element, "row"):
            name = row.get("name")

            if name is not None:
                rows.append(name)

    return rows


def _parse_perturbations(content: bytes) -> Tuple[Dict[str, Any], List[Dict[str, str]]]:
    root = ElementTree.fromstring(content.strip())
    perturbations = {}
    users = []

    for element in root.iter():
        element_name = _local_name(element.tag)

        if element_name == "user":
            users.append(dict(element.attrib))
            continue

        if element_name != "mutant":
            continue

        name = element.get("name")

        if name is None:
            continue

        perturbations[name] = [
            {
                "type": _local_name(child.tag),
                "attributes": dict(child.attrib),
            }
            for child in list(element)
        ]

    return perturbations, users


def _parse_state_assignment(
    value: str,
    nodes: Mapping[str, Any],
    component_names: Mapping[str, str],
) -> Dict[str, int]:
    state = {}

    for item in value.split():
        if ";" not in item:
            raise ValueError(f"invalid state assignment item {item!r}")

        component, value_text = item.split(";", maxsplit=1)

        if component not in nodes:
            raise ValueError(f"unknown state component {component!r}")

        if not value_text.isdigit():
            raise ValueError(
                f"unsupported state value {value_text!r}: only Boolean "
                "or non-negative multi-valued levels can be represented as "
                "Hypercube values"
            )

        component_value = int(value_text)
        component_metadata = nodes.get(component)
        boolean_component = component_names.get(component, component)
        maxvalue = 1

        if isinstance(component_metadata, Mapping):
            maxvalue = int(component_metadata.get("maxvalue", 1))

        if component_value > maxvalue:
            raise ValueError(
                f"unsupported state value {value_text!r} for {component!r}: "
                f"expected a value between 0 and {maxvalue}"
            )

        if maxvalue == 1:
            state[boolean_component] = component_value
            continue

        for threshold in range(1, maxvalue + 1):
            state[_threshold_component(boolean_component, threshold)] = int(
                component_value >= threshold
            )

    return state


def _children(
    element: ElementTree.Element,
    name: str,
) -> List[ElementTree.Element]:
    return [child for child in list(element) if _local_name(child.tag) == name]


def _ginml_xml_metadata(element: ElementTree.Element) -> Dict[str, Any]:
    metadata: Dict[str, Any] = {
        "tag": _local_name(element.tag),
        "attributes": {
            _local_name(key): value for key, value in element.attrib.items()
        },
    }
    text = (element.text or "").strip()
    children = [_ginml_xml_metadata(child) for child in list(element)]

    if text:
        metadata["text"] = text

    if children:
        metadata["children"] = children

    return metadata


def _local_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", maxsplit=1)[-1]

    return tag


def _parse_node_order(nodeorder: Optional[str]) -> List[str]:
    if nodeorder is None:
        return []

    return [node for node in nodeorder.split() if node]


def _parse_int_attribute(
    element: ElementTree.Element,
    name: str,
    *,
    default: int,
) -> int:
    value = element.get(name)

    if value is None:
        return default

    return int(value)


def _parse_bool_attribute(
    element: ElementTree.Element,
    name: str,
    *,
    default: bool,
) -> bool:
    value = element.get(name)

    if value is None:
        return default

    return value.lower() == "true"


def _ginml_edge_signs(edge: Mapping[str, Any]) -> List[int]:
    sign = edge.get("sign")

    if sign in ["positive", "+", "1"]:
        return [1]

    if sign in ["negative", "-", "-1"]:
        return [-1]

    if sign == "dual":
        return [1, -1]

    effects = edge.get("attributes", {}).get("effects")

    if isinstance(effects, str):
        signs = []

        for effect in effects.split():
            if ":" not in effect:
                continue

            _, effect_sign = effect.split(":", maxsplit=1)
            signs.extend(_ginml_edge_signs({"sign": effect_sign}))

        return sorted(set(signs), reverse=True)

    return []


_TOKEN_PATTERN = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*(?::[0-9]+)?\b")
