#!/usr/bin/env python

"""
GINML and ZGINML readers.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union, cast
from xml.etree import ElementTree

from ..boolean_algebra import Hypercube
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph
from ._executable_model import ExecutableModel


def read_ginml(path: Union[str, Path]) -> ExecutableModel:
    """
    Read a GINML logical model.

    GINML import is conservative. Signed interactions are converted to an
    `InfluenceGraph` when possible. Boolean rules are converted to a
    `BooleanNetwork` only when the logical parameters are safely representable.
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

    nodes = _parse_ginml_nodes(graph_element)
    edges = _parse_ginml_edges(graph_element)
    metadata = _ginml_metadata(
        source=source,
        graph_element=graph_element,
        nodes=nodes,
        edges=edges,
    )

    influence_graph = _build_influence_graph(nodes, edges, metadata)
    boolean_network = _build_boolean_network(nodes, metadata)

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

            elif child_name not in ["annotation", "nodevisualsetting"]:
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
            "unsupported_children": unsupported_children,
        }

    return nodes


def _parse_ginml_edges(graph_element: ElementTree.Element) -> List[Dict[str, Any]]:
    edges = []

    for element in _children(graph_element, "edge"):
        source = element.get("from")
        target = element.get("to")

        if source is None or target is None:
            raise ValueError("invalid GINML file: edge without 'from' or 'to'")

        edges.append(
            {
                "id": element.get("id"),
                "source": source,
                "target": target,
                "sign": element.get("sign"),
                "attributes": dict(element.attrib),
                "children": [
                    ElementTree.tostring(child, encoding="unicode")
                    for child in list(element)
                ],
            }
        )

    return edges


def _ginml_metadata(
    *,
    source: str,
    graph_element: ElementTree.Element,
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
        "node_order": _parse_node_order(graph_element.get("nodeorder")),
        "nodes": {
            node_id: {
                "name": node["name"],
                "maxvalue": node["maxvalue"],
                "input": node["input"],
                "attributes": node["attributes"],
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
        "unsupported_sections": [],
    }


def _build_influence_graph(
    nodes: Mapping[str, Mapping[str, Any]],
    edges: List[Dict[str, Any]],
    metadata: Dict[str, Any],
) -> Optional[InfluenceGraph]:
    graph = InfluenceGraph()

    for node_id, node in nodes.items():
        graph.add_node(
            node_id,
            name=node["name"],
            maxvalue=node["maxvalue"],
            input=node["input"],
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


def _build_boolean_network(
    nodes: Mapping[str, Mapping[str, Any]],
    metadata: Dict[str, Any],
) -> Optional[BooleanNetwork]:
    try:
        rules = _boolean_rules_from_ginml_nodes(nodes)
        return BooleanNetwork(rules)

    except ValueError as error:
        metadata["boolean_network_reason"] = str(error)
        return None


def _boolean_rules_from_ginml_nodes(
    nodes: Mapping[str, Mapping[str, Any]],
) -> Dict[str, str]:
    rules = {}

    for node_id, node in nodes.items():
        maxvalue = int(node["maxvalue"])

        if maxvalue == 1:
            rules[node_id] = _boolean_rule_for_native_node(node_id, node, nodes)
            continue

        for threshold in range(1, maxvalue + 1):
            rules[_threshold_component(node_id, threshold)] = (
                _boolean_rule_for_threshold_node(node_id, threshold, node, nodes)
            )

    return rules


def _boolean_rule_for_native_node(
    node_id: str,
    node: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
) -> str:
    if bool(node["input"]):
        return node_id

    constant = _constant_rule_from_parameters(node)

    if constant is not None:
        if constant not in ["0", "1"]:
            raise ValueError(
                "unsupported logical parameters: native Boolean components "
                "must have constant values 0 or 1"
            )

        return constant

    expressions = _expressions_for_value(node, 1)

    if not expressions:
        return "0"

    converted = [
        _convert_ginml_expression(expression, nodes)
        for expression in expressions
    ]

    return _join_or(converted)


def _boolean_rule_for_threshold_node(
    node_id: str,
    threshold: int,
    node: Mapping[str, Any],
    nodes: Mapping[str, Mapping[str, Any]],
) -> str:
    threshold_component = _threshold_component(node_id, threshold)

    if bool(node["input"]):
        if threshold == 1:
            return threshold_component

        return (
            f"{_threshold_component(node_id, threshold - 1)} "
            f"& {threshold_component}"
        )

    constant = _constant_rule_from_parameters(node)

    if constant is not None:
        return "1" if int(constant) >= threshold else "0"

    expressions = []

    for value in node["values"]:
        if int(value["val"]) >= threshold:
            expressions.extend(value["expressions"])

    if not expressions:
        return "0"

    converted = [
        _convert_ginml_expression(expression, nodes)
        for expression in expressions
    ]

    return _join_or(converted)


def _constant_rule_from_parameters(node: Mapping[str, Any]) -> Optional[str]:
    parameters = list(node["parameters"])

    if not parameters:
        return None

    if len(parameters) != 1 or set(parameters[0]) != {"val"}:
        raise ValueError(
            "unsupported logical parameters: only constant parameter values "
            "can be converted safely"
        )

    return str(int(parameters[0]["val"]))


def _expressions_for_value(node: Mapping[str, Any], value: int) -> List[str]:
    expressions = []

    for value_data in node["values"]:
        if int(value_data["val"]) == value:
            expressions.extend(value_data["expressions"])

    return expressions


def _convert_ginml_expression(
    expression: str,
    nodes: Mapping[str, Mapping[str, Any]],
) -> str:
    def replace_token(match: re.Match[str]) -> str:
        token = match.group(0)

        if token in ["true", "True"]:
            return "1"

        if token in ["false", "False"]:
            return "0"

        if ":" in token:
            component, threshold_text = token.split(":", maxsplit=1)
            threshold = int(threshold_text)

            if component not in nodes:
                raise ValueError(
                    f"unsupported logical parameters: unknown component "
                    f"{component!r}"
                )

            maxvalue = int(nodes[component]["maxvalue"])

            if threshold < 1 or threshold > maxvalue:
                raise ValueError(
                    f"unsupported logical parameters: threshold {token!r} is "
                    f"outside component range 1..{maxvalue}"
                )

            return _threshold_condition(component, threshold, maxvalue)

        if token not in nodes:
            raise ValueError(
                f"unsupported logical parameters: unknown component {token!r}"
            )

        maxvalue = int(nodes[token]["maxvalue"])

        if maxvalue > 1:
            return _threshold_component(token, 1)

        return token

    converted = _TOKEN_PATTERN.sub(replace_token, expression)
    return converted.replace("!", "~")


def _threshold_condition(component: str, threshold: int, maxvalue: int) -> str:
    terms = [
        _threshold_component(component, level)
        for level in range(1, threshold + 1)
    ]

    if threshold == 1 and maxvalue == 1:
        return component

    return " & ".join(terms)


def _threshold_component(component: str, threshold: int) -> str:
    return f"{component}_b{threshold}"


def _join_or(expressions: Iterable[str]) -> str:
    expressions = [expression for expression in expressions if expression]

    if not expressions:
        return "0"

    if len(expressions) == 1:
        return expressions[0]

    return " | ".join(f"({expression})" for expression in expressions)


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
            parsed, unparsed, state_metadata = _parse_initial_states(
                content,
                model.metadata.get("nodes", {}),
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
            initial_states[name] = Hypercube(_parse_state_assignment(value, nodes))
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

    return metadata


def _parse_model_simplifier(content: bytes) -> Dict[str, Any]:
    root = ElementTree.fromstring(content.strip())
    simplifications: List[Dict[str, Any]] = []

    for element in root.iter():
        if _local_name(element.tag) != "simplificationConfig":
            continue

        attributes: Dict[str, Any] = dict(element.attrib)

        if "removeList" in attributes:
            attributes["removeList"] = _parse_node_order(attributes["removeList"])

        simplifications.append(attributes)

    return {"simplifications": simplifications}


def _parse_simulation_parameters(content: bytes) -> Dict[str, Any]:
    root = ElementTree.fromstring(content.strip())
    metadata: Dict[str, Any] = dict(root.attrib)

    if "nodeOrder" in metadata:
        metadata["nodeOrder"] = _parse_node_order(metadata["nodeOrder"])

    parameters: Dict[str, Dict[str, Any]] = {}

    for element in _children(root, "parameter"):
        name = element.get("name")

        if name is None:
            continue

        parameters[name] = {
            "attributes": dict(element.attrib),
            "initstates": _row_names(element, "initstates"),
            "inputs": _row_names(element, "inputs"),
        }

    metadata["parameters"] = parameters

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
) -> Dict[str, int]:
    state = {}

    for item in value.split():
        if ";" not in item:
            raise ValueError(f"invalid state assignment item {item!r}")

        component, value_text = item.split(";", maxsplit=1)

        if not value_text.isdigit():
            raise ValueError(
                f"unsupported state value {value_text!r}: only Boolean "
                "or non-negative multi-valued levels can be represented as "
                "Hypercube values"
            )

        component_value = int(value_text)
        component_metadata = nodes.get(component)
        maxvalue = 1

        if isinstance(component_metadata, Mapping):
            maxvalue = int(component_metadata.get("maxvalue", 1))

        if component_value > maxvalue:
            raise ValueError(
                f"unsupported state value {value_text!r} for {component!r}: "
                f"expected a value between 0 and {maxvalue}"
            )

        if maxvalue == 1:
            state[component] = component_value
            continue

        for threshold in range(1, maxvalue + 1):
            state[_threshold_component(component, threshold)] = int(
                component_value >= threshold
            )

    return state


def _children(
    element: ElementTree.Element,
    name: str,
) -> List[ElementTree.Element]:
    return [child for child in list(element) if _local_name(child.tag) == name]


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
