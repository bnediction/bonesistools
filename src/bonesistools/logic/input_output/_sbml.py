#!/usr/bin/env python

"""
SBML Level 3 Qual reader.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union, cast
from xml.etree import ElementTree

from ..boolean_algebra import Hypercube
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph
from ._booleanization import (
    _boolean_component,
    _booleanize_logical_model,
    _encode_state,
    _join_and,
    _join_or,
    _level_range_rule,
    _LogicalModel,
    _negate_rule,
    _normalize_component_names,
)
from ._executable_model import ExecutableModel

_UNSUPPORTED_SBML_DOCUMENT = (
    "unsupported SBML document "
    "(expected an SBML Level 3 document using the Qual package)"
)


def read_sbml(file: Union[str, Path]) -> ExecutableModel:
    """
    Read an SBML Level 3 Qual logical model.

    SBML documents without the Level 3 Qual package are rejected explicitly.
    Qualitative species and transitions are converted to a Boolean network,
    with multi-valued species represented by monotone threshold components.

    Examples
    --------
    >>> model = read_sbml("model.sbml")
    >>> network = model.get("boolean_network")

    Parameters
    ----------
    file: str or Path
        SBML Level 3 Qual file.

    Returns
    -------
    ExecutableModel
        Parsed executable logical model.

    Raises
    ------
    ValueError
        If the document is not SBML Level 3 Qual or uses an unsupported Qual
        or MathML construction.
    """

    file = Path(file)
    root = ElementTree.fromstring(file.read_bytes())
    qual_namespace = _validate_sbml_qual_document(root)
    model_element = _first_child(root, "model")

    if model_element is None:
        raise ValueError(_UNSUPPORTED_SBML_DOCUMENT)

    species = _parse_qualitative_species(model_element, qual_namespace)
    component_names = _normalize_component_names(species)
    transitions = _parse_qual_transitions(
        model_element,
        qual_namespace,
        species,
        component_names,
    )
    logical_model = _logical_model_from_sbml(
        species,
        transitions,
        component_names,
    )
    boolean_network = _booleanize_logical_model(logical_model)
    layout = _parse_sbml_layout(model_element)
    influence_graph = _build_sbml_influence_graph(
        boolean_network,
        species,
        transitions,
        component_names,
        layout,
        model_element,
    )
    initial_states = _sbml_initial_states(species, component_names)
    renamed_components = {
        source: target for source, target in component_names.items() if source != target
    }
    metadata = {
        "format": "sbml",
        "source": str(file),
        "sbml_level": 3,
        "sbml_version": int(root.get("version", "1")),
        "qual_namespace": qual_namespace,
        "model": _local_attributes(model_element),
        "qualitative_species": species,
        "transitions": transitions,
        "layout": layout,
        "annotations": _sbml_annotations(model_element),
    }

    if renamed_components:
        metadata["boolean_component_names"] = renamed_components

    return ExecutableModel(
        boolean_network=boolean_network,
        influence_graph=influence_graph,
        initial_states=initial_states,
        metadata=metadata,
    )


def _validate_sbml_qual_document(root: ElementTree.Element) -> str:
    """Validate the SBML level and locate its Qual namespace."""

    if (
        _local_name(root.tag) != "sbml"
        or root.get("level") != "3"
        or "/sbml/level3/" not in _namespace(root.tag)
    ):
        raise ValueError(_UNSUPPORTED_SBML_DOCUMENT)

    qual_namespaces = {
        _namespace(element.tag)
        for element in root.iter()
        if _local_name(element.tag)
        in {
            "listOfQualitativeSpecies",
            "qualitativeSpecies",
            "listOfTransitions",
            "transition",
        }
        and "/qual/" in _namespace(element.tag)
    }

    if not qual_namespaces:
        raise ValueError(_UNSUPPORTED_SBML_DOCUMENT)

    return sorted(qual_namespaces)[0]


def _parse_qualitative_species(
    model_element: ElementTree.Element,
    qual_namespace: str,
) -> Dict[str, Dict[str, Any]]:
    """Parse Qual species and their activity domains."""

    species = {}

    for element in model_element.iter():
        if element.tag != f"{{{qual_namespace}}}qualitativeSpecies":
            continue

        attributes = _local_attributes(element)
        component = attributes.get("id")

        if not component:
            raise ValueError("unsupported SBML Qual species without an identifier")

        max_level = int(attributes.get("maxLevel", 1))

        if max_level < 1:
            raise ValueError(
                f"unsupported maximum level {max_level!r} for SBML Qual "
                f"species {component!r}"
            )

        initial_level = attributes.get("initialLevel")

        if initial_level is not None:
            initial_level = int(initial_level)

            if initial_level < 0 or initial_level > max_level:
                raise ValueError(
                    f"unsupported initial level {initial_level!r} for SBML "
                    f"Qual species {component!r}"
                )

        species[component] = {
            "id": component,
            "name": attributes.get("name", component),
            "max_level": max_level,
            "constant": _as_boolean(attributes.get("constant", "false")),
            "initial_level": initial_level,
            "attributes": attributes,
            "annotations": _element_annotations(element),
        }

    if not species:
        raise ValueError(_UNSUPPORTED_SBML_DOCUMENT)

    return species


def _parse_qual_transitions(
    model_element: ElementTree.Element,
    qual_namespace: str,
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> List[Dict[str, Any]]:
    """Parse Qual transitions and convert their MathML conditions."""

    transitions = []

    for element in model_element.iter():
        if element.tag != f"{{{qual_namespace}}}transition":
            continue

        attributes = _local_attributes(element)
        inputs = [
            _local_attributes(child)
            for child in element.iter()
            if child.tag == f"{{{qual_namespace}}}input"
        ]
        outputs = [
            _local_attributes(child)
            for child in element.iter()
            if child.tag == f"{{{qual_namespace}}}output"
        ]

        if not outputs:
            raise ValueError(
                f"unsupported SBML Qual transition {attributes.get('id')!r} "
                "without an output"
            )

        default_terms = [
            child
            for child in element.iter()
            if child.tag == f"{{{qual_namespace}}}defaultTerm"
        ]

        if len(default_terms) > 1:
            raise ValueError(
                f"unsupported SBML Qual transition {attributes.get('id')!r} "
                "with multiple default terms"
            )

        default_level = (
            int(_local_attributes(default_terms[0]).get("resultLevel", 0))
            if default_terms
            else 0
        )
        function_terms = []

        for term in element.iter():
            if term.tag != f"{{{qual_namespace}}}functionTerm":
                continue

            term_attributes = _local_attributes(term)

            if "resultLevel" not in term_attributes:
                raise ValueError(
                    f"unsupported SBML Qual function term in transition "
                    f"{attributes.get('id')!r} without a result level"
                )

            math = next(
                (child for child in list(term) if _local_name(child.tag) == "math"),
                None,
            )

            if math is None:
                raise ValueError(
                    f"unsupported SBML Qual function term in transition "
                    f"{attributes.get('id')!r} without MathML"
                )

            function_terms.append(
                {
                    "result_level": int(term_attributes["resultLevel"]),
                    "condition": _parse_mathml_condition(
                        math,
                        species,
                        component_names,
                    ),
                    "attributes": term_attributes,
                }
            )

        transitions.append(
            {
                "id": attributes.get("id"),
                "attributes": attributes,
                "inputs": inputs,
                "outputs": outputs,
                "default_level": default_level,
                "function_terms": function_terms,
                "annotations": _element_annotations(element),
            }
        )

    return transitions


def _logical_model_from_sbml(
    species: Mapping[str, Mapping[str, Any]],
    transitions: Iterable[Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> _LogicalModel:
    """Build the shared logical representation from Qual transitions."""

    transitions_by_output = {}

    for transition in transitions:
        for output in transition["outputs"]:
            component = output.get("qualitativeSpecies")

            if component not in species:
                raise ValueError(
                    f"unsupported SBML Qual transition output {component!r}: "
                    "unknown qualitative species"
                )

            if output.get("transitionEffect", "assignmentLevel") not in [
                "assignmentLevel",
            ]:
                raise ValueError(
                    f"unsupported SBML Qual transition effect "
                    f"{output.get('transitionEffect')!r} for {component!r}"
                )

            if component in transitions_by_output:
                raise ValueError(
                    f"unsupported multiple SBML Qual transitions targeting "
                    f"{component!r}"
                )

            transitions_by_output[component] = transition

    max_levels = {
        component: int(data["max_level"]) for component, data in species.items()
    }
    input_components = frozenset(
        component
        for component, data in species.items()
        if bool(data["constant"]) or component not in transitions_by_output
    )
    threshold_rules = {}

    for component, transition in transitions_by_output.items():
        max_level = max_levels[component]
        default_level = int(transition["default_level"])
        terms = transition["function_terms"]

        if default_level < 0 or default_level > max_level:
            raise ValueError(
                f"unsupported default level {default_level!r} for SBML Qual "
                f"species {component!r}"
            )

        for term in terms:
            result_level = int(term["result_level"])

            if result_level < 0 or result_level > max_level:
                raise ValueError(
                    f"unsupported result level {result_level!r} for SBML Qual "
                    f"species {component!r}"
                )

        all_term_conditions = _join_or(term["condition"] for term in terms)
        default_condition = _negate_rule(all_term_conditions)
        threshold_rules[component] = {}

        for threshold in range(1, max_level + 1):
            conditions = [
                term["condition"]
                for term in terms
                if int(term["result_level"]) >= threshold
            ]

            if default_level >= threshold:
                conditions.append(default_condition)

            threshold_rules[component][threshold] = _join_or(conditions)

    return _LogicalModel(
        max_levels=max_levels,
        input_components=input_components,
        threshold_rules=threshold_rules,
        component_names=component_names,
    )


def _parse_mathml_condition(
    math: ElementTree.Element,
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    """Convert a supported MathML Boolean condition to a Boolean rule."""

    children = list(math)

    if len(children) != 1:
        raise ValueError("unsupported SBML Qual MathML condition")

    return _parse_mathml_expression(children[0], species, component_names)


def _parse_mathml_expression(
    element: ElementTree.Element,
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    """Recursively convert one MathML expression."""

    name = _local_name(element.tag)

    if name in ["true", "false"]:
        return "1" if name == "true" else "0"

    if name != "apply":
        raise ValueError(f"unsupported SBML Qual MathML element {name!r}")

    children = list(element)

    if not children:
        raise ValueError("unsupported empty SBML Qual MathML application")

    operator = _local_name(children[0].tag)
    operands = children[1:]

    if operator in ["and", "or"]:
        rules = [
            _parse_mathml_expression(operand, species, component_names)
            for operand in operands
        ]
        return _join_and(rules) if operator == "and" else _join_or(rules)

    if operator == "not":
        if len(operands) != 1:
            raise ValueError("unsupported SBML Qual MathML 'not' arity")

        return _negate_rule(
            _parse_mathml_expression(operands[0], species, component_names)
        )

    if operator == "xor":
        if len(operands) < 2:
            raise ValueError("unsupported SBML Qual MathML 'xor' arity")

        result = _parse_mathml_expression(operands[0], species, component_names)

        for operand in operands[1:]:
            right = _parse_mathml_expression(operand, species, component_names)
            result = _join_or(
                [
                    _join_and([result, _negate_rule(right)]),
                    _join_and([_negate_rule(result), right]),
                ]
            )

        return result

    if operator in ["eq", "neq", "lt", "leq", "gt", "geq"]:
        if len(operands) != 2:
            raise ValueError(
                f"unsupported SBML Qual MathML {operator!r} comparison arity"
            )

        return _mathml_comparison_rule(
            operator,
            _parse_mathml_operand(operands[0], species),
            _parse_mathml_operand(operands[1], species),
            species,
            component_names,
        )

    raise ValueError(f"unsupported SBML Qual MathML operator {operator!r}")


def _parse_mathml_operand(
    element: ElementTree.Element,
    species: Mapping[str, Mapping[str, Any]],
) -> Tuple[str, Union[str, int]]:
    """Parse a MathML comparison operand."""

    name = _local_name(element.tag)
    text = (element.text or "").strip()

    if name == "ci":
        if text not in species:
            raise ValueError(f"unsupported SBML Qual MathML species reference {text!r}")

        return "component", text

    if name == "cn":
        try:
            return "number", int(text)
        except ValueError:
            raise ValueError(
                f"unsupported non-integer SBML Qual MathML constant {text!r}"
            )

    raise ValueError(f"unsupported SBML Qual MathML operand {name!r}")


def _mathml_comparison_rule(
    operator: str,
    left: Tuple[str, Union[str, int]],
    right: Tuple[str, Union[str, int]],
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    """Expand a MathML comparison over discrete component levels."""

    left_values = _mathml_operand_values(left, species)
    right_values = _mathml_operand_values(right, species)
    conditions = []

    for left_value in left_values:
        for right_value in right_values:
            if not _compare_levels(left_value, right_value, operator):
                continue

            terms = []

            if left[0] == "component":
                terms.append(
                    _exact_level_rule(
                        cast(str, left[1]),
                        left_value,
                        species,
                        component_names,
                    )
                )

            if right[0] == "component":
                terms.append(
                    _exact_level_rule(
                        cast(str, right[1]),
                        right_value,
                        species,
                        component_names,
                    )
                )

            conditions.append(_join_and(terms))

    return _join_or(conditions)


def _mathml_operand_values(
    operand: Tuple[str, Union[str, int]],
    species: Mapping[str, Mapping[str, Any]],
) -> range:
    """Return all values represented by a MathML operand."""

    if operand[0] == "number":
        value = cast(int, operand[1])
        return range(value, value + 1)

    component = cast(str, operand[1])
    return range(0, int(species[component]["max_level"]) + 1)


def _compare_levels(left: int, right: int, operator: str) -> bool:
    """Evaluate one discrete MathML comparison."""

    if operator == "eq":
        return left == right
    if operator == "neq":
        return left != right
    if operator == "lt":
        return left < right
    if operator == "leq":
        return left <= right
    if operator == "gt":
        return left > right

    return left >= right


def _exact_level_rule(
    source_component: str,
    level: int,
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> str:
    """Encode one exact qualitative-species level."""

    return _level_range_rule(
        component_names[source_component],
        minimum=level,
        maximum=level,
        max_level=int(species[source_component]["max_level"]),
    )


def _sbml_initial_states(
    species: Mapping[str, Mapping[str, Any]],
    component_names: Mapping[str, str],
) -> Dict[str, Hypercube]:
    """Convert explicit initial levels to a named Boolean hypercube."""

    state = {
        component: int(data["initial_level"])
        for component, data in species.items()
        if data["initial_level"] is not None
    }

    if not state:
        return {}

    encoded = _encode_state(
        state,
        max_levels={
            component: int(data["max_level"]) for component, data in species.items()
        },
        component_names=component_names,
    )
    return {"default": Hypercube(encoded)}


def _build_sbml_influence_graph(
    boolean_network: BooleanNetwork,
    species: Mapping[str, Mapping[str, Any]],
    transitions: Iterable[Mapping[str, Any]],
    component_names: Mapping[str, str],
    layout: Mapping[str, Mapping[str, Any]],
    model_element: ElementTree.Element,
) -> InfluenceGraph:
    """Build an attributed influence graph aligned with the Boolean network."""

    transitions = list(transitions)
    graph = InfluenceGraph()
    graph.graph.update(
        {
            f"sbml_{_attribute_name(key)}": value
            for key, value in _local_attributes(model_element).items()
        }
    )
    component_origins = {}
    transition_targets = {
        output.get("qualitativeSpecies")
        for transition in transitions
        for output in transition["outputs"]
    }

    for source_component, data in species.items():
        component = component_names[source_component]
        max_level = int(data["max_level"])

        for threshold in range(1, max_level + 1):
            boolean_component = _boolean_component(
                component,
                threshold=threshold,
                max_level=max_level,
            )
            component_origins[boolean_component] = source_component
            attributes = {
                "name": (data["name"] if max_level == 1 else boolean_component),
                "maxvalue": 1,
                "input": bool(data["constant"])
                or source_component not in transition_targets,
                **{
                    f"sbml_{_attribute_name(key)}": value
                    for key, value in data["attributes"].items()
                    if key not in ["id", "name", "maxLevel", "constant"]
                },
            }

            if boolean_component != source_component or max_level > 1:
                attributes["sbml_component"] = source_component

            if max_level > 1:
                attributes["sbml_name"] = data["name"]
                attributes["sbml_max_level"] = max_level
                attributes["sbml_threshold"] = threshold

            attributes.update(
                {
                    f"sbml_{_attribute_name(key)}": value
                    for key, value in layout.get(source_component, {}).items()
                }
            )
            graph.add_node(boolean_component, **attributes)

    inputs_by_pair = _sbml_inputs_by_pair(transitions)

    for source, target, sign in sorted(boolean_network.influences()):
        source_component = component_origins[source]
        target_component = component_origins[target]
        candidates = inputs_by_pair.get((source_component, target_component), [])
        attributes = {}

        if candidates:
            attributes["sbml_inputs"] = candidates
        else:
            attributes["sbml_generated"] = True

        graph.add_edge(source, target, sign=cast(Any, sign), **attributes)

    return graph


def _sbml_inputs_by_pair(
    transitions: Iterable[Mapping[str, Any]],
) -> Dict[Tuple[str, str], List[Dict[str, Any]]]:
    """Index source-format transition inputs by regulator and target."""

    indexed: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}

    for transition in transitions:
        targets = [output.get("qualitativeSpecies") for output in transition["outputs"]]

        for input_data in transition["inputs"]:
            source = input_data.get("qualitativeSpecies")

            for target in targets:
                if source is None or target is None:
                    continue

                indexed.setdefault((source, target), []).append(dict(input_data))

    return indexed


def _parse_sbml_layout(
    model_element: ElementTree.Element,
) -> Dict[str, Dict[str, Any]]:
    """Extract species positions and dimensions from the Layout package."""

    layout = {}

    for glyph in model_element.iter():
        if _local_name(glyph.tag) not in ["generalGlyph", "speciesGlyph"]:
            continue

        attributes = _local_attributes(glyph)
        reference = attributes.get("reference", attributes.get("species"))

        if reference is None:
            continue

        visual = {}

        for child in glyph.iter():
            child_name = _local_name(child.tag)
            child_attributes = _local_attributes(child)

            if child_name == "position":
                visual.update(
                    {
                        key: value
                        for key, value in child_attributes.items()
                        if key in ["x", "y", "z"]
                    }
                )
            elif child_name == "dimensions":
                visual.update(
                    {
                        key: value
                        for key, value in child_attributes.items()
                        if key in ["width", "height", "depth"]
                    }
                )

        layout[reference] = visual

    return layout


def _sbml_annotations(model_element: ElementTree.Element) -> List[str]:
    """Preserve model-level notes and annotations as XML strings."""

    return [
        ElementTree.tostring(child, encoding="unicode")
        for child in list(model_element)
        if _local_name(child.tag) in ["notes", "annotation"]
    ]


def _element_annotations(element: ElementTree.Element) -> List[str]:
    """Preserve notes and annotations attached to one SBML element."""

    return [
        ElementTree.tostring(child, encoding="unicode")
        for child in list(element)
        if _local_name(child.tag) in ["notes", "annotation"]
    ]


def _first_child(
    element: ElementTree.Element,
    name: str,
) -> Optional[ElementTree.Element]:
    """Return the first direct child with the requested local name."""

    return next(
        (child for child in list(element) if _local_name(child.tag) == name),
        None,
    )


def _namespace(tag: str) -> str:
    """Return an expanded XML tag namespace."""

    if tag.startswith("{") and "}" in tag:
        return tag[1:].split("}", maxsplit=1)[0]

    return ""


def _local_name(tag: str) -> str:
    """Return an XML tag or attribute local name."""

    return tag.split("}", maxsplit=1)[-1]


def _local_attributes(element: ElementTree.Element) -> Dict[str, str]:
    """Return XML attributes indexed by local name."""

    return {_local_name(key): value for key, value in element.attrib.items()}


def _attribute_name(name: str) -> str:
    """Normalize an SBML attribute name for graph metadata."""

    return "".join(character if character.isalnum() else "_" for character in name)


def _as_boolean(value: Any) -> bool:
    """Parse an SBML Boolean attribute."""

    if isinstance(value, bool):
        return value

    normalized = str(value).strip().lower()

    if normalized in ["true", "1"]:
        return True
    if normalized in ["false", "0"]:
        return False

    raise ValueError(f"unsupported SBML Boolean value {value!r}")
