#!/usr/bin/env python

import inspect
from collections.abc import Mapping as MappingInstance
from typing import Any, Callable, Dict, Mapping, Tuple


def frequency_edge_style(frequency: float) -> Dict[str, str]:
    """
    Return pydot edge attributes from an ensemble occurrence frequency.

    Parameters
    ----------
    frequency: float
        Fraction of Boolean networks in which the signed influence is observed.
        Expected to range between 0 and 1.

    Returns
    -------
    Dict[str, str]
        Mapping of pydot edge attributes.
    """

    if frequency < 0.25:
        style = "dotted"

    elif frequency < 0.5:
        style = "dashed"

    elif frequency < 0.75:
        style = "solid"

    else:
        style = "bold"

    return {
        "style": style,
        "penwidth": "2",
    }


def count_node_style(function_count: int) -> Dict[str, str]:
    """
    Return pydot node attributes from the number of distinct rule structures.

    `function_count` counts how many distinct Boolean rule structures were
    observed for a node across an ensemble. A value of 1 means all Boolean
    networks assign the same rule structure to the node. Larger values indicate
    more alternative inferred functions. The visual style therefore highlights
    low counts with stronger fill colors and high counts with lighter or dotted
    styles.

    Parameters
    ----------
    function_count: int
        Number of distinct Boolean rule structures observed for the node.

    Returns
    -------
    Dict[str, str]
        Pydot node attributes.
    """

    if function_count == 1:
        fillcolor = "darkgoldenrod2"
        style = "rounded,filled,bold"

    elif function_count == 2:
        fillcolor = "lightgoldenrod1"
        style = "rounded,filled"

    elif function_count == 3:
        fillcolor = "cornsilk"
        style = "rounded,filled"

    elif function_count < 10:
        fillcolor = "white"
        style = "rounded,filled"

    else:
        fillcolor = "white"
        style = "rounded,filled,dotted"

    return {
        "fillcolor": fillcolor,
        "style": style,
        "shape": "oval",
    }


def stability_node_style(function_stability: float) -> Dict[str, str]:
    """
    Return pydot node attributes from rule-structure stability.

    `function_stability` is the frequency of the most common Boolean rule
    structure for a node across an ensemble. A value of 1 means all Boolean
    networks assign the same rule structure to the node. Lower values indicate
    that the node has several competing inferred functions. The visual style
    therefore highlights stable nodes with stronger fill colors and unstable
    nodes with lighter or dotted styles.

    Parameters
    ----------
    function_stability: float
        Frequency of the most common Boolean rule structure for the node.

    Returns
    -------
    Dict[str, str]
        Pydot node attributes.
    """

    if function_stability == 1:
        fillcolor = "darkgoldenrod2"
        style = "rounded,filled,bold"

    elif function_stability >= 0.75:
        fillcolor = "lightgoldenrod1"
        style = "rounded,filled"

    elif function_stability >= 0.5:
        fillcolor = "cornsilk"
        style = "rounded,filled"

    else:
        fillcolor = "white"
        style = "rounded,filled,dotted"

    return {
        "fillcolor": fillcolor,
        "style": style,
        "shape": "oval",
    }


def _style_callable_parameters(
    style_callable: Callable[..., Mapping[str, Any]],
    style_name: str,
) -> Tuple[inspect.Parameter, ...]:
    """
    Return validated parameters for a style callable.
    """

    try:
        parameters = tuple(inspect.signature(style_callable).parameters.values())
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{style_name} callable must have an inspectable signature"
        ) from error

    unsupported = [
        parameter
        for parameter in parameters
        if parameter.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
    ]

    if unsupported:
        names = ", ".join(repr(parameter.name) for parameter in unsupported)
        raise TypeError(
            f"{style_name} callable must use explicitly named parameters; "
            f"unsupported parameter(s): {names}"
        )

    return parameters


def _evaluate_style_from_attributes(
    style_callable: Callable[..., Mapping[str, Any]],
    attributes: Mapping[str, Any],
    parameters: Tuple[Any, ...],
    *,
    style_name: str,
    element_name: str,
    element: Any,
) -> Mapping[str, Any]:
    """
    Evaluate a style callable from graph element attributes.
    """

    kwargs: Dict[str, Any] = {}
    missing_attributes = []

    for parameter in parameters:
        name = parameter.name

        if name in attributes:
            kwargs[name] = attributes[name]
            continue

        if parameter.default is inspect.Parameter.empty:
            missing_attributes.append(name)

    if missing_attributes:
        names = ", ".join(repr(name) for name in missing_attributes)
        raise ValueError(
            f"{style_name} callable requested missing {element_name} "
            f"attribute(s) for {element_name} {element!r}: {names}"
        )

    result = style_callable(**kwargs)

    if not isinstance(result, MappingInstance):
        raise TypeError(
            f"{style_name} callable must return a mapping of pydot attributes, "
            f"but returned {type(result)}"
        )

    return result
