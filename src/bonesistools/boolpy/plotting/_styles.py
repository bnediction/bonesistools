#!/usr/bin/env python

from typing import (
    Any,
    Dict,
    Mapping,
)


def ratio_edge_style(ratio: float) -> Dict[str, str]:
    """
    Return pydot edge attributes from an ensemble occurrence ratio.

    Parameters
    ----------
    ratio: float
        Fraction of Boolean networks in which the signed influence is observed.
        Expected to range between 0 and 1.

    Returns
    -------
    Dict[str, str]
        Mapping of pydot edge attributes.
    """

    if ratio < 0.25:
        style = "dotted"

    elif ratio < 0.5:
        style = "dashed"

    elif ratio < 0.75:
        style = "solid"

    else:
        style = "bold"

    return {
        "style": style,
        "penwidth": "2",
    }


def count_node_style(data: Mapping[str, Any]) -> Dict[str, str]:
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
    data: Mapping[str, Any]
        Node attributes containing the `function_count` value.

    Returns
    -------
    Dict[str, str]
        Pydot node attributes.
    """

    count = data["function_count"]

    if count == 1:
        fillcolor = "darkgoldenrod2"
        style = "rounded,filled,bold"

    elif count == 2:
        fillcolor = "lightgoldenrod1"
        style = "rounded,filled"

    elif count == 3:
        fillcolor = "cornsilk"
        style = "rounded,filled"

    elif count < 10:
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


def stability_node_style(data: Mapping[str, Any]) -> Dict[str, str]:
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
    data: Mapping[str, Any]
        Node attributes containing the `function_stability` value.

    Returns
    -------
    Dict[str, str]
        Pydot node attributes.
    """

    stability = data["function_stability"]

    if stability == 1:
        fillcolor = "darkgoldenrod2"
        style = "rounded,filled,bold"

    elif stability >= 0.75:
        fillcolor = "lightgoldenrod1"
        style = "rounded,filled"

    elif stability >= 0.5:
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
