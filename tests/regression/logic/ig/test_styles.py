#!/usr/bin/env python

import bonesistools as bt
from bonesistools.logic.influence_graph._styles import (
    count_node_style,
    frequency_edge_style,
    stability_node_style,
)


def test_frequency_edge_style_thresholds():
    assert frequency_edge_style(0.0)["style"] == "dotted"
    assert frequency_edge_style(0.25)["style"] == "dashed"
    assert frequency_edge_style(0.5)["style"] == "solid"
    assert frequency_edge_style(0.75)["style"] == "bold"


def test_count_and_stability_node_style_thresholds():
    assert count_node_style(1) == {
        "fillcolor": "darkgoldenrod2",
        "style": "rounded,filled,bold",
        "shape": "oval",
    }
    assert count_node_style(2)["fillcolor"] == "lightgoldenrod1"
    assert count_node_style(3)["fillcolor"] == "cornsilk"
    assert count_node_style(4)["style"] == "rounded,filled"
    assert count_node_style(10)["style"] == "rounded,filled,dotted"

    assert stability_node_style(1) == {
        "fillcolor": "darkgoldenrod2",
        "style": "rounded,filled,bold",
        "shape": "oval",
    }
    assert stability_node_style(0.75)["fillcolor"] == "lightgoldenrod1"
    assert stability_node_style(0.5)["fillcolor"] == "cornsilk"
    assert stability_node_style(0.49)["style"] == "rounded,filled,dotted"


def test_style_helpers_are_private_influence_graph_helpers():
    for name in ["count_node_style", "frequency_edge_style", "stability_node_style"]:
        assert name not in dir(bt.logic.ig)
        assert not hasattr(bt.logic.ig, name)

    assert "frequency_edge_style" not in dir(bt.logic.bn)
    assert not hasattr(bt.logic.bn, "frequency_edge_style")
