#!/usr/bin/env python

from bonesistools.boolpy.plotting import (
    count_node_style,
    ratio_edge_style,
    stability_node_style,
)


def test_ratio_edge_style_thresholds():
    assert ratio_edge_style(0.0)["style"] == "dotted"
    assert ratio_edge_style(0.25)["style"] == "dashed"
    assert ratio_edge_style(0.5)["style"] == "solid"
    assert ratio_edge_style(0.75)["style"] == "bold"


def test_count_and_stability_node_style_thresholds():
    assert count_node_style({"function_count": 1}) == {
        "fillcolor": "darkgoldenrod2",
        "style": "rounded,filled,bold",
        "shape": "oval",
    }
    assert count_node_style({"function_count": 2})["fillcolor"] == "lightgoldenrod1"
    assert count_node_style({"function_count": 3})["fillcolor"] == "cornsilk"
    assert count_node_style({"function_count": 4})["style"] == "rounded,filled"
    assert count_node_style({"function_count": 10})["style"] == "rounded,filled,dotted"

    assert stability_node_style({"function_stability": 1}) == {
        "fillcolor": "darkgoldenrod2",
        "style": "rounded,filled,bold",
        "shape": "oval",
    }
    assert (
        stability_node_style({"function_stability": 0.75})["fillcolor"]
        == "lightgoldenrod1"
    )
    assert stability_node_style({"function_stability": 0.5})["fillcolor"] == "cornsilk"
    assert (
        stability_node_style({"function_stability": 0.49})["style"]
        == "rounded,filled,dotted"
    )
