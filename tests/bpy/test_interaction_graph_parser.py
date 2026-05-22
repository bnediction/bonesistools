#!/usr/bin/env python

import math

import networkx as nx
import pytest

import bonesistools as bt


def test_read_interaction_graph_keeps_edge_attributes(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text(
        "source,target,sign,confidence\n"
        "A,B,1,high\n"
        "B,C,-1,low\n"
        "C,A,,unknown\n",
    )

    graph = bt.bpy.ig.read_interaction_graph(infile)

    assert isinstance(graph, nx.MultiDiGraph)
    assert set(graph.nodes) == {"A", "B", "C"}
    assert graph["A"]["B"][0]["sign"] == 1.0
    assert graph["A"]["B"][0]["confidence"] == "high"
    assert graph["B"]["C"][0]["sign"] == -1.0
    assert math.isnan(graph["C"]["A"][0]["sign"])


def test_read_interaction_graph_supports_custom_separator(tmp_path):
    infile = tmp_path / "graph.tsv"
    infile.write_text(
        "source\ttarget\tsign\n" "A\tB\t1\n",
    )

    graph = bt.bpy.ig.read_interaction_graph(infile, sep="\t")

    assert list(graph.edges(data=True)) == [("A", "B", {"sign": 1.0})]


def test_read_interaction_graph_validates_file_columns_and_signs(tmp_path):
    missing = tmp_path / "missing.csv"

    with pytest.raises(FileNotFoundError):
        bt.bpy.ig.read_interaction_graph(missing)

    missing_column = tmp_path / "missing_column.csv"
    missing_column.write_text("source,target\nA,B\n")

    with pytest.raises(ValueError, match="missing columns"):
        bt.bpy.ig.read_interaction_graph(missing_column)

    invalid_sign = tmp_path / "invalid_sign.csv"
    invalid_sign.write_text("source,target,sign\nA,B,0\n")

    with pytest.raises(ValueError, match="unsupported sign values"):
        bt.bpy.ig.read_interaction_graph(invalid_sign)


def test_read_interaction_graph_rejects_invalid_genesyn(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        bt.bpy.ig.read_interaction_graph(infile, genesyn=object())
