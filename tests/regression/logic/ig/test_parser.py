#!/usr/bin/env python

import math
from typing import Any, cast

import networkx as nx
import pytest

import bonesistools as bt
from bonesistools.logic.input_output import _influence_graph


def _edge_data(graph, source, target, key=0):
    return cast(Any, graph[source][target])[key]


def test_read_influence_graph_keeps_edge_attributes(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text(
        "source,target,sign,confidence\n"
        "A,B,1,high\n"
        "B,C,-1,low\n"
        "C,A,,unknown\n",
    )

    graph = bt.logic.io.read_influence_graph(infile)

    assert isinstance(graph, nx.MultiDiGraph)
    assert set(graph.nodes) == {"A", "B", "C"}
    assert _edge_data(graph, "A", "B")["sign"] == 1.0
    assert _edge_data(graph, "A", "B")["confidence"] == "high"
    assert _edge_data(graph, "B", "C")["sign"] == -1.0
    assert math.isnan(_edge_data(graph, "C", "A")["sign"])


def test_read_influence_graph_supports_custom_separator(tmp_path):
    infile = tmp_path / "graph.tsv"
    infile.write_text(
        "source\ttarget\tsign\n" "A\tB\t1\n",
    )

    graph = bt.logic.io.read_influence_graph(infile, sep="\t")

    assert list(graph.edges(data=True)) == [("A", "B", {"sign": 1.0})]


def test_read_influence_graph_options_are_keyword_only(tmp_path):
    file = tmp_path / "graph.csv"
    file.write_text("source,target,sign\nA,B,1\n")

    with pytest.raises(TypeError):
        cast(Any, bt.logic.io.read_influence_graph)(file, None)


def test_read_influence_graph_accepts_named_file(tmp_path):
    file = tmp_path / "graph.csv"
    file.write_text("source,target,sign\nA,B,1\n")

    graph = bt.logic.io.read_influence_graph(file=file)

    assert list(graph.edges(data=True)) == [("A", "B", {"sign": 1.0})]


def test_read_influence_graph_validates_file_columns_and_signs(tmp_path):
    missing = tmp_path / "missing.csv"

    with pytest.raises(FileNotFoundError):
        bt.logic.io.read_influence_graph(missing)

    missing_column = tmp_path / "missing_column.csv"
    missing_column.write_text("source,target\nA,B\n")

    with pytest.raises(ValueError, match="missing columns"):
        bt.logic.io.read_influence_graph(missing_column)

    invalid_sign = tmp_path / "invalid_sign.csv"
    invalid_sign.write_text("source,target,sign\nA,B,0\n")

    with pytest.raises(ValueError, match="unsupported sign values"):
        bt.logic.io.read_influence_graph(invalid_sign)


def test_read_influence_graph_rejects_invalid_genesyn(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        bt.logic.io.read_influence_graph(infile, genesyn=cast(Any, object()))


def test_read_influence_graph_applies_genesyn(monkeypatch, tmp_path):
    class FakeGeneSynonyms:
        def __init__(self):
            self.calls = []

        def __call__(self, graph, **kwargs):
            self.calls.append((graph, kwargs))

    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")
    genesyn = FakeGeneSynonyms()

    monkeypatch.setattr(_influence_graph, "GeneSynonyms", FakeGeneSynonyms)

    graph = bt.logic.io.read_influence_graph(
        infile,
        genesyn=cast(Any, genesyn),
        input_identifier_type="gene_id",
        output_identifier_type="ensembl_id",
    )

    assert genesyn.calls == [
        (
            graph,
            {
                "input_identifier_type": "gene_id",
                "output_identifier_type": "ensembl_id",
                "copy": False,
            },
        )
    ]


def test_deprecated_read_influence_graph_routes_to_io(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")

    with pytest.warns(FutureWarning, match="bt.logic.ig.read_influence_graph"):
        graph = cast(nx.MultiDiGraph, bt.logic.ig.read_influence_graph(infile))

    expected = bt.logic.io.read_influence_graph(infile)
    assert list(graph.edges(data=True)) == list(expected.edges(data=True))
