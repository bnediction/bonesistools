#!/usr/bin/env python

import math
from typing import Any, cast

import networkx as nx
import pytest

import bonesistools as bt


def _edge_data(graph, source, target, key=0):
    return cast(Any, graph[source][target])[key]


def test_read_influence_graph_keeps_edge_attributes(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text(
        "source,target,sign,confidence\nA,B,1,high\nB,C,-1,low\nC,A,,unknown\n",
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
        "source\ttarget\tsign\nA\tB\t1\n",
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


def test_read_influence_graph_rejects_invalid_identifiers(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")

    with pytest.raises(TypeError, match="unsupported argument type for 'identifiers'"):
        bt.logic.io.read_influence_graph(infile, identifiers=cast(Any, object()))


def test_read_influence_graph_applies_identifiers(tmp_path):
    class FakeGeneIdentifiers:
        def __init__(self):
            self.calls = []

        def __call__(self, graph, **kwargs):
            self.calls.append((graph, kwargs))

    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")
    identifiers = FakeGeneIdentifiers()

    graph = bt.logic.io.read_influence_graph(
        infile,
        identifiers=cast(Any, identifiers),
        input_type="gene_id",
        output_type="ensembl_id",
    )

    assert identifiers.calls == [
        (
            graph,
            {
                "input_type": "gene_id",
                "output_type": "ensembl_id",
                "copy": False,
            },
        )
    ]


def test_read_influence_graph_accepts_deprecated_identifier_arguments(
    tmp_path,
):
    class FakeGeneIdentifiers:
        def __call__(self, graph, **kwargs):
            self.kwargs = kwargs

    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")
    identifiers = FakeGeneIdentifiers()

    with pytest.warns(FutureWarning) as warning_records:
        cast(Any, bt.logic.io.read_influence_graph)(
            infile,
            genesyn=cast(Any, identifiers),
            input_identifier_type="gene_id",
            output_identifier_type="ensembl_id",
        )

    messages = [str(record.message) for record in warning_records]
    assert any(
        "genesyn" in message and "identifiers" in message for message in messages
    )
    assert any("input_identifier_type" in message for message in messages)
    assert any("output_identifier_type" in message for message in messages)
    assert identifiers.kwargs == {
        "input_type": "gene_id",
        "output_type": "ensembl_id",
        "copy": False,
    }


def test_deprecated_read_influence_graph_routes_to_io(tmp_path):
    infile = tmp_path / "graph.csv"
    infile.write_text("source,target,sign\nA,B,1\n")

    with pytest.warns(FutureWarning, match="bt.logic.ig.read_influence_graph"):
        graph = cast(nx.MultiDiGraph, bt.logic.ig.read_influence_graph(infile))

    expected = bt.logic.io.read_influence_graph(infile)
    assert list(graph.edges(data=True)) == list(expected.edges(data=True))
