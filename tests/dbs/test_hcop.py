#!/usr/bin/env python

from typing import Any, cast

import networkx as nx
import pandas as pd
import pytest

import bonesistools as bt
from bonesistools.databases import hcop
from bonesistools.databases.hcop import _orthologs as _hcop_orthologs


def test_hcop_translate_translates_requested_columns(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda target_organism, version="bundled": pd.DataFrame(
                {
                    "human_symbol": ["A", "B", "C"],
                    "mouse_symbol": ["A_mouse", "B_mouse", "C_mouse"],
                    "support": [
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl",
                    ],
                }
            )
        ),
    )

    net = pd.DataFrame(
        {
            "source": ["A", "C", "A"],
            "target": ["B", "B", "missing"],
            "sign": [1, -1, 1],
        }
    )

    translated = hcop.orthologs(target_organism="mouse").translate_df(
        net,
        columns=["source", "target"],
        keep_if_missing=False,
    )

    assert translated.to_dict("records") == [
        {"source": "A_mouse", "target": "B_mouse", "sign": 1}
    ]


def test_hcop_translate_prioritizes_evidence_then_frequency(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda target_organism, version="bundled": pd.DataFrame(
                {
                    "human_symbol": ["A", "B", "B", "B", "T"],
                    "mouse_symbol": ["a", "b1", "b2", "b2", "t"],
                    "support": [
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI,NCBI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                    ],
                }
            )
        ),
    )

    net = pd.DataFrame(
        {
            "source": ["A_B"],
            "target": ["T"],
            "row": [1],
        }
    )

    orthologs = hcop.orthologs(target_organism="mouse")
    assert orthologs.translate("B") == ["b1", "b2"]

    translated = orthologs.translate_df(
        net,
        columns=["source", "target"],
    )

    assert translated.to_dict("records") == [
        {"source": "a_b1", "target": "t", "row": 1},
    ]


def test_hcop_translate_df_expands_one_to_many_like_decoupler(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda target_organism, version="bundled": pd.DataFrame(
                {
                    "human_symbol": ["A", "A", "B", "B", "B", "D", "T"],
                    "mouse_symbol": ["a1", "a2", "b1", "b2", "b3", "d", "t"],
                    "support": [
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                        "Ensembl,HGNC,MGI",
                    ],
                }
            )
        ),
    )

    net = pd.DataFrame(
        {
            "source": ["A", "A_D", "A_B"],
            "target": ["T", "T", "T"],
            "row": [1, 2, 3],
        }
    )

    translated = hcop.orthologs(target_organism="mouse").translate_df(
        net,
        columns=["source", "target"],
        one_to_many=2,
    )

    assert translated.to_dict("records") == [
        {"source": "a1", "target": "t", "row": 1},
        {"source": "a2", "target": "t", "row": 1},
        {"source": "a1_d", "target": "t", "row": 2},
        {"source": "a2_d", "target": "t", "row": 2},
    ]


def test_hcop_translates_sequences_and_dispatches_supported_objects():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B", "C"],
                "target_symbol": ["a", "b", "c"],
                "support": [
                    "Ensembl,HGNC,MGI",
                    "Ensembl,HGNC,MGI",
                    "Ensembl,HGNC,MGI",
                ],
                "evidence": [3, 3, 3],
            }
        ),
        target_organism="mouse",
    )

    assert orthologs("A") == ["a"]
    assert orthologs(("A", "missing")) == ("a", "missing")
    assert orthologs.translate_sequence(("A", "missing"), keep_if_missing=False) == (
        "a",
    )
    assert orthologs.translate_sequence({"A", "B"}) == {"a", "b"}
    assert orthologs.translate("A_B") == ["a_b"]

    interactions = [("A", "B", {"sign": 1}), ("A", "missing", {"sign": -1})]
    assert orthologs(interactions, keep_if_missing=False) == [
        ("a", "b", {"sign": 1})
    ]

    df = pd.DataFrame({"source": ["A"], "target": ["B"]})
    assert orthologs(df, columns=["source", "target"]).to_dict("records") == [
        {"source": "a", "target": "b"}
    ]


def test_hcop_translate_graph_preserves_attributes_and_copy_semantics():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B"],
                "target_symbol": ["a", "b"],
                "support": ["Ensembl,HGNC,MGI", "Ensembl,HGNC,MGI"],
                "evidence": [3, 3],
            }
        ),
        target_organism="mouse",
    )

    graph = nx.Graph()
    graph.add_node("A", role="source")
    graph.add_node("B", role="target")
    graph.add_node("missing", role="unknown")
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "missing", sign=-1)

    translated = orthologs.translate_graph(graph, keep_if_missing=False, copy=True)

    assert translated is not None
    assert set(translated.nodes) == {"a", "b"}
    assert translated.nodes["a"]["role"] == "source"
    assert translated["a"]["b"]["sign"] == 1
    assert set(graph.nodes) == {"A", "B", "missing"}

    assert orthologs.translate_graph(graph, copy=False) is None
    assert set(graph.nodes) == {"a", "b", "missing"}


def test_hcop_translate_boolean_network_copy_and_mapping_behaviors():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B"],
                "target_symbol": ["a", "b"],
                "support": ["Ensembl,HGNC,MGI", "Ensembl,HGNC,MGI"],
                "evidence": [3, 3],
            }
        ),
        target_organism="mouse",
    )

    bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": 1})
    copied = orthologs.translate_bn(bn, copy=True)

    assert isinstance(copied, bt.bpy.bn.BooleanNetwork)
    assert copied.rules == {"a": "b", "b": "1"}
    assert bn.rules == {"A": "B", "B": "1"}

    mapping_bn = {"A": "B", "B": 1}
    translated_mapping = orthologs.translate_bn(mapping_bn, copy=True)

    assert translated_mapping == {"a": "b", "b": "1"}

    with pytest.raises(TypeError, match="copy=True"):
        orthologs.translate_bn(mapping_bn, copy=False)

    with pytest.raises(ValueError, match="unmapped component"):
        orthologs.translate_bn(
            bt.bpy.bn.BooleanNetwork({"A": "missing", "missing": 1}),
            keep_if_missing=False,
        )


def test_hcop_translate_rejects_invalid_arguments():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A"],
                "target_symbol": ["a"],
                "support": ["Ensembl,HGNC,MGI"],
                "evidence": [3],
            }
        ),
        target_organism="mouse",
        min_evidence=3,
    )

    with pytest.raises(TypeError, match="unsupported argument type for 'df'"):
        orthologs.translate_df(cast(Any, "not a dataframe"))

    with pytest.raises(ValueError, match="output_organism"):
        hcop.orthologs(target_organism="frog")

    with pytest.raises(ValueError, match="columns"):
        orthologs.translate_df(pd.DataFrame({"other": ["A"]}), columns=["source"])

    with pytest.raises(ValueError, match="one_to_many"):
        orthologs.translate_df(
            pd.DataFrame({"source": ["A"]}),
            columns=["source"],
            one_to_many=0,
        )

    with pytest.raises(TypeError, match="copy=True"):
        orthologs.translate_df(
            pd.DataFrame({"source": ["A"]}),
            columns=["source"],
            copy=False,
            one_to_many=1,
        )


def test_hcop_organisms_lists_supported_targets():
    assert "mouse" in hcop.organisms()
    assert "rat" in hcop.organisms()
