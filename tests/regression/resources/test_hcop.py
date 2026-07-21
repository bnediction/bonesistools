#!/usr/bin/env python

from types import MappingProxyType
from typing import Any, cast

import networkx as nx
import pandas as pd
import pytest

import bonesistools as bt
from bonesistools.resources import hcop
from bonesistools.resources.hcop import _orthologs as _hcop_orthologs


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

    translated = hcop.orthologs(target_organism="mouse").translate_dataframe(
        net,
        columns=["source", "target"],
        on_unmapped="drop",
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

    translated = orthologs.translate_dataframe(
        net,
        columns=["source", "target"],
    )

    assert translated.to_dict("records") == [
        {"source": "a_b1", "target": "t", "row": 1},
    ]


def test_hcop_translate_dataframe_expands_one_to_many_like_decoupler(monkeypatch):
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

    translated = hcop.orthologs(target_organism="mouse").translate_dataframe(
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


def test_hcop_one_to_many_mapping_is_refreshed_after_reset():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "A", "T"],
                "target_symbol": ["a1", "a2", "t"],
                "support": ["Ensembl,HGNC,MGI"] * 3,
                "evidence": [3, 3, 3],
            }
        ),
        target_organism="mouse",
    )
    dataframe = pd.DataFrame({"source": ["A"], "target": ["T"], "row": [1]})

    before_reset = orthologs.translate_dataframe(
        dataframe,
        columns=["source", "target"],
        one_to_many=2,
    )

    orthologs.reset(
        target_organism="mouse",
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "T", "T"],
                "target_symbol": ["new_a", "new_t1", "new_t2"],
                "support": ["Ensembl,HGNC,MGI"] * 3,
                "evidence": [3, 3, 3],
            }
        ),
    )
    after_reset = orthologs.translate_dataframe(
        dataframe,
        columns=["source", "target"],
        one_to_many=2,
    )

    assert before_reset.to_dict("records") == [
        {"source": "a1", "target": "t", "row": 1},
        {"source": "a2", "target": "t", "row": 1},
    ]
    assert after_reset.to_dict("records") == [
        {"source": "new_a", "target": "new_t1", "row": 1},
        {"source": "new_a", "target": "new_t2", "row": 1},
    ]


def test_hcop_dataframe_translation_preserves_complex_and_missing_values():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B"],
                "target_symbol": ["a", "b"],
                "support": ["source"] * 2,
                "evidence": [1, 1],
            }
        ),
        target_organism="mouse",
    )
    dataframe = pd.DataFrame({"genesymbol": ["A", "A_B", "A_missing", "missing", 1]})

    with pytest.raises(ValueError, match="no target ortholog"):
        orthologs.translate_dataframe(dataframe)

    kept = orthologs.translate_dataframe(dataframe, on_unmapped="keep")
    removed = orthologs.translate_dataframe(dataframe, on_unmapped="drop")
    with_nan = orthologs.translate_dataframe(dataframe, on_unmapped="nan")

    assert kept["genesymbol"].tolist() == [
        "a",
        "a_b",
        "a_missing",
        "missing",
        "1",
    ]
    assert removed["genesymbol"].tolist() == ["a", "a_b"]
    assert with_nan["genesymbol"].iloc[:2].tolist() == ["a", "a_b"]
    assert with_nan["genesymbol"].iloc[2:].isna().all()


def test_hcop_dataframe_translation_drops_unmapped_rows_in_place():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A"],
                "target_symbol": ["a"],
                "support": ["source"],
                "evidence": [1],
            }
        ),
        target_organism="mouse",
    )
    dataframe = pd.DataFrame(
        {
            "genesymbol": ["A", "missing"],
            "row": [1, 2],
        },
        index=[4, 8],
    )

    result = orthologs.translate_dataframe(
        dataframe,
        on_unmapped="drop",
        copy=False,
    )

    assert result is None
    assert dataframe.to_dict("records") == [{"genesymbol": "a", "row": 1}]
    assert dataframe.index.tolist() == [0]


def test_hcop_dataframe_one_to_many_honors_unmapped_policy():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "A", "B", "B", "B"],
                "target_symbol": ["a1", "a2", "b1", "b2", "b3"],
                "support": ["source"] * 5,
                "evidence": [1] * 5,
            }
        ),
        target_organism="mouse",
    )
    dataframe = pd.DataFrame(
        {
            "genesymbol": ["A", "missing", "A_missing", "B"],
            "row": [1, 2, 3, 4],
        }
    )

    with pytest.raises(ValueError, match="no target ortholog"):
        orthologs.translate_dataframe(
            dataframe,
            one_to_many=2,
        )

    kept = orthologs.translate_dataframe(
        dataframe,
        on_unmapped="keep",
        one_to_many=2,
    )
    dropped = orthologs.translate_dataframe(
        dataframe,
        on_unmapped="drop",
        one_to_many=2,
    )
    with_nan = orthologs.translate_dataframe(
        dataframe,
        on_unmapped="nan",
        one_to_many=2,
    )

    assert kept.to_dict("records") == [
        {"genesymbol": "a1", "row": 1},
        {"genesymbol": "a2", "row": 1},
        {"genesymbol": "missing", "row": 2},
        {"genesymbol": "a1_missing", "row": 3},
        {"genesymbol": "a2_missing", "row": 3},
    ]
    assert dropped.to_dict("records") == [
        {"genesymbol": "a1", "row": 1},
        {"genesymbol": "a2", "row": 1},
    ]
    assert with_nan["row"].tolist() == [1, 1, 2, 3]
    assert with_nan["genesymbol"].iloc[:2].tolist() == ["a1", "a2"]
    assert with_nan["genesymbol"].iloc[2:].isna().all()


def test_hcop_one_to_many_preserves_source_row_order_and_duplicates():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "A", "A"],
                "target_symbol": ["a2", "a1", "a2"],
                "support": ["source"] * 3,
                "evidence": [1, 1, 1],
            }
        ),
        target_organism="mouse",
    )

    translated = orthologs.translate_dataframe(
        pd.DataFrame({"genesymbol": ["A"]}),
        one_to_many=3,
    )

    assert translated["genesymbol"].tolist() == ["a2", "a1", "a2"]


def test_hcop_best_mapping_is_refreshed_after_reset():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A"],
                "target_symbol": ["a"],
                "support": ["source"],
                "evidence": [1],
            }
        ),
        target_organism="mouse",
    )
    assert orthologs.translate_sequence(["A"]) == ["a"]

    orthologs.reset(
        table=pd.DataFrame(
            {
                "human_symbol": ["A"],
                "target_symbol": ["new_a"],
                "support": ["source"],
                "evidence": [1],
            }
        )
    )

    assert orthologs.translate_sequence(["A"]) == ["new_a"]


def test_hcop_mapping_order_is_deterministic():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["B", "A", "B"],
                "target_symbol": ["b2", "a", "b1"],
                "support": ["source"] * 3,
                "evidence": [1, 1, 2],
            }
        ),
        target_organism="mouse",
    )

    assert list(orthologs.to_dict().items()) == [
        ("A", ["a"]),
        ("B", ["b1", "b2"]),
    ]


def test_hcop_translates_between_non_human_organisms_with_scores(monkeypatch):
    def read_hcop_table(output_organism, version="bundled"):
        if output_organism == "mouse":
            return pd.DataFrame(
                {
                    "human_symbol": ["H1", "H2", "-"],
                    "mouse_symbol": ["M", "M", "placeholder_mouse"],
                    "support": ["a,b,c,d", "a,b", "a,b,c"],
                }
            )
        if output_organism == "rat":
            return pd.DataFrame(
                {
                    "human_symbol": ["H1", "H1", "H2", "H2", "-"],
                    "rat_symbol": ["R1", "R2", "R2", "R3", "placeholder_rat"],
                    "support": [
                        "a,b,c,d",
                        "a,b,c,d,e",
                        "a,b,c,d",
                        "a,b,c",
                        "a,b,c",
                    ],
                }
            )
        raise AssertionError(output_organism)

    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(read_hcop_table),
    )

    orthologs = hcop.orthologs(
        input_organism="mouse",
        output_organism="rat",
        min_evidence=1,
    )

    assert orthologs.translate("M") == ["R2", "R1", "R3"]
    assert orthologs.translate("placeholder_mouse", on_unmapped="drop") == []

    paths = orthologs.to_dataframe()
    r2_paths = paths.loc[paths["target_symbol"] == "R2"]
    assert set(r2_paths["human_symbol"]) == {"H1", "H2"}
    assert set(r2_paths["evidence"]) == {2, 4}
    assert set(r2_paths["best_evidence"]) == {4}
    assert set(r2_paths["paths"]) == {2}

    r1_paths = paths.loc[paths["target_symbol"] == "R1"]
    assert set(r1_paths["best_evidence"]) == {4}
    assert set(r1_paths["paths"]) == {1}

    translated = orthologs.translate_dataframe(
        pd.DataFrame({"source": ["M"], "row": [1]}),
        columns="source",
        one_to_many=3,
    )
    assert translated.to_dict("records") == [
        {"source": "R2", "row": 1},
        {"source": "R1", "row": 1},
        {"source": "R3", "row": 1},
    ]


def test_hcop_translates_non_human_genes_to_human(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda output_organism, version="bundled": pd.DataFrame(
                {
                    "human_symbol": ["H1", "H2", "-"],
                    "mouse_symbol": ["M", "M", "placeholder_mouse"],
                    "support": ["a,b,c,d", "a,b", "a,b,c"],
                }
            )
        ),
    )

    orthologs = hcop.orthologs(
        input_organism="mouse",
        output_organism="human",
        min_evidence=1,
    )

    assert orthologs.translate("M") == ["H1", "H2"]
    assert orthologs.translate("placeholder_mouse", on_unmapped="drop") == []
    assert list(orthologs.to_dataframe().columns) == [
        "source_symbol",
        "target_symbol",
        "support",
        "evidence",
    ]


def test_hcop_filters_both_branches_of_non_human_mapping(monkeypatch):
    def read_hcop_table(output_organism, version="bundled"):
        if output_organism == "mouse":
            return pd.DataFrame(
                {
                    "human_symbol": ["H1", "H2"],
                    "mouse_symbol": ["M", "M"],
                    "support": ["a,b,c,d", "a,b"],
                }
            )
        return pd.DataFrame(
            {
                "human_symbol": ["H1", "H1", "H2", "H2"],
                "rat_symbol": ["R1", "R2", "R2", "R3"],
                "support": ["a,b,c", "a,b,c,d,e", "a,b,c,d", "a,b,c"],
            }
        )

    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(read_hcop_table),
    )

    orthologs = hcop.orthologs(
        input_organism="mouse",
        output_organism="rat",
        min_evidence=3,
    )

    assert orthologs.translate("M") == ["R2", "R1"]
    assert set(orthologs.to_dataframe()["human_symbol"]) == {"H1"}


def test_hcop_accepts_direct_custom_non_human_mapping():
    orthologs = _hcop_orthologs.Orthologs(
        input_organism="mouse",
        output_organism="rat",
        table=pd.DataFrame(
            {
                "source_symbol": ["M"],
                "target_symbol": ["R"],
                "support": ["custom"],
                "evidence": [1],
            }
        ),
    )

    assert orthologs.translate("M") == ["R"]


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
    assert orthologs(("A", "missing"), on_unmapped="keep") == ("a", "missing")
    assert orthologs.translate_sequence(("A", "missing"), on_unmapped="drop") == ("a",)
    assert orthologs.translate_sequence({"A", "B"}) == {"a", "b"}
    assert orthologs.translate("A_B") == ["a_b"]

    interactions = [("A", "B", {"sign": 1}), ("A", "missing", {"sign": -1})]
    assert orthologs(interactions, on_unmapped="drop") == [("a", "b", {"sign": 1})]

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

    translated = orthologs.translate_graph(graph, on_unmapped="drop", copy=True)

    assert translated is not None
    assert set(translated.nodes) == {"a", "b"}
    assert translated.nodes["a"]["role"] == "source"
    assert translated["a"]["b"]["sign"] == 1
    assert set(graph.nodes) == {"A", "B", "missing"}

    assert orthologs.translate_graph(graph, on_unmapped="keep", copy=False) is None
    assert set(graph.nodes) == {"a", "b", "missing"}


def test_hcop_translate_hypercube_copy_dispatch_and_missing_behaviors():
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
    hypercube = bt.logic.ba.Hypercube({"A": 1, "B": 0, "missing": 1})

    copied = orthologs.translate_hypercube(hypercube, on_unmapped="keep")
    dispatched = orthologs(hypercube, on_unmapped="keep")

    assert type(copied) is type(hypercube)
    assert copied == {"a": 1, "b": 0, "missing": 1}
    assert dispatched == copied
    assert hypercube == {"A": 1, "B": 0, "missing": 1}

    assert (
        orthologs.translate_hypercube(
            hypercube,
            on_unmapped="keep",
            copy=False,
        )
        is None
    )
    assert hypercube == {"a": 1, "b": 0, "missing": 1}

    with pytest.raises(ValueError, match="cannot remove an unmapped component"):
        orthologs.translate_hypercube(
            bt.logic.ba.Hypercube({"A": 1, "missing": 0}),
            on_unmapped="raise",
        )


def test_hcop_translate_hypercube_preserves_hypercube_like_types():
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
    mapping = {"A": 1, "B": "*"}
    immutable_mapping = MappingProxyType({"A": 0})

    translated_mapping = orthologs.translate_hypercube(mapping)
    translated_immutable = orthologs.translate_hypercube(immutable_mapping)

    assert type(translated_mapping) is dict
    assert translated_mapping == {"a": 1, "b": "*"}
    assert mapping == {"A": 1, "B": "*"}
    assert type(translated_immutable) is type(immutable_mapping)
    assert dict(translated_immutable) == {"a": 0}
    assert dict(immutable_mapping) == {"A": 0}

    assert orthologs.translate_hypercube(mapping, copy=False) is None
    assert mapping == {"a": 1, "b": "*"}

    with pytest.raises(TypeError, match="mutable hypercube-like mapping"):
        orthologs.translate_hypercube(immutable_mapping, copy=False)


def test_hcop_translate_hypercube_rejects_component_collisions():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B"],
                "target_symbol": ["shared", "shared"],
                "support": ["Ensembl,HGNC,MGI", "Ensembl,HGNC,MGI"],
                "evidence": [3, 3],
            }
        ),
        target_organism="mouse",
    )

    with pytest.raises(ValueError, match="merge hypercube components"):
        orthologs.translate_hypercube(bt.logic.ba.Hypercube({"A": 1, "B": 0}))


def test_hcop_translate_boolean_network_copy_and_mapping_behaviors(monkeypatch):
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

    relabel_calls = []
    original_relabel = bt.logic.bn.BooleanNetwork.relabel

    def record_relabel(network, mapping):
        relabel_calls.append(dict(mapping))
        return original_relabel(network, mapping)

    monkeypatch.setattr(
        bt.logic.bn.BooleanNetwork,
        "relabel",
        record_relabel,
    )

    bn = bt.logic.bn.BooleanNetwork({"A": "B", "B": 1})
    copied = orthologs.translate_boolean_network(bn, copy=True)

    assert isinstance(copied, bt.logic.bn.BooleanNetwork)
    assert copied.rules == {"a": "b", "b": "1"}
    assert bn.rules == {"A": "B", "B": "1"}
    assert relabel_calls == [{"A": "a", "B": "b"}]

    mapping_bn = {"A": "B", "B": 1}
    translated_mapping = orthologs.translate_boolean_network(mapping_bn, copy=True)

    assert translated_mapping == {"a": "b", "b": "1"}

    with pytest.raises(TypeError, match="copy=True"):
        orthologs.translate_boolean_network(mapping_bn, copy=False)

    with pytest.raises(ValueError, match="unmapped component"):
        orthologs.translate_boolean_network(
            bt.logic.bn.BooleanNetwork({"A": "missing", "missing": 1}),
            on_unmapped="raise",
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
        orthologs.translate_dataframe(cast(Any, "not a dataframe"))

    with pytest.raises(TypeError, match="unsupported argument type for 'hypercube'"):
        orthologs.translate_hypercube(cast(Any, [("A", 1)]))

    with pytest.raises(ValueError, match="output_organism"):
        hcop.orthologs(target_organism="frog")

    with pytest.raises(ValueError, match="input_organism"):
        hcop.orthologs(input_organism="frog")

    with pytest.raises(ValueError, match="must differ"):
        hcop.orthologs(input_organism="mouse", output_organism="mouse")

    with pytest.raises(ValueError, match="columns"):
        orthologs.translate_dataframe(
            pd.DataFrame({"other": ["A"]}),
            columns=["source"],
        )

    with pytest.raises(ValueError, match="one_to_many"):
        orthologs.translate_dataframe(
            pd.DataFrame({"source": ["A"]}),
            columns=["source"],
            one_to_many=0,
        )

    with pytest.raises(TypeError, match="copy=True"):
        orthologs.translate_dataframe(
            pd.DataFrame({"source": ["A"]}),
            columns=["source"],
            copy=False,
            one_to_many=1,
        )

    with pytest.raises(ValueError, match="on_unmapped"):
        orthologs.translate_dataframe(
            pd.DataFrame({"source": ["A"]}),
            on_unmapped=cast(Any, "invalid"),
        )

    with pytest.raises(ValueError, match="on_unmapped"):
        orthologs.translate_hypercube(
            bt.logic.ba.Hypercube({"A": 1}),
            on_unmapped=cast(Any, "nan"),
        )

    with pytest.raises(TypeError):
        cast(Any, orthologs)("A", False)


def test_hcop_keep_if_missing_argument_is_deprecated():
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
    )
    legacy_translate = cast(Any, orthologs.translate_dataframe)
    legacy_translate_sequence = cast(Any, orthologs.translate_sequence)

    with pytest.warns(FutureWarning, match="on_unmapped"):
        translated = legacy_translate(
            pd.DataFrame({"source": ["A", "missing"]}),
            keep_if_missing=False,
        )
    with pytest.warns(FutureWarning, match="on_unmapped"):
        translated_sequence = legacy_translate_sequence(
            ("A", "missing"),
            False,
        )

    assert translated.to_dict("records") == [{"source": "a"}]
    assert translated_sequence == ("a",)


def test_hcop_deprecated_translation_aliases_warn():
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
    )

    with pytest.warns(FutureWarning, match="translate_dataframe"):
        dataframe = orthologs.translate_df(pd.DataFrame({"source": ["A"]}))
    assert dataframe.to_dict("records") == [{"source": "a"}]

    network = bt.logic.bn.BooleanNetwork({"A": 1})
    with pytest.warns(FutureWarning, match="translate_boolean_network"):
        translated = orthologs.translate_bn(network, copy=True)
    assert isinstance(translated, bt.logic.bn.BooleanNetwork)
    assert translated.rules == {"a": "1"}


def test_hcop_organisms_lists_supported_targets():
    assert "mouse" in hcop.organisms()
    assert "rat" in hcop.organisms()


def test_hcop_deprecated_ortholog_versions_warns():
    deprecated_versions = getattr(hcop.orthologs, "versions")

    with pytest.warns(FutureWarning, match="hcop.versions"):
        result = deprecated_versions("mouse")

    assert result == hcop.versions()


def test_hcop_orthologs_loads_hcop_mapping(monkeypatch):
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
                        "Ensembl,HGNC",
                        "Ensembl,HGNC,MGI,NCBI",
                    ],
                }
            )
        ),
    )

    orthologs = hcop.orthologs(target_organism="mouse", min_evidence=3)

    assert orthologs.to_dataframe().to_dict("records") == [
        {
            "human_symbol": "A",
            "target_symbol": "A_mouse",
            "support": "Ensembl,HGNC,MGI",
            "evidence": 3,
        },
        {
            "human_symbol": "C",
            "target_symbol": "C_mouse",
            "support": "Ensembl,HGNC,MGI,NCBI",
            "evidence": 4,
        },
    ]


def test_hcop_orthologs_dataframe_and_dict_are_deep_copies(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda target_organism, version="bundled": pd.DataFrame(
                {
                    "human_symbol": ["A", "B"],
                    "mouse_symbol": ["A_mouse", "B_mouse"],
                    "support": ["Ensembl,HGNC,MGI", "Ensembl,HGNC,MGI"],
                }
            )
        ),
    )

    orthologs = hcop.orthologs(target_organism="mouse")

    table = orthologs.to_dataframe()
    table.loc[0, "target_symbol"] = "mutated"
    mapping = orthologs.to_dict()
    mapping["A"].append("mutated")

    assert orthologs.to_dataframe().loc[0, "target_symbol"] == "A_mouse"
    assert orthologs.to_dict() == {"A": ["A_mouse"], "B": ["B_mouse"]}


def test_hcop_orthologs_rejects_invalid_table(monkeypatch):
    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(
            lambda target_organism, version="bundled": pd.DataFrame(
                {"human_symbol": ["A"]}
            )
        ),
    )

    with pytest.raises(ValueError, match="missing required columns"):
        hcop.orthologs(target_organism="mouse")


def test_hcop_orthologs_loads_bundled_snapshot(tmp_path, monkeypatch):
    bundled_file = tmp_path / "human_mouse_hcop.tsv"
    bundled_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\nA\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(_hcop_orthologs, "HCOP_DIR", tmp_path)

    orthologs = hcop.orthologs(target_organism="mouse")

    assert orthologs.version == "bundled"
    assert orthologs.translate("A") == ["A_mouse"]
    assert hcop.versions() == ["bundled", "latest"]


def test_hcop_orthologs_loads_custom_file(tmp_path):
    hcop_file = tmp_path / "custom_hcop.tsv"
    hcop_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\nA\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    orthologs = hcop.orthologs(target_organism="mouse", version=hcop_file)

    assert orthologs.version == str(hcop_file)
    assert orthologs.translate("A") == ["A_mouse"]


def test_hcop_latest_table_uses_local_download_cache(tmp_path, monkeypatch):
    latest_file = tmp_path / "human_mouse_hcop.tsv"
    latest_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\nA\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )
    calls = []

    def cached_download(url, **kwargs):
        calls.append((url, kwargs))
        return latest_file

    monkeypatch.setattr(_hcop_orthologs, "_cached_download", cached_download)

    table = _hcop_orthologs.Orthologs._read_hcop_table(
        "mouse",
        version="latest",
    )

    assert calls == [
        (
            "https://storage.googleapis.com/public-download-files/hcop/"
            "human_mouse_hcop_fifteen_column.txt.gz",
            {
                "resource": "hcop",
                "category": "tables",
                "max_age": 72 * 60 * 60,
            },
        )
    ]
    assert table.to_dict("records") == [
        {
            "human_symbol": "A",
            "mouse_symbol": "A_mouse",
            "support": "Ensembl,HGNC,MGI",
        }
    ]


def test_hcop_non_human_mapping_rejects_single_custom_hcop_file(tmp_path):
    hcop_file = tmp_path / "custom_hcop.tsv"
    hcop_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\nA\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="requires two HCOP tables"):
        hcop.orthologs(
            input_organism="mouse",
            output_organism="rat",
            version=hcop_file,
        )


def test_hcop_orthologs_reset_reconfigures_existing_converter():
    orthologs = _hcop_orthologs.Orthologs(
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "B"],
                "target_symbol": ["A_mouse", "B_mouse"],
                "support": ["Ensembl,HGNC,MGI", "Ensembl,HGNC,MGI"],
                "evidence": [3, 3],
            }
        ),
        target_organism="mouse",
        min_evidence=3,
    )

    orthologs.reset(
        target_organism="rat",
        min_evidence=2,
        version="",
        table=pd.DataFrame(
            {
                "human_symbol": ["A", "C"],
                "target_symbol": ["A_rat", "C_rat"],
                "support": ["Ensembl,HGNC", "Ensembl,HGNC"],
                "evidence": [2, 2],
            }
        ),
    )

    assert orthologs.input_organism == "human"
    assert orthologs.output_organism == "rat"
    assert orthologs.target_organism == "rat"
    assert orthologs.min_evidence == 2
    assert orthologs.version == "bundled"
    assert orthologs.translate("A") == ["A_rat"]
    assert orthologs.translate("B", on_unmapped="drop") == []
    assert orthologs.find("A", "B", "C") == ["A", "C"]
    assert orthologs.contains("A", "B", "C") == [True, False, True]


def test_hcop_orthologs_reads_special_target_symbol_columns(monkeypatch):
    def read_hcop_table(target_organism, version="bundled"):
        if target_organism == "anole_lizard":
            return pd.DataFrame(
                {
                    "human_symbol": ["A"],
                    "anole lizard_symbol": ["A_anole"],
                    "support": ["Ensembl,HGNC,MGI"],
                }
            )
        if target_organism == "fruitfly":
            return pd.DataFrame(
                {
                    "human_symbol": ["B"],
                    "fruit fly_symbol": ["B_fly"],
                    "support": ["Ensembl,HGNC,MGI"],
                }
            )
        raise AssertionError(target_organism)

    monkeypatch.setattr(
        _hcop_orthologs.Orthologs,
        "_read_hcop_table",
        staticmethod(read_hcop_table),
    )

    anole = hcop.orthologs(target_organism="anole_lizard")
    fruitfly = hcop.orthologs(target_organism="fruitfly")

    assert anole.to_dataframe().to_dict("records") == [
        {
            "human_symbol": "A",
            "target_symbol": "A_anole",
            "support": "Ensembl,HGNC,MGI",
            "evidence": 3,
        }
    ]
    assert fruitfly.to_dataframe().to_dict("records") == [
        {
            "human_symbol": "B",
            "target_symbol": "B_fly",
            "support": "Ensembl,HGNC,MGI",
            "evidence": 3,
        }
    ]


def test_hcop_orthologs_missing_version_lists_versions():
    with pytest.raises(
        FileNotFoundError,
        match="Available versions: bundled, latest",
    ):
        hcop.orthologs(target_organism="mouse", version="missing_hcop.tsv")


def test_hcop_orthologs_missing_bundled_file_fails_clearly(tmp_path, monkeypatch):
    monkeypatch.setattr(_hcop_orthologs, "HCOP_DIR", tmp_path)
    with pytest.raises(
        FileNotFoundError,
        match="bundled HCOP file not found",
    ):
        hcop.orthologs(target_organism="mouse")
