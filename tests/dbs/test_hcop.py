#!/usr/bin/env python

import pandas as pd
import pytest

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
        orthologs.translate_df("not a dataframe")

    with pytest.raises(ValueError, match="output_organism"):
        hcop.orthologs(target_organism="frog")

    with pytest.raises(ValueError, match="columns"):
        orthologs.translate_df(pd.DataFrame({"other": ["A"]}), columns=["source"])

    with pytest.raises(TypeError, match="one_to_many"):
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
