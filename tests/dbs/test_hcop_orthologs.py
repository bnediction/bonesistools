#!/usr/bin/env python

import pandas as pd
import pytest

from bonesistools.databases import hcop
from bonesistools.databases.hcop import _orthologs as _hcop_orthologs


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
        "human_symbol\tmouse_symbol\tsupport\n" "A\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(_hcop_orthologs, "HCOP_DIR", tmp_path)

    orthologs = hcop.orthologs(target_organism="mouse")

    assert orthologs.version == "bundled"
    assert orthologs.translate("A") == ["A_mouse"]
    assert hcop.orthologs.versions("mouse") == ["bundled", "latest"]


def test_hcop_orthologs_loads_custom_file(tmp_path):
    hcop_file = tmp_path / "custom_hcop.tsv"
    hcop_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\n" "A\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    orthologs = hcop.orthologs(target_organism="mouse", version=hcop_file)

    assert orthologs.version == str(hcop_file)
    assert orthologs.translate("A") == ["A_mouse"]


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
    assert orthologs.translate("B") == []
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
