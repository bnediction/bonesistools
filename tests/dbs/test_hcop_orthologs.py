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
    bundled_dir = tmp_path / "bundled"
    bundled_dir.mkdir(parents=True)
    bundled_file = bundled_dir / "human_mouse_hcop_fifteen_column.txt"
    bundled_file.write_text(
        "human_symbol\tmouse_symbol\tsupport\n" "A\tA_mouse\tEnsembl,HGNC,MGI\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(_hcop_orthologs, "HCOP_BUNDLED_DIR", bundled_dir)

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


def test_hcop_orthologs_missing_version_lists_versions():
    with pytest.raises(
        FileNotFoundError,
        match="Available versions: bundled, latest",
    ):
        hcop.orthologs(target_organism="mouse", version="missing_hcop.tsv")


def test_hcop_orthologs_missing_bundled_file_fails_clearly(tmp_path, monkeypatch):
    monkeypatch.setattr(_hcop_orthologs, "HCOP_BUNDLED_DIR", tmp_path)
    with pytest.raises(
        FileNotFoundError,
        match="bundled HCOP file not found",
    ):
        hcop.orthologs(target_organism="mouse")
