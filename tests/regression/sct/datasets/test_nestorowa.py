#!/usr/bin/env python

import gzip
from pathlib import Path
from typing import Any, Dict, List, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools import sct
from bonesistools.sctools.datasets import _nestorowa as nestorowa_module
from bonesistools.sctools.datasets import _registry as registry_module


def _write_gzip(path: Path, text: str) -> None:

    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _all_cell_type_columns() -> List[str]:

    columns = []
    for group_columns in nestorowa_module._LABEL_GROUPS.values():
        columns.extend(group_columns)

    return columns


def _write_nestorowa_fixture(directory: Path) -> None:

    directory.mkdir(parents=True, exist_ok=True)
    _write_gzip(
        directory / "GSE81682_HTSeq_counts.txt.gz",
        "\n".join(
            [
                "gene\tc1\tc2\tc3",
                "ENSG1\t1\t0\t3",
                "ENSG2\t0\t2\t4",
                "NON1\t9\t9\t9",
                "",
            ]
        ),
    )
    cell_types = pd.DataFrame(
        0,
        index=["c1", "c2", "c3"],
        columns=_all_cell_type_columns(),
    )
    cell_types.loc["c1", "LTHSC"] = 1
    cell_types.loc["c2", "CMP"] = 1
    cell_types.loc["c2", "GMP"] = 1
    cell_types.to_csv(directory / "all_cell_types.txt", sep="\t")
    (directory / "cluster_ids.txt").write_text("c1 purple\nc2 gold\nc3 purple\n")


def _patch_gene_mapping(monkeypatch) -> None:

    def fake_map_ensembl_to_official_names(adata: ad.AnnData) -> None:

        adata.var_names = ["GeneA", "GeneB"]
        adata.var["symbol"] = adata.var_names.astype(str)

    monkeypatch.setattr(
        nestorowa_module,
        "_map_ensembl_to_official_names",
        fake_map_ensembl_to_official_names,
    )


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)

    return cast(csr_matrix, matrix)


def test_nestorowa_reads_source_files(monkeypatch, tmp_path):
    source_dir = tmp_path / "source"
    _write_nestorowa_fixture(source_dir)
    _patch_gene_mapping(monkeypatch)

    adata = nestorowa_module._read_full_nestorowa(
        nestorowa_module._nestorowa_files(source_dir)
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.shape == (3, 2)
    assert adata.obs_names.tolist() == ["c1", "c2", "c3"]
    assert adata.var_names.tolist() == ["GeneA", "GeneB"]
    assert adata.var["ensembl"].tolist() == ["ENSG1", "ENSG2"]
    assert adata.obs["label"].astype(str).tolist() == ["HSC", "CMP", "HSC"]
    assert adata.obs["clusters"].astype(int).tolist() == [1, 0, 1]
    assert adata.uns["nestorowa"]["labels"] == {
        "unique": 1,
        "multiple": 1,
        "missing": 1,
    }

    np.testing.assert_array_equal(
        _as_csr(adata.X).toarray(),
        np.array(
            [
                [1, 0],
                [0, 2],
                [3, 4],
            ],
        ),
    )
    assert len(adata.layers) == 0


def test_datasets_load_nestorowa_uses_cache(monkeypatch, tmp_path):
    _patch_gene_mapping(monkeypatch)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache-home"))
    monkeypatch.setitem(
        registry_module._DATASETS,
        "nestorowa",
        {
            **registry_module._DATASETS["nestorowa"],
            "cells": 3,
            "genes": 2,
        },
    )
    fixture_dir = tmp_path / "fixtures"
    _write_nestorowa_fixture(fixture_dir)
    downloaded: List[str] = []
    download_directories: List[Path] = []

    def fake_download(url: str, output: Path, quiet: bool = False) -> None:

        source = fixture_dir / output.name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(source.read_bytes())
        downloaded.append(url)
        download_directories.append(output.parent)

    monkeypatch.setattr(nestorowa_module, "_download", fake_download)

    first = bt.sct.datasets.load("nestorowa", quiet=True)
    second = bt.sct.datasets.load("nestorowa", quiet=True)

    assert first.shape == (3, 2)
    assert second.shape == (3, 2)
    assert first.obs["label"].astype(str).tolist() == ["HSC", "CMP", "HSC"]
    np.testing.assert_array_equal(
        _as_csr(first.X).toarray(),
        _as_csr(second.X).toarray(),
    )
    expected_urls = [url for _, url in nestorowa_module._NESTOROWA_URLS.values()]
    assert downloaded == expected_urls
    assert set(download_directories) == {
        registry_module._dataset_cache_dir("nestorowa")
    }
    assert registry_module._dataset_cache_dir("nestorowa").exists()


def test_nestorowa_deprecated_alias_calls_registry(monkeypatch, tmp_path):
    _patch_gene_mapping(monkeypatch)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache-home"))
    monkeypatch.setitem(
        registry_module._DATASETS,
        "nestorowa",
        {
            **registry_module._DATASETS["nestorowa"],
            "cells": 3,
            "genes": 2,
        },
    )
    fixture_dir = tmp_path / "fixtures"
    _write_nestorowa_fixture(fixture_dir)

    def fake_download(url: str, output: Path, quiet: bool = False) -> None:

        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes((fixture_dir / output.name).read_bytes())

    monkeypatch.setattr(nestorowa_module, "_download", fake_download)

    with pytest.warns(DeprecationWarning, match="datasets.load"):
        adata = bt.sct.datasets.nestorowa(quiet=True)
    loaded = bt.sct.datasets.load("nestorowa", quiet=True)

    assert adata.shape == (3, 2)
    assert (
        adata.obs["label"].astype(str).tolist()
        == loaded.obs["label"].astype(str).tolist()
    )
    np.testing.assert_array_equal(
        _as_csr(adata.X).toarray(),
        _as_csr(loaded.X).toarray(),
    )


def test_datasets_load_unknown_lists_available_datasets():
    with pytest.raises(ValueError, match="nestorowa.*pbmc3k|pbmc3k.*nestorowa"):
        bt.sct.datasets.load("unknown")


def test_nestorowa_gene_mapping_merges_duplicated_symbols(monkeypatch):
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 2.0, 3.0],
                [5.0, 6.0, 7.0],
            ],
        ),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["ENSG1", "ENSG2", "ENSG3"]),
    )
    adata.var["ensembl"] = adata.var_names.astype(str)

    def fake_convert_gene_identifiers(scdata, **_):

        scdata.var_names = ["GeneA", "GeneA", "GeneB"]
        return None

    def fake_standardize_gene_identifiers(scdata, **_):

        return None

    monkeypatch.setattr(
        sct.pp,
        "convert_gene_identifiers",
        fake_convert_gene_identifiers,
    )
    monkeypatch.setattr(
        sct.pp,
        "standardize_gene_identifiers",
        fake_standardize_gene_identifiers,
    )

    nestorowa_module._map_ensembl_to_official_names(adata)

    assert adata.var_names.tolist() == ["GeneA", "GeneB"]
    assert adata.var["symbol"].tolist() == ["GeneA", "GeneB"]
    assert adata.var["ensembl"].tolist() == ["ENSG1", "ENSG3"]
    np.testing.assert_array_equal(
        np.asarray(adata.X),
        np.array(
            [
                [3.0, 3.0],
                [11.0, 7.0],
            ],
        ),
    )
    assert len(adata.layers) == 0


def test_nestorowa_file_paths_use_given_directory(tmp_path):
    paths: Dict[str, Path] = nestorowa_module._nestorowa_files(tmp_path)

    assert set(paths) == {
        "read_counts",
        "cell_types",
        "cluster_ids",
    }
    assert all(path.parent == tmp_path for path in paths.values())
