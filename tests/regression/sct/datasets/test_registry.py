#!/usr/bin/env python

import tarfile
from pathlib import Path
from typing import Any, List, cast

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.datasets import _pbmc3k as pbmc3k_module
from bonesistools.sctools.datasets import _registry as registry_module

_REQUIRED_METADATA_FIELDS = {
    "description",
    "organism",
    "tissue",
    "technology",
    "source",
    "license",
    "url",
    "citation",
}


def _write_pbmc3k_fixture(directory: Path) -> None:

    data_dir = directory / "filtered_gene_bc_matrices" / "hg19"
    data_dir.mkdir(parents=True)
    (data_dir / "matrix.mtx").write_text(
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "3 2 4",
                "1 1 1",
                "2 1 2",
                "1 2 3",
                "3 2 4",
                "",
            ]
        )
    )
    (data_dir / "barcodes.tsv").write_text("cell1\ncell2\n")
    (data_dir / "genes.tsv").write_text(
        "gene_id_1\tGeneA\n" "gene_id_2\tGeneB\n" "gene_id_3\tGeneC\n"
    )
    with tarfile.open(directory / "pbmc3k.tar.gz", "w:gz") as archive:
        archive.add(
            directory / "filtered_gene_bc_matrices",
            arcname="filtered_gene_bc_matrices",
        )


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)

    return cast(csr_matrix, matrix)


def test_datasets_load_pbmc3k_returns_raw_counts(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache-home"))
    monkeypatch.setitem(
        registry_module._DATASETS,
        "pbmc3k",
        {
            **registry_module._DATASETS["pbmc3k"],
            "cells": 2,
            "genes": 3,
        },
    )
    fixture_dir = tmp_path / "fixtures"
    fixture_dir.mkdir()
    _write_pbmc3k_fixture(fixture_dir)
    downloaded: List[str] = []

    def fake_download(url: str, output: Path, quiet: bool = False) -> None:

        if output.exists():
            return

        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes((fixture_dir / "pbmc3k.tar.gz").read_bytes())
        downloaded.append(url)

    monkeypatch.setattr(pbmc3k_module, "_download", fake_download)

    first = bt.sct.datasets.load("pbmc3k", quiet=True)
    second = bt.sct.datasets.load("pbmc3k", quiet=True)

    assert first.shape == (2, 3)
    assert first.obs_names.tolist() == ["cell1", "cell2"]
    assert first.var_names.tolist() == ["GeneA", "GeneB", "GeneC"]
    assert _as_csr(first.X).nnz > 0
    assert len(first.obs) > 0
    assert len(first.var) > 0
    assert "cell_type" not in first.obs
    assert first.uns["pbmc3k"]["source"]["url"] == pbmc3k_module._PBMC3K_URL
    np.testing.assert_array_equal(
        _as_csr(first.X).toarray(),
        np.array([[1, 2, 0], [3, 0, 4]]),
    )
    np.testing.assert_array_equal(
        _as_csr(first.X).toarray(),
        _as_csr(second.X).toarray(),
    )
    assert downloaded == [pbmc3k_module._PBMC3K_URL]


def test_datasets_load_rejects_invalid_cached_shape(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache-home"))

    def fake_loader(_cache_dir: Path, _quiet: bool) -> ad.AnnData:

        return ad.AnnData(X=np.ones((1, 1)))

    monkeypatch.setitem(
        registry_module._DATASET_LOADERS,
        "pbmc3k",
        fake_loader,
    )

    with pytest.raises(
        ValueError,
        match=(
            "invalid cached dataset 'pbmc3k'.*"
            "expected 2700 observations but found 1.*"
            "expected 32738 features but found 1.*"
            "bt\\.sct\\.datasets\\.clear\\('pbmc3k'\\)"
        ),
    ):
        bt.sct.datasets.load("pbmc3k", quiet=True)


def test_datasets_info_pbmc3k_contains_source_license_url_and_citation():
    metadata = bt.sct.datasets.info("PBMC3K")

    assert metadata["name"] == "pbmc3k"
    assert metadata["organism"] == "Homo sapiens"
    assert metadata["tissue"] == "peripheral blood"
    assert metadata["technology"] == "10x Chromium Single Cell 3'"
    assert metadata["cells"] == 2700
    assert metadata["genes"] == 32738
    assert metadata["source"] == "10x Genomics"
    assert metadata["license"] == "CC BY 4.0"
    assert metadata["url"] == (
        "https://www.10xgenomics.com/datasets/"
        "3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0"
    )
    assert "Zheng" in str(metadata["citation"])
    assert "10.1038/ncomms14049" in str(metadata["citation"])
    assert len(str(metadata["description"])) > 0


def test_datasets_info_nestorowa_contains_source_license_and_counts_url():
    metadata = bt.sct.datasets.info("nestorowa")

    assert metadata["organism"] == "Mus musculus"
    assert metadata["tissue"] == "bone marrow"
    assert metadata["technology"] == "Smart-seq2"
    assert metadata["cells"] == 1656
    assert metadata["genes"] == 46078
    assert metadata["source"] == "NCBI GEO (GSE81682)"
    assert metadata["license"] == "Not specified"
    assert str(metadata["url"]).endswith("acc=GSE81682")


def test_datasets_available_lists_registered_datasets():
    datasets = bt.sct.datasets.available()

    assert "nestorowa" in datasets.index
    assert "pbmc3k" in datasets.index
    assert "description" in datasets.columns
    assert "organism" in datasets.columns
    assert "cells" in datasets.columns
    assert "genes" in datasets.columns
    assert "url" in datasets.columns


def test_dataset_metadata_entries_follow_documented_conventions():
    datasets = bt.sct.datasets.available()

    for name in datasets.index:
        metadata = bt.sct.datasets.info(str(name))

        assert _REQUIRED_METADATA_FIELDS <= set(metadata)
        assert isinstance(metadata["description"], str)
        assert isinstance(metadata["organism"], str)
        tissue = metadata["tissue"]
        assert isinstance(tissue, str)
        assert tissue == tissue.lower()
        assert isinstance(metadata["technology"], str)
        assert isinstance(metadata["source"], str)
        assert isinstance(metadata["license"], str)
        url = metadata["url"]
        assert isinstance(url, str)
        assert url.startswith(("http://", "https://"))
        assert isinstance(metadata["citation"], str)
        assert "cells" in metadata
        assert "genes" in metadata or "peaks" in metadata

    assert not str(bt.sct.datasets.info("pbmc3k")["url"]).endswith(".tar.gz")


def test_datasets_info_unknown_lists_available_datasets():
    with pytest.raises(ValueError, match="nestorowa.*pbmc3k|pbmc3k.*nestorowa"):
        bt.sct.datasets.info("unknown")


def test_datasets_clear_removes_named_and_full_cache(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache-home"))
    root = registry_module._dataset_cache_root()
    pbmc3k_dir = registry_module._dataset_cache_dir("pbmc3k")
    nestorowa_dir = registry_module._dataset_cache_dir("nestorowa")
    (pbmc3k_dir / "data").mkdir(parents=True)
    (pbmc3k_dir / "data" / "file.txt").write_text("pbmc")
    (nestorowa_dir / "data").mkdir(parents=True)
    (nestorowa_dir / "data" / "file.txt").write_text("nestorowa")

    bt.sct.datasets.clear("pbmc3k")
    bt.sct.datasets.clear("pbmc3k")

    assert not pbmc3k_dir.exists()
    assert nestorowa_dir.exists()

    bt.sct.datasets.clear()
    bt.sct.datasets.clear()

    assert not root.exists()


def test_datasets_clear_unknown_lists_available_datasets():
    with pytest.raises(ValueError, match="nestorowa.*pbmc3k|pbmc3k.*nestorowa"):
        bt.sct.datasets.clear("unknown")
