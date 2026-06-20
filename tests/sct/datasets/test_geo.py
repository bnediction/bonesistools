#!/usr/bin/env python

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

import bonesistools as bt
from bonesistools.sctools.datasets import _geo as geo


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def read(self):
        return self.payload


def _fake_urlopen_from_html(monkeypatch, html):
    opened = []

    def fake_urlopen(url):
        opened.append(url)
        return FakeResponse(html.encode())

    monkeypatch.setattr(geo, "urlopen", fake_urlopen)
    return opened


def _write_gzip(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _write_10x_files(directory: Path) -> None:
    _write_gzip(
        directory / "GSM5492245_RNA_matrix.mtx.gz",
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "3 2 3",
                "1 1 1",
                "2 1 2",
                "3 2 3",
                "",
            ]
        ),
    )
    _write_gzip(
        directory / "GSM5492245_RNA_barcodes.tsv.gz",
        "cell1\ncell2\n",
    )
    _write_gzip(
        directory / "GSM5492245_RNA_genes.tsv.gz",
        "gene_id_1\tGene1\n"
        "gene_id_2\tGene2\n"
        "gene_id_3\tGene3\n",
    )


def _write_geo_index(directory: Path) -> None:
    (directory / "index.html").write_text(
        "\n".join(
            [
                '<a href="GSM5492245_RNA_matrix.mtx.gz">matrix</a>',
                '<a href="GSM5492245_RNA_barcodes.tsv.gz">barcodes</a>',
                '<a href="GSM5492245_RNA_genes.tsv.gz">genes</a>',
                '<a href="GSM5492245_signal.bw">ignored</a>',
                '<a href="GSM5492245_fragments.tsv.gz">ignored</a>',
            ]
        )
    )


def test_geo_validates_gsm_accessions():
    assert geo._as_gsm("GSM12345") == "GSM12345"

    with pytest.raises(ValueError):
        geo._as_gsm("GSE12345")

    with pytest.raises(ValueError):
        geo._as_gsm("foo")

    with pytest.raises(ValueError):
        bt.sct.datasets.from_geo("GSE12345", quiet=True)

    with pytest.raises(ValueError):
        bt.sct.datasets.from_geo("foo", quiet=True)


def test_geo_gsm_path_generation():
    assert geo._gsm_series_dir("GSM4138110") == "GSM4138nnn"
    assert geo._gsm_base_url("GSM4138110").endswith(
        "/samples/GSM4138nnn/GSM4138110/suppl/"
    )


def test_geo_html_parsing_detects_10x_supplementary_files(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_sample_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_sample_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_sample_genes.tsv.gz">genes</a>
        <a href="GSM4138110_sample_features.tsv.gz">features</a>
        <a href="GSM4138110_fragments.tsv.gz">ignored</a>
        """,
    )

    urls = geo._gsm_urls("GSM4138110")

    assert urls["matrix"].endswith("GSM4138110_sample_matrix.mtx.gz")
    assert urls["barcodes"].endswith("GSM4138110_sample_barcodes.tsv.gz")
    assert urls["genes"].endswith("GSM4138110_sample_features.tsv.gz")
    assert urls["source"].endswith("/samples/GSM4138nnn/GSM4138110/suppl/")


def test_geo_html_parsing_rejects_multiple_matrix_files(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_a_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_b_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_genes.tsv.gz">genes</a>
        """,
    )

    with pytest.raises(ValueError, match="multiple matrix files"):
        geo._gsm_urls("GSM4138110")


def test_geo_html_parsing_rejects_missing_matrix_file(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_genes.tsv.gz">genes</a>
        """,
    )

    with pytest.raises(ValueError, match="matrix file not found"):
        geo._gsm_urls("GSM4138110")


def test_geo_html_parsing_rejects_missing_barcode_file(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_genes.tsv.gz">genes</a>
        """,
    )

    with pytest.raises(ValueError, match="barcode file not found"):
        geo._gsm_urls("GSM4138110")


def test_geo_html_parsing_rejects_missing_gene_or_feature_file(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_barcodes.tsv.gz">barcodes</a>
        """,
    )

    with pytest.raises(ValueError, match="gene/feature file not found"):
        geo._gsm_urls("GSM4138110")


def test_geo_from_geo_uses_cache_without_downloading(monkeypatch, tmp_path):
    sample_dir = tmp_path / "cache" / "GSM5492245"
    _write_gzip(
        sample_dir / "matrix.mtx.gz",
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "1 1 1",
                "1 1 7",
                "",
            ]
        ),
    )
    _write_gzip(sample_dir / "barcodes.tsv.gz", "cell\n")
    _write_gzip(sample_dir / "genes.tsv.gz", "gene_id\tGene\n")

    def fail_urlopen(_):
        raise AssertionError("cache hit should not open GEO URLs")

    monkeypatch.setattr(geo, "urlopen", fail_urlopen)

    adata = bt.sct.datasets.from_geo(
        "GSM5492245",
        cache_dir=tmp_path / "cache",
        quiet=True,
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.shape == (1, 1)
    assert adata.X.toarray().tolist() == [[7]]


def test_geo_from_gsm_downloads_and_reads_local_geo_10x_dataset(tmp_path):
    suppl = tmp_path / "geo" / "samples" / "GSM5492nnn" / "GSM5492245" / "suppl"
    suppl.mkdir(parents=True)
    _write_10x_files(suppl)
    _write_geo_index(suppl)

    adata = geo._from_gsm(
        "GSM5492245",
        cache_dir=tmp_path / "cache",
        geo_ftp_base=(tmp_path / "geo").as_uri(),
        quiet=True,
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.shape == (2, 3)
    assert adata.obs_names.tolist() == ["cell1", "cell2"]
    assert adata.var_names.tolist() == ["Gene1", "Gene2", "Gene3"]
    assert adata.var["Accession"].tolist() == ["gene_id_1", "gene_id_2", "gene_id_3"]
    np.testing.assert_array_equal(
        adata.X.toarray(),
        np.array([[1, 2, 0], [0, 0, 3]]),
    )
    assert adata.uns["geo"] == {
        "accession": "GSM5492245",
        "type": "GSM",
        "source": f"{(tmp_path / 'geo').as_uri()}/samples/"
        "GSM5492nnn/GSM5492245/suppl/",
    }
    assert (tmp_path / "cache" / "GSM5492245" / "matrix.mtx.gz").exists()
    assert (tmp_path / "cache" / "GSM5492245" / "barcodes.tsv.gz").exists()
    assert (tmp_path / "cache" / "GSM5492245" / "genes.tsv.gz").exists()
