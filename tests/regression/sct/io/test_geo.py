#!/usr/bin/env python

import gzip
import io
import warnings
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

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


class FakeDownloadResponse:
    def __init__(self, chunks, headers=None):
        self.chunks = list(chunks)
        self.headers = {} if headers is None else headers

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def read(self, _size):
        if not self.chunks:
            return b""
        return self.chunks.pop(0)


class NonClosingStringIO(io.StringIO):
    def close(self):
        pass


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
        "gene_id_1\tGene1\n" "gene_id_2\tGene2\n" "gene_id_3\tGene3\n",
    )


def _as_csr(matrix: Any) -> csr_matrix:

    assert sparse.isspmatrix_csr(matrix)
    return cast(csr_matrix, matrix)


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

    with pytest.raises(TypeError):
        geo._as_gsm(cast(Any, 12345))

    with pytest.raises(TypeError):
        bt.sct.io.from_geo(cast(Any, 12345), quiet=True)

    with pytest.raises(ValueError):
        geo._as_gsm("GSE12345")

    with pytest.raises(ValueError):
        geo._as_gsm("foo")

    with pytest.raises(ValueError):
        bt.sct.io.from_geo("GSE12345", quiet=True)

    with pytest.raises(ValueError):
        bt.sct.io.from_geo("foo", quiet=True)


def test_geo_gsm_path_generation():
    assert geo._gsm_series_dir("GSM12") == "GSMnnn"
    assert geo._gsm_series_dir("GSM4138110") == "GSM4138nnn"
    assert geo._gsm_base_url("GSM4138110").endswith(
        "/samples/GSM4138nnn/GSM4138110/suppl/"
    )


def test_geo_link_parser_keeps_decoded_anchor_hrefs_only():
    parser = geo._LinkParser()

    parser.feed("""
        <span href="ignored">ignored</span>
        <a>missing</a>
        <a href="GSM1%20matrix.mtx.gz">matrix</a>
        """)

    assert parser.links == ["GSM1 matrix.mtx.gz"]


def test_geo_html_parsing_detects_10x_supplementary_files(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_sample_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_sample_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_sample_genes.tsv.gz">genes</a>
        <a href="GSM4138110_sample_features.tsv.gz">features</a>
        <a href="GSM9999999_sample_features.tsv.gz">wrong sample</a>
        <a href="GSM4138110_fragments.tsv.gz">ignored</a>
        """,
    )

    urls = geo._gsm_urls("GSM4138110")

    assert urls["matrix"].endswith("GSM4138110_sample_matrix.mtx.gz")
    assert urls["barcodes"].endswith("GSM4138110_sample_barcodes.tsv.gz")
    assert urls["features"].endswith("GSM4138110_sample_features.tsv.gz")
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


def test_geo_html_parsing_rejects_multiple_gene_or_feature_files(monkeypatch):
    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_a_features.tsv.gz">features</a>
        <a href="GSM4138110_b_features.tsv.gz">features</a>
        """,
    )

    with pytest.raises(ValueError, match="multiple feature files"):
        geo._gsm_urls("GSM4138110")

    _fake_urlopen_from_html(
        monkeypatch,
        """
        <a href="GSM4138110_matrix.mtx.gz">matrix</a>
        <a href="GSM4138110_barcodes.tsv.gz">barcodes</a>
        <a href="GSM4138110_a_genes.tsv.gz">genes</a>
        <a href="GSM4138110_b_genes.tsv.gz">genes</a>
        """,
    )

    with pytest.raises(ValueError, match="multiple gene files"):
        geo._gsm_urls("GSM4138110")


def test_geo_formatting_helpers_produce_human_readable_progress():
    assert geo._format_bytes(512) == "512B"
    assert geo._format_bytes(2048) == "2.00K"
    assert geo._format_duration(2.5) == "2.5s"
    assert geo._format_duration(62.5) == "1m02.5s"
    assert geo._format_duration(3723.0) == "1h02m03.0s"

    unknown_total = io.StringIO()
    geo._print_download_progress(
        "matrix.mtx.gz",
        downloaded=2048,
        total=None,
        elapsed=2.0,
        file=unknown_total,
    )
    assert "2.00K" in unknown_total.getvalue()
    assert "1.00K/s" in unknown_total.getvalue()

    known_total = io.StringIO()
    geo._print_download_progress(
        "matrix.mtx.gz",
        downloaded=50,
        total=100,
        elapsed=5.0,
        file=known_total,
    )
    assert " 50%" in known_total.getvalue()
    assert "[============>...........]" in known_total.getvalue()

    aligned = []
    for filename in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
        progress = io.StringIO()
        geo._print_download_progress(
            filename,
            downloaded=100,
            total=100,
            elapsed=1.0,
            file=progress,
        )
        aligned.append(progress.getvalue())

    assert len({progress.index("%[") for progress in aligned}) == 1
    assert len({progress.index(" in ") for progress in aligned}) == 1


def test_geo_cache_dir_uses_xdg_cache_home(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))

    assert geo._geo_cache_dir("GSM1", None) == (
        tmp_path / "xdg" / "bonesistools" / "geo" / "GSM1"
    )


def test_geo_download_skips_existing_output(monkeypatch, tmp_path):
    output = tmp_path / "matrix.mtx.gz"
    output.write_bytes(b"cached")

    def fail_urlopen(_):
        raise AssertionError("existing output should not be downloaded")

    monkeypatch.setattr(geo, "urlopen", fail_urlopen)

    geo._download("https://example.org/matrix.mtx.gz", output, quiet=True)

    assert output.read_bytes() == b"cached"


def test_geo_download_writes_temporary_file_and_reports_progress(monkeypatch, tmp_path):
    output = tmp_path / "matrix.mtx.gz"
    progress = NonClosingStringIO()
    times = iter([0.0, 0.05, 0.2, 0.4])

    def fake_urlopen(_request):
        return FakeDownloadResponse(
            [b"abc", b"def"],
            headers={"Content-Length": "6"},
        )

    monkeypatch.setattr(geo, "urlopen", fake_urlopen)
    monkeypatch.setattr(geo, "_open_progress_handle", lambda: progress)
    monkeypatch.setattr(geo.time, "monotonic", lambda: next(times))

    geo._download("https://example.org/matrix.mtx.gz", output, quiet=False)

    assert output.read_bytes() == b"abcdef"
    assert not Path(f"{output}.tmp").exists()
    assert "100%" in progress.getvalue()
    assert progress.getvalue().endswith("\n")


def test_geo_open_progress_handle_prefers_tty_and_falls_back_to_stderr(
    monkeypatch,
):
    class ExistingTTYPath:
        def __init__(self, path):
            self.path = path

        def exists(self):
            return True

    class MissingTTYPath:
        def __init__(self, path):
            self.path = path

        def exists(self):
            return False

    handle = NonClosingStringIO()
    monkeypatch.setattr(geo, "Path", ExistingTTYPath)
    monkeypatch.setattr(geo, "open", lambda *_args, **_kwargs: handle, raising=False)

    assert geo._open_progress_handle() is handle

    monkeypatch.setattr(
        geo,
        "open",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError()),
        raising=False,
    )
    monkeypatch.setattr(geo.sys.stderr, "isatty", lambda: True)

    assert geo._open_progress_handle() is geo.sys.stderr

    monkeypatch.setattr(geo, "Path", MissingTTYPath)
    monkeypatch.setattr(geo.sys.stderr, "isatty", lambda: False)

    assert geo._open_progress_handle() is None


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
    _write_gzip(sample_dir / "features.tsv.gz", "gene_id\tGene\n")

    def fail_urlopen(_):
        raise AssertionError("cache hit should not open GEO URLs")

    monkeypatch.setattr(geo, "urlopen", fail_urlopen)

    adata = bt.sct.io.from_geo(
        "GSM5492245",
        cache_dir=tmp_path / "cache",
        quiet=True,
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.shape == (1, 1)
    assert _as_csr(adata.X).toarray().tolist() == [[7]]


def test_geo_from_geo_makes_duplicate_var_names_unique_before_anndata_warning(
    tmp_path,
):
    sample_dir = tmp_path / "cache" / "GSM5492245"
    _write_gzip(
        sample_dir / "matrix.mtx.gz",
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "2 1 2",
                "1 1 7",
                "2 1 11",
                "",
            ]
        ),
    )
    _write_gzip(sample_dir / "barcodes.tsv.gz", "cell\n")
    _write_gzip(
        sample_dir / "features.tsv.gz",
        "gene_id_1\tGene\n" "gene_id_2\tGene\n",
    )

    with warnings.catch_warnings(record=True) as records:
        warnings.simplefilter("always")
        adata = bt.sct.io.from_geo(
            "GSM5492245",
            cache_dir=tmp_path / "cache",
            quiet=True,
        )

    assert not any(
        "Variable names are not unique" in str(record.message) for record in records
    )
    assert adata.var_names.tolist() == ["Gene", "Gene-1"]
    assert adata.var["symbol"].tolist() == ["Gene", "Gene"]
    np.testing.assert_array_equal(_as_csr(adata.X).toarray(), np.array([[7, 11]]))


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
        _as_csr(adata.X).toarray(),
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
    assert (tmp_path / "cache" / "GSM5492245" / "features.tsv.gz").exists()


def test_geo_read_10x_mtx_rejects_incompatible_shapes(tmp_path):
    sample_dir = tmp_path / "cache" / "GSM5492245"
    _write_gzip(
        sample_dir / "matrix.mtx.gz",
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "3 3 1",
                "1 1 1",
                "",
            ]
        ),
    )
    _write_gzip(sample_dir / "barcodes.tsv.gz", "cell1\ncell2\n")
    _write_gzip(sample_dir / "features.tsv.gz", "gene_id\tGene\n")

    with pytest.raises(ValueError, match="matrix shape is incompatible"):
        geo._read_10x_mtx(
            sample_dir / "matrix.mtx.gz",
            sample_dir / "barcodes.tsv.gz",
            sample_dir / "features.tsv.gz",
        )


def test_geo_read_10x_mtx_preserves_feature_type_column(tmp_path):
    sample_dir = tmp_path / "cache" / "GSM5492245"
    _write_gzip(
        sample_dir / "matrix.mtx.gz",
        "\n".join(
            [
                "%%MatrixMarket matrix coordinate integer general",
                "%",
                "1 1 1",
                "1 1 1",
                "",
            ]
        ),
    )
    _write_gzip(sample_dir / "barcodes.tsv.gz", "cell\n")
    _write_gzip(sample_dir / "features.tsv.gz", "gene_id\tGene\tGene Expression\n")

    adata = geo._read_10x_mtx(
        sample_dir / "matrix.mtx.gz",
        sample_dir / "barcodes.tsv.gz",
        sample_dir / "features.tsv.gz",
    )

    assert adata.var["feature_type"].tolist() == ["Gene Expression"]


def test_geo_read_10x_mtx_supports_scipy_without_spmatrix_keyword(
    monkeypatch,
    tmp_path,
):
    sample_dir = tmp_path / "cache" / "GSM5492245"
    sample_dir.mkdir(parents=True)
    _write_gzip(sample_dir / "barcodes.tsv.gz", "cell\n")
    _write_gzip(sample_dir / "features.tsv.gz", "gene_id\tGene\n")
    matrix_path = sample_dir / "matrix.mtx.gz"
    matrix_path.write_text("unused", encoding="utf-8")
    calls = []

    def fake_mmread(path, **kwargs):
        calls.append((path, kwargs))
        if "spmatrix" in kwargs:
            raise TypeError("old SciPy")
        return csr_matrix([[7]])

    monkeypatch.setattr(geo.scipy.io, "mmread", fake_mmread)

    adata = geo._read_10x_mtx(
        matrix_path,
        sample_dir / "barcodes.tsv.gz",
        sample_dir / "features.tsv.gz",
    )

    assert calls == [
        (str(matrix_path), {"spmatrix": True}),
        (str(matrix_path), {}),
    ]
    np.testing.assert_array_equal(_as_csr(adata.X).toarray(), np.array([[7]]))
