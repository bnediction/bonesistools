#!/usr/bin/env python

"""
Utilities for loading GEO datasets as AnnData objects.
"""

from __future__ import annotations

import gzip
import os
import re
import sys
import time
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple, Union, cast
from urllib.error import URLError
from urllib.parse import unquote, urljoin
from urllib.request import urlopen

import pandas as pd
import scipy.io
from anndata import AnnData
from anndata import utils as ad_utils
from scipy import sparse

from ..._validation import _as_boolean

PathInput = Union[str, os.PathLike]

_GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
_GSM_PATTERN = re.compile(r"^GSM[0-9]+$")
_DOWNLOAD_CHUNK_SIZE = 1024 * 1024
_SUPPLEMENTARY_SUFFIXES = (
    "_matrix.mtx.gz",
    "_barcodes.tsv.gz",
    "_genes.tsv.gz",
    "_features.tsv.gz",
)


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: List[str] = []

    def handle_starttag(
        self,
        tag: str,
        attrs: List[Tuple[str, Optional[str]]],
    ) -> None:
        if tag != "a":
            return

        for name, value in attrs:
            if name == "href" and value is not None:
                self.links.append(unquote(value))


def from_geo(
    accession: str,
    cache_dir: Optional[PathInput] = None,
    quiet: bool = False,
) -> AnnData:
    """
    Download a GEO dataset and return it as an AnnData object.

    Parameters
    ----------
    accession: str
        GEO accession, such as `"GSM4138110"`.
    cache_dir: path-like, optional
        Directory used to cache downloaded files. If None, use
        `~/.cache/bonesistools/geo`.
    quiet: bool (default: False)
        Whether to suppress download progress messages.

    Returns
    -------
    AnnData
        Annotated count matrix loaded from the GEO supplementary files.

        GEO metadata is stored in:

        - `adata.uns["geo"]`: accession, accession type and source URL.

    Notes
    -----
    Currently only GSM supplementary datasets following the standard 10x
    matrix/barcodes/features layout are supported.

    Examples
    --------
    >>> adata = bt.sct.datasets.from_geo("GSM4138110")
    """

    if not isinstance(accession, str):
        raise TypeError(
            f"unsupported argument type for 'accession': "
            f"expected {str} but received {type(accession)}"
        )

    quiet = _as_boolean(quiet, "quiet")
    if accession.startswith("GSM"):
        return _from_gsm(
            accession,
            cache_dir=cache_dir,
            quiet=quiet,
        )

    raise ValueError(
        f"invalid argument value for 'accession': "
        f"currently expected a GSM accession but received {accession!r}"
    )


def _from_gsm(
    gsm: str,
    cache_dir: Optional[PathInput] = None,
    geo_ftp_base: str = _GEO_FTP_BASE,
    quiet: bool = False,
) -> AnnData:

    gsm = _as_gsm(gsm)
    quiet = _as_boolean(quiet, "quiet")
    sample_dir = _geo_cache_dir(gsm, cache_dir)
    matrix_path = sample_dir / "matrix.mtx.gz"
    barcodes_path = sample_dir / "barcodes.tsv.gz"
    feature_path = sample_dir / "features.tsv.gz"
    source_url = _gsm_base_url(gsm, geo_ftp_base=geo_ftp_base)

    if not (matrix_path.exists() and barcodes_path.exists() and feature_path.exists()):
        urls = _gsm_urls(gsm, geo_ftp_base=geo_ftp_base)
        source_url = urls["source"]
        _download(
            urls["matrix"],
            matrix_path,
            quiet=quiet,
        )
        _download(
            urls["barcodes"],
            barcodes_path,
            quiet=quiet,
        )
        _download(
            urls["features"],
            feature_path,
            quiet=quiet,
        )

    adata = _read_10x_mtx(matrix_path, barcodes_path, feature_path)
    adata.uns["geo"] = {
        "accession": gsm,
        "type": "GSM",
        "source": source_url,
    }

    return adata


def _as_gsm(gsm: str) -> str:

    if not isinstance(gsm, str):
        raise TypeError(
            f"unsupported argument type for 'gsm': "
            f"expected {str} but received {type(gsm)}"
        )

    if _GSM_PATTERN.match(gsm) is None:
        raise ValueError(
            f"invalid argument value for 'gsm': "
            f"expected a GSM accession but received {gsm!r}"
        )

    return gsm


def _gsm_series_dir(gsm: str) -> str:

    gsm = _as_gsm(gsm)
    digits = gsm[3:]
    if len(digits) <= 3:
        return "GSMnnn"
    return f"GSM{digits[:-3]}nnn"


def _gsm_base_url(gsm: str, geo_ftp_base: str = _GEO_FTP_BASE) -> str:

    gsm = _as_gsm(gsm)
    return f"{geo_ftp_base.rstrip('/')}/samples/" f"{_gsm_series_dir(gsm)}/{gsm}/suppl/"


def _gsm_urls(
    gsm: str,
    geo_ftp_base: str = _GEO_FTP_BASE,
) -> Dict[str, str]:

    gsm = _as_gsm(gsm)
    base_url = _gsm_base_url(gsm, geo_ftp_base=geo_ftp_base)
    parser = _LinkParser()
    try:
        response = urlopen(base_url)
    except (IsADirectoryError, URLError):
        response = urlopen(urljoin(base_url, "index.html"))

    with response:
        parser.feed(response.read().decode("utf-8", errors="replace"))

    urls = []
    for link in parser.links:
        filename = link.rsplit("/", 1)[-1]
        if not filename.startswith(gsm):
            continue
        if filename.endswith(_SUPPLEMENTARY_SUFFIXES):
            urls.append(urljoin(base_url, link))

    matrix_urls = _select_urls(urls, "_matrix.mtx.gz")
    barcode_urls = _select_urls(urls, "_barcodes.tsv.gz")
    gene_urls = _select_urls(urls, "_genes.tsv.gz")
    feature_urls = _select_urls(urls, "_features.tsv.gz")

    matrix_url = _require_single_url(matrix_urls, "matrix", gsm, base_url)
    barcodes_url = _require_single_url(barcode_urls, "barcode", gsm, base_url)
    feature_url = _select_gene_or_feature_url(gene_urls, feature_urls, gsm, base_url)

    return {
        "matrix": matrix_url,
        "barcodes": barcodes_url,
        "features": feature_url,
        "source": base_url,
    }


def _download(
    url: str,
    output: Path,
    quiet: bool = False,
) -> None:

    if output.exists():
        return

    quiet = _as_boolean(quiet, "quiet")
    progress_handle = None if quiet else _open_progress_handle()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = Path(f"{output}.tmp")
    try:
        with urlopen(url) as response:
            total_header = response.headers.get("Content-Length")
            total = int(total_header) if total_header is not None else None
            downloaded = 0
            start_time = time.monotonic()
            last_update = 0.0

            with temporary_output.open("wb") as handle:
                while True:
                    chunk = response.read(_DOWNLOAD_CHUNK_SIZE)
                    if not chunk:
                        break

                    handle.write(chunk)
                    downloaded += len(chunk)

                    if progress_handle is None:
                        continue

                    now = time.monotonic()
                    if now - last_update < 0.1:
                        continue

                    last_update = now
                    _print_download_progress(
                        output.name,
                        downloaded,
                        total,
                        now - start_time,
                        file=progress_handle,
                    )

        if not quiet:
            if progress_handle is not None:
                _print_download_progress(
                    output.name,
                    downloaded,
                    total,
                    time.monotonic() - start_time,
                    file=progress_handle,
                )
                print(file=progress_handle)

        temporary_output.replace(output)
    finally:
        if progress_handle is not None and progress_handle is not sys.stderr:
            progress_handle.close()


def _open_progress_handle() -> Optional[TextIO]:

    if Path("/dev/tty").exists():
        try:
            return open("/dev/tty", "w")
        except OSError:
            pass
    if sys.stderr.isatty():
        return sys.stderr
    return None


def _print_download_progress(
    filename: str,
    downloaded: int,
    total: Optional[int],
    elapsed: float,
    file: TextIO,
) -> None:

    speed = downloaded / elapsed if elapsed > 0 else 0
    if total is None or total == 0:
        print(
            f"\rdownloading {filename} {_format_bytes(downloaded):>8} "
            f"{_format_bytes(speed)}/s",
            end="",
            file=file,
            flush=True,
        )
        return

    percent = min(100, int(downloaded * 100 / total))
    width = 24
    filled = percent * width // 100
    if filled >= width:
        bar = "=" * width
    else:
        bar = "=" * filled + ">" + "." * max(0, width - filled - 1)
    print(
        f"\rdownloading {filename} {percent:3d}%[{bar}] "
        f"{_format_bytes(downloaded):>8} {_format_bytes(speed)}/s "
        f"in {_format_duration(elapsed)}",
        end="",
        file=file,
        flush=True,
    )


def _format_bytes(size: float) -> str:

    units = ("B", "K", "M", "G", "T")
    value = float(size)
    for unit in units:
        if abs(value) < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{value:.0f}{unit}"
            return f"{value:.2f}{unit}"
        value /= 1024

    return f"{value:.2f}T"


def _format_duration(seconds: float) -> str:

    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, remaining_seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{remaining_seconds:04.1f}s"
    hours, remaining_minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(remaining_minutes):02d}m{remaining_seconds:04.1f}s"


def _geo_cache_dir(accession: str, cache_dir: Optional[PathInput]) -> Path:

    if cache_dir is None:
        root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
        cache = root / "bonesistools" / "geo"
    else:
        cache = Path(cache_dir)
    return cache / accession


def _select_urls(urls: List[str], suffix: str) -> List[str]:

    return [url for url in urls if url.rsplit("/", 1)[-1].endswith(suffix)]


def _require_single_url(
    urls: List[str],
    kind: str,
    gsm: str,
    base_url: str,
) -> str:

    if len(urls) == 0:
        raise ValueError(f"{kind} file not found for {gsm}: {base_url}")
    if len(urls) > 1:
        raise ValueError(f"multiple {kind} files found for {gsm}")
    return urls[0]


def _select_gene_or_feature_url(
    gene_urls: List[str],
    feature_urls: List[str],
    gsm: str,
    base_url: str,
) -> str:

    if len(feature_urls) > 1:
        raise ValueError(f"multiple feature files found for {gsm}")
    if len(gene_urls) > 1:
        raise ValueError(f"multiple gene files found for {gsm}")
    if len(feature_urls) == 1:
        return feature_urls[0]
    if len(gene_urls) == 1:
        return gene_urls[0]
    raise ValueError(f"gene/feature file not found for {gsm}: {base_url}")


def _read_10x_mtx(
    matrix_path: Path,
    barcodes_path: Path,
    feature_path: Path,
) -> AnnData:

    counts = sparse.csr_matrix(scipy.io.mmread(str(matrix_path)))
    obs_names = pd.Index(
        _read_table(barcodes_path).iloc[:, 0].astype(str),
        name=None,
    )
    obs_names = ad_utils.make_index_unique(obs_names)
    obs_names.name = None

    var = _make_var(_read_table(feature_path))
    var.index = ad_utils.make_index_unique(var.index)
    var.index.name = None
    if counts.shape == (len(var), len(obs_names)):
        counts = counts.T.tocsr()
    elif counts.shape != (len(obs_names), len(var)):
        raise ValueError(
            "matrix shape is incompatible with barcodes and genes "
            f"(matrix={counts.shape}, cells={len(obs_names)}, genes={len(var)})"
        )

    adata = AnnData(
        X=counts,
        obs=pd.DataFrame(index=obs_names),
        var=var,
    )
    adata.obs.index.name = None
    adata.var.index.name = None

    return adata


def _read_table(path: Path) -> pd.DataFrame:

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t", header=None, dtype=str)


def _make_var(features: pd.DataFrame) -> pd.DataFrame:

    symbols = (
        features.iloc[:, 1] if features.shape[1] > 1 else features.iloc[:, 0]
    ).astype(str)
    feature_ids = features.iloc[:, 0].astype(str)
    var_names = pd.Index(symbols, name=None)
    var = pd.DataFrame(index=var_names)
    var["Accession"] = feature_ids.to_numpy()
    var["symbol"] = symbols.to_numpy()
    if features.shape[1] > 2:
        var["feature_type"] = (
            cast(pd.Series, features.iloc[:, 2]).astype(str).to_numpy()
        )

    return var
