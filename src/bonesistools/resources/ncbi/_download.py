#!/usr/bin/env python

from __future__ import annotations

import gzip
import shutil
import urllib.error
import urllib.request
from pathlib import Path

from ..cache import _cached_download

_LATEST_CACHE_MAX_AGE = 72 * 60 * 60


def _download_gene_info(
    url: str,
    outfile: Path,
    *,
    cache_latest: bool,
) -> None:
    """Download and decompress an NCBI gene_info file atomically."""

    outfile.parent.mkdir(parents=True, exist_ok=True)
    temporary_gzip = Path(f"{outfile}.gz.tmp")
    temporary_outfile = Path(f"{outfile}.tmp")

    try:
        if cache_latest:
            compressed_file = _cached_download(
                url,
                resource="ncbi",
                category="gene_info",
                max_age=_LATEST_CACHE_MAX_AGE,
                suffix=".gene_info.gz",
            )
        else:
            with urllib.request.urlopen(url) as response:
                with open(temporary_gzip, "wb") as file:
                    shutil.copyfileobj(response, file)
            compressed_file = temporary_gzip

        with gzip.open(compressed_file, "rb") as source:
            with open(temporary_outfile, "wb") as destination:
                shutil.copyfileobj(source, destination)

        temporary_outfile.replace(outfile)
    except (OSError, urllib.error.URLError) as error:
        raise RuntimeError("failed to download NCBI gene_info file") from error
    finally:
        _unlink(temporary_gzip)
        _unlink(temporary_outfile)


def _unlink(file: Path) -> None:

    try:
        file.unlink()
    except FileNotFoundError:
        pass
