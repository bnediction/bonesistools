#!/usr/bin/env python

"""
Utilities for loading the PBMC3k tutorial dataset.
"""

from __future__ import annotations

import os
import shutil
import tarfile
from pathlib import Path, PurePosixPath
from typing import Dict, List

from anndata import AnnData

from ..._validation import _as_boolean
from ._geo import _download, _read_10x_mtx

_PBMC3K_URL = (
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
_PBMC3K_ARCHIVE = "pbmc3k_filtered_gene_bc_matrices.tar.gz"
_PBMC3K_10X_DIR = PurePosixPath("filtered_gene_bc_matrices") / "hg19"


def _load_pbmc3k(
    cache_dir: Path,
    quiet: bool = False,
) -> AnnData:
    """
    Load raw PBMC3k counts from the 10x Genomics tutorial dataset.
    """

    quiet = _as_boolean(quiet, "quiet")
    paths = _pbmc3k_files(cache_dir)
    _download(_PBMC3K_URL, paths["archive"], quiet=quiet)
    _extract_pbmc3k_archive(paths["archive"], paths["directory"])

    adata = _read_10x_mtx(
        paths["matrix"],
        paths["barcodes"],
        paths["features"],
    )
    adata.uns["pbmc3k"] = {
        "source": {
            "counts": "10x Genomics PBMC3k",
            "url": _PBMC3K_URL,
        },
        "license": "CC BY 4.0",
    }

    return adata


def _pbmc3k_files(directory: Path) -> Dict[str, Path]:

    data_dir = directory / "filtered_gene_bc_matrices" / "hg19"

    return {
        "archive": directory / _PBMC3K_ARCHIVE,
        "directory": directory,
        "matrix": data_dir / "matrix.mtx",
        "barcodes": data_dir / "barcodes.tsv",
        "features": data_dir / "genes.tsv",
    }


def _extract_pbmc3k_archive(
    archive_path: Path,
    output_dir: Path,
) -> None:

    paths = _pbmc3k_files(output_dir)
    if (
        paths["matrix"].exists()
        and paths["barcodes"].exists()
        and paths["features"].exists()
    ):
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    expected_members = [
        str(_PBMC3K_10X_DIR / "matrix.mtx"),
        str(_PBMC3K_10X_DIR / "barcodes.tsv"),
        str(_PBMC3K_10X_DIR / "genes.tsv"),
    ]

    with tarfile.open(archive_path, "r:gz") as archive:
        members = [_get_archive_member(archive, name) for name in expected_members]
        _validate_archive_members(output_dir, members)
        for member in members:
            _extract_archive_member(archive, member, output_dir)


def _get_archive_member(
    archive: tarfile.TarFile,
    name: str,
) -> tarfile.TarInfo:

    try:
        return archive.getmember(name)
    except KeyError as exc:
        raise ValueError(f"PBMC3k archive is missing {name!r}") from exc


def _validate_archive_members(
    output_dir: Path,
    members: List[tarfile.TarInfo],
) -> None:

    output = output_dir.resolve()
    for member in members:
        target = (output_dir / member.name).resolve()
        if os.path.commonpath([str(output), str(target)]) != str(output):
            raise ValueError(f"unsafe PBMC3k archive member: {member.name!r}")


def _extract_archive_member(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
    output_dir: Path,
) -> None:

    if not member.isfile():
        raise ValueError(f"PBMC3k archive member is not a file: {member.name!r}")

    source = archive.extractfile(member)
    if source is None:
        raise ValueError(f"cannot extract PBMC3k archive member: {member.name!r}")

    target = output_dir / member.name
    target.parent.mkdir(parents=True, exist_ok=True)
    with source, target.open("wb") as output:
        shutil.copyfileobj(source, output)
