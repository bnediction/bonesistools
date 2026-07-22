#!/usr/bin/env python

from __future__ import annotations

import gzip
import hashlib
import os
import pickle
import tempfile
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, Mapping

from ..cache import _unlink
from ..cache import path as cache_path

_MAPPING_CACHE_SCHEMA = 1
_PICKLE_PROTOCOL = 4
_REQUIRED_MAPPING_KEYS = frozenset(
    {
        "databases",
        "ensembl_id",
        "gene_id",
        "name",
    }
)


def _cached_identifier_mappings(
    source_file: Path,
    *,
    source_key: str,
    factory: Callable[[], Dict[str, Any]],
) -> Dict[str, Any]:
    """Return mappings derived from an NCBI source, rebuilding when needed."""

    source_digest = _file_digest(source_file)
    cache_file = _mapping_cache_file(source_key)

    try:
        return _read_mappings(
            cache_file,
            source_key=source_key,
            source_digest=source_digest,
        )
    except FileNotFoundError:
        pass
    except Exception as error:
        _discard_invalid_cache(cache_file, error)

    mappings = factory()
    try:
        _write_mappings(
            cache_file,
            source_key=source_key,
            source_digest=source_digest,
            mappings=mappings,
        )
    except Exception as error:
        warnings.warn(
            f"unable to store NCBI identifier cache at {str(cache_file)!r}: "
            f"{error}",
            UserWarning,
            stacklevel=2,
        )
    return mappings


def _mapping_cache_file(source_key: str) -> Path:

    digest = hashlib.sha256(source_key.encode("utf-8")).hexdigest()
    return cache_path("ncbi") / "mappings" / f"{digest}.pickle.gz"


def _file_digest(file: Path) -> str:

    digest = hashlib.sha256()
    with file.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_mappings(
    cache_file: Path,
    *,
    source_key: str,
    source_digest: str,
) -> Dict[str, Any]:

    with gzip.open(cache_file, "rb") as handle:
        payload = pickle.load(handle)

    if not isinstance(payload, Mapping):
        raise ValueError("invalid cache payload")
    if payload.get("schema") != _MAPPING_CACHE_SCHEMA:
        raise ValueError("unsupported cache schema")
    if payload.get("source_key") != source_key:
        raise ValueError("cache source mismatch")
    if payload.get("source_digest") != source_digest:
        raise ValueError("cached NCBI source is obsolete")

    mappings = payload.get("mappings")
    if not isinstance(mappings, dict):
        raise ValueError("invalid cached mappings")
    if not _REQUIRED_MAPPING_KEYS <= mappings.keys():
        raise ValueError("incomplete cached mappings")
    if not all(isinstance(mappings[key], dict) for key in _REQUIRED_MAPPING_KEYS):
        raise ValueError("invalid cached mapping indexes")
    return mappings


def _write_mappings(
    cache_file: Path,
    *,
    source_key: str,
    source_digest: str,
    mappings: Dict[str, Any],
) -> None:

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{cache_file.name}.",
        suffix=".tmp",
        dir=str(cache_file.parent),
    )
    temporary_file = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as raw_handle:
            with gzip.GzipFile(
                fileobj=raw_handle,
                mode="wb",
                compresslevel=1,
                mtime=0,
            ) as compressed_handle:
                pickle.dump(
                    {
                        "schema": _MAPPING_CACHE_SCHEMA,
                        "source_key": source_key,
                        "source_digest": source_digest,
                        "mappings": mappings,
                    },
                    compressed_handle,
                    protocol=_PICKLE_PROTOCOL,
                )
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        temporary_file.replace(cache_file)
    finally:
        _unlink(temporary_file)


def _discard_invalid_cache(cache_file: Path, error: Exception) -> None:

    _unlink(cache_file)
    warnings.warn(
        "removing incompatible, corrupt, or obsolete NCBI identifier cache "
        f"at {str(cache_file)!r}: {error}; mappings will be rebuilt",
        UserWarning,
        stacklevel=3,
    )
