#!/usr/bin/env python

"""
Management of locally cached biological resources.

Current remote resources are reused for up to 72 hours before revalidation.
Versioned archives remain cached until explicitly removed.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path
from typing import Dict, List, Mapping, Optional, cast
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from typing_extensions import Protocol

_CACHE_RESOURCES = ("hcop", "ncbi", "omnipath")
_DOWNLOAD_CHUNK_SIZE = 1024 * 1024
_DOWNLOAD_HEADERS = {"User-Agent": "bonesistools"}

__all__ = [
    "clear",
    "path",
]


class _DownloadResponse(Protocol):
    headers: Mapping[str, str]

    def read(self, size: int = -1) -> bytes: ...


def path(resource: Optional[str] = None) -> Path:
    """
    Return the local cache path for downloaded biological resources.

    Examples
    --------
    >>> import bonesistools as bt
    >>> cache_dir = bt.resources.cache.path("omnipath")
    >>> cache_dir.name
    'omnipath'

    Parameters
    ----------
    resource: {"hcop", "ncbi", "omnipath"}, optional
        Resource cache to locate. If `None`, return the bonesistools cache
        root.

    Returns
    -------
    pathlib.Path
        Platform-specific user cache path.
    """

    root = _user_cache_path()
    if resource is None:
        return root

    return root / _as_cache_resource(resource)


def clear(*resources: str) -> None:
    """
    Remove locally cached biological resources.

    Examples
    --------
    >>> import bonesistools as bt
    >>> bt.resources.cache.clear("omnipath")
    >>> bt.resources.cache.clear()

    Parameters
    ----------
    *resources: {"hcop", "ncbi", "omnipath"}
        Resource caches to remove. If omitted, remove every biological
        resource cache managed by bonesistools.
    """

    selected = resources or _CACHE_RESOURCES
    for resource in selected:
        shutil.rmtree(path(_as_cache_resource(resource)), ignore_errors=True)


def __dir__() -> List[str]:

    return sorted(__all__)


def _cached_download(
    url: str,
    *,
    resource: str,
    category: str,
    max_age: Optional[float],
    suffix: Optional[str] = None,
) -> Path:
    """Download a URL atomically or return its locally cached response."""

    resource_name = _as_cache_resource(resource)
    cache_dir = path(resource_name) / category
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{_url_key(url)}{suffix or _url_suffix(url)}"
    metadata_file = cache_file.with_name(f"{cache_file.name}.json")

    if _is_fresh(cache_file, max_age):
        return cache_file

    metadata = _read_metadata(metadata_file, url)
    headers = dict(_DOWNLOAD_HEADERS)
    if cache_file.exists():
        if "etag" in metadata:
            headers["If-None-Match"] = metadata["etag"]
        if "last_modified" in metadata:
            headers["If-Modified-Since"] = metadata["last_modified"]

    request = Request(url, headers=headers)
    try:
        with urlopen(request) as raw_response:
            response = cast(_DownloadResponse, raw_response)
            _write_response(cache_file, response)
            _write_metadata(
                metadata_file,
                {
                    "url": url,
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                },
            )
    except HTTPError as error:
        if error.code == 304 and cache_file.exists():
            cache_file.touch()
            return cache_file
        if cache_file.exists() and error.code >= 500:
            _warn_stale_cache(resource_name, url, error)
            return cache_file
        raise
    except (OSError, URLError) as error:
        if cache_file.exists():
            _warn_stale_cache(resource_name, url, error)
            return cache_file
        raise

    return cache_file


def _user_cache_path() -> Path:

    return _platform_cache_path(
        platform=sys.platform,
        environment=os.environ,
        home=Path.home(),
    )


def _platform_cache_path(
    *,
    platform: str,
    environment: Mapping[str, str],
    home: Path,
) -> Path:

    if platform == "win32":
        local_app_data = environment.get("LOCALAPPDATA")
        root = Path(local_app_data) if local_app_data else home / "AppData" / "Local"
        return root / "bonesistools" / "Cache"

    xdg_cache_home = environment.get("XDG_CACHE_HOME")
    if xdg_cache_home:
        return Path(xdg_cache_home).expanduser() / "bonesistools"

    if platform == "darwin":
        return home / "Library" / "Caches" / "bonesistools"

    return home / ".cache" / "bonesistools"


def _as_cache_resource(resource: str) -> str:

    if not isinstance(resource, str):
        raise TypeError(
            f"unsupported argument type for 'resource': "
            f"expected {str} but received {type(resource)}"
        )

    resource_name = resource.strip().lower()
    if resource_name not in _CACHE_RESOURCES:
        available = ", ".join(repr(name) for name in _CACHE_RESOURCES)
        raise ValueError(
            f"unknown resource cache {resource!r}; available resources are {available}"
        )

    return resource_name


def _url_key(url: str) -> str:

    return hashlib.sha256(url.encode("utf-8")).hexdigest()


def _url_suffix(url: str) -> str:

    suffixes = Path(urlsplit(url).path).suffixes
    return "".join(suffixes) or ".cache"


def _is_fresh(cache_file: Path, max_age: Optional[float]) -> bool:

    if not cache_file.exists():
        return False
    if max_age is None:
        return True

    return time.time() - cache_file.stat().st_mtime <= max_age


def _read_metadata(metadata_file: Path, url: str) -> Dict[str, str]:

    try:
        loaded_metadata: object = json.loads(metadata_file.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return {}

    if not isinstance(loaded_metadata, dict):
        return {}

    metadata = cast(Dict[object, object], loaded_metadata)
    if metadata.get("url") != url:
        return {}

    result: Dict[str, str] = {}
    for key in ("etag", "last_modified"):
        value = metadata.get(key)
        if isinstance(value, str) and value:
            result[key] = value
    return result


def _write_response(cache_file: Path, response: _DownloadResponse) -> None:

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{cache_file.name}.",
        suffix=".tmp",
        dir=str(cache_file.parent),
    )
    temporary_file = Path(temporary_name)
    expected_size = _content_length(response)
    downloaded = 0
    try:
        with os.fdopen(descriptor, "wb") as handle:
            while True:
                chunk = response.read(_DOWNLOAD_CHUNK_SIZE)
                if not chunk:
                    break
                handle.write(chunk)
                downloaded += len(chunk)

        if expected_size is not None and downloaded != expected_size:
            raise OSError(
                f"incomplete download: expected {expected_size} bytes but "
                f"received {downloaded}"
            )

        temporary_file.replace(cache_file)
    finally:
        _unlink(temporary_file)


def _discard_cached_download(cache_file: Path) -> None:
    """Remove a downloaded response after its derived cache is committed."""

    metadata_file = cache_file.with_name(f"{cache_file.name}.json")
    for file in (cache_file, metadata_file):
        try:
            _unlink(file)
        except OSError:
            pass


def _content_length(response: _DownloadResponse) -> Optional[int]:

    value = response.headers.get("Content-Length")
    return int(value) if value is not None else None


def _write_metadata(metadata_file: Path, metadata: Dict[str, Optional[str]]) -> None:

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{metadata_file.name}.",
        suffix=".tmp",
        dir=str(metadata_file.parent),
    )
    temporary_file = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(
                {key: value for key, value in metadata.items() if value is not None},
                handle,
                sort_keys=True,
            )
        temporary_file.replace(metadata_file)
    finally:
        _unlink(temporary_file)


def _unlink(file: Path) -> None:

    try:
        file.unlink()
    except FileNotFoundError:
        pass


def _warn_stale_cache(resource: str, url: str, error: BaseException) -> None:

    warnings.warn(
        f"could not refresh cached {resource} response from {url!r}; "
        f"using the previous local response ({error})",
        UserWarning,
        stacklevel=3,
    )
