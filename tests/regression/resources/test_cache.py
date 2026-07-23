#!/usr/bin/env python

import io
import os
import sys
from email.message import Message
from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

import bonesistools as bt
from bonesistools.resources import cache


class _Response(io.BytesIO):
    def __init__(self, content, headers=None):
        super().__init__(content)
        self.headers = headers or {}


@pytest.fixture
def isolated_cache(monkeypatch, tmp_path):
    root = tmp_path / "bonesistools"
    monkeypatch.setattr(cache, "_user_cache_path", lambda: root)
    return root


def test_cache_path_uses_current_platform_directory():
    if sys.platform == "win32":
        local_app_data = os.environ.get("LOCALAPPDATA")
        base = (
            Path(local_app_data)
            if local_app_data
            else Path.home() / "AppData" / "Local"
        )
        expected = base / "bonesistools" / "Cache"
    else:
        xdg_cache_home = os.environ.get("XDG_CACHE_HOME")
        if xdg_cache_home:
            expected = Path(xdg_cache_home).expanduser() / "bonesistools"
        elif sys.platform == "darwin":
            expected = Path.home() / "Library" / "Caches" / "bonesistools"
        else:
            expected = Path.home() / ".cache" / "bonesistools"

    assert bt.resources.cache.path() == expected


@pytest.mark.parametrize(
    ("platform", "environment", "home", "expected"),
    [
        (
            "linux",
            {},
            Path("/home/user"),
            Path("/home/user/.cache/bonesistools"),
        ),
        (
            "freebsd14",
            {},
            Path("/home/user"),
            Path("/home/user/.cache/bonesistools"),
        ),
        (
            "linux",
            {"XDG_CACHE_HOME": "/cache"},
            Path("/home/user"),
            Path("/cache/bonesistools"),
        ),
        (
            "darwin",
            {},
            Path("/Users/user"),
            Path("/Users/user/Library/Caches/bonesistools"),
        ),
        (
            "darwin",
            {"XDG_CACHE_HOME": "/cache"},
            Path("/Users/user"),
            Path("/cache/bonesistools"),
        ),
        (
            "win32",
            {"LOCALAPPDATA": "/local-app-data"},
            Path("/Users/user"),
            Path("/local-app-data/bonesistools/Cache"),
        ),
        (
            "win32",
            {},
            Path("/Users/user"),
            Path("/Users/user/AppData/Local/bonesistools/Cache"),
        ),
    ],
)
def test_platform_cache_path_matches_supported_platform_conventions(
    platform,
    environment,
    home,
    expected,
):
    assert (
        cache._platform_cache_path(
            platform=platform,
            environment=environment,
            home=home,
        )
        == expected
    )


def test_cached_response_is_replaced_atomically_on_current_platform(
    monkeypatch,
    tmp_path,
):
    if sys.platform == "win32":
        monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    else:
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

    responses = iter((b"first", b"second"))

    def urlopen(request):
        content = next(responses)
        return _Response(content, {"Content-Length": str(len(content))})

    monkeypatch.setattr(cache, "urlopen", urlopen)
    url = "https://example.test/interactions"

    cached = cache._cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=0,
        suffix=".tsv",
    )
    os.utime(str(cached), (0, 0))
    refreshed = cache._cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=0,
        suffix=".tsv",
    )

    assert refreshed == cached
    assert refreshed.read_bytes() == b"second"
    assert list(refreshed.parent.glob("*.tmp")) == []

    bt.resources.cache.clear("omnipath")
    assert not bt.resources.cache.path("omnipath").exists()


def test_cache_path_is_platform_specific_and_does_not_create_it(isolated_cache):
    assert bt.resources.cache.path() == isolated_cache
    assert bt.resources.cache.path("hcop") == isolated_cache / "hcop"
    assert bt.resources.cache.path("ncbi") == isolated_cache / "ncbi"
    assert bt.resources.cache.path("omnipath") == isolated_cache / "omnipath"
    assert not isolated_cache.exists()


def test_cache_public_namespace_is_restricted_to_cache_operations():
    assert dir(bt.resources.cache) == ["clear", "path"]


def test_cache_clear_removes_selected_or_all_resource_caches(isolated_cache):
    hcop = isolated_cache / "hcop"
    ncbi = isolated_cache / "ncbi"
    omnipath = isolated_cache / "omnipath"
    unrelated = isolated_cache / "datasets"
    hcop.mkdir(parents=True)
    ncbi.mkdir()
    omnipath.mkdir(parents=True)
    unrelated.mkdir()

    bt.resources.cache.clear("omnipath")

    assert not omnipath.exists()
    assert hcop.exists()
    assert ncbi.exists()
    assert unrelated.exists()

    omnipath.mkdir()
    bt.resources.cache.clear()

    assert not hcop.exists()
    assert not ncbi.exists()
    assert not omnipath.exists()
    assert unrelated.exists()


def test_cache_rejects_unknown_resources(isolated_cache):
    with pytest.raises(ValueError, match="unknown resource cache"):
        bt.resources.cache.path("unknown")

    with pytest.raises(TypeError, match="unsupported argument type"):
        bt.resources.cache.clear(1)  # type: ignore[arg-type]


def test_immutable_download_is_reused_without_network_access(
    isolated_cache,
    monkeypatch,
):
    calls = []

    def urlopen(request):
        calls.append(request.full_url)
        return _Response(b"source\ttarget\nA\tB\n", {"Content-Length": "18"})

    monkeypatch.setattr(cache, "urlopen", urlopen)
    url = "https://example.test/archive.tsv"

    first = cache._cached_download(
        url,
        resource="omnipath",
        category="archives",
        max_age=None,
    )
    second = cache._cached_download(
        url,
        resource="omnipath",
        category="archives",
        max_age=None,
    )

    assert first == second
    assert first.read_bytes() == b"source\ttarget\nA\tB\n"
    assert calls == [url]


def test_stale_download_is_revalidated_with_http_metadata(
    isolated_cache,
    monkeypatch,
):
    requests = []

    def first_urlopen(request):
        requests.append(request)
        return _Response(
            b"cached response",
            {
                "Content-Length": "15",
                "ETag": '"resource-version"',
                "Last-Modified": "Tue, 21 Jul 2026 10:00:00 GMT",
            },
        )

    monkeypatch.setattr(cache, "urlopen", first_urlopen)
    url = "https://example.test/interactions"
    cached = cache._cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=0,
        suffix=".tsv",
    )
    os.utime(str(cached), (0, 0))

    def revalidate(request):
        requests.append(request)
        raise HTTPError(request.full_url, 304, "Not Modified", Message(), None)

    monkeypatch.setattr(cache, "urlopen", revalidate)
    result = cache._cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=0,
        suffix=".tsv",
    )

    assert result == cached
    assert result.read_bytes() == b"cached response"
    assert requests[1].get_header("If-none-match") == '"resource-version"'
    assert (
        requests[1].get_header("If-modified-since") == "Tue, 21 Jul 2026 10:00:00 GMT"
    )
    assert result.stat().st_mtime > 0


def test_stale_download_is_available_when_refresh_fails(
    isolated_cache,
    monkeypatch,
):
    monkeypatch.setattr(
        cache,
        "urlopen",
        lambda request: _Response(b"cached", {"Content-Length": "6"}),
    )
    url = "https://example.test/interactions"
    cached = cache._cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=0,
        suffix=".tsv",
    )
    os.utime(str(cached), (0, 0))

    def unavailable(request):
        raise URLError("offline")

    monkeypatch.setattr(cache, "urlopen", unavailable)

    with pytest.warns(UserWarning, match="using the previous local response"):
        result = cache._cached_download(
            url,
            resource="omnipath",
            category="queries",
            max_age=0,
            suffix=".tsv",
        )

    assert result.read_bytes() == b"cached"


def test_incomplete_download_never_replaces_the_cache(
    isolated_cache,
    monkeypatch,
):
    monkeypatch.setattr(
        cache,
        "urlopen",
        lambda request: _Response(b"short", {"Content-Length": "10"}),
    )

    with pytest.raises(OSError, match="incomplete download"):
        cache._cached_download(
            "https://example.test/archive.tsv",
            resource="omnipath",
            category="archives",
            max_age=None,
        )

    archive_dir = isolated_cache / "omnipath" / "archives"
    assert list(archive_dir.glob("*.tsv")) == []
    assert list(archive_dir.glob("*.tmp")) == []
