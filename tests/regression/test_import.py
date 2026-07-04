#!/usr/bin/env python

from pathlib import Path

import pytest

import bonesistools as bt
from bonesistools import _metadata


def test_version_matches_pyproject():
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    expected = None
    in_project = False
    for line in pyproject.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped == "[project]":
            in_project = True
            continue
        if in_project and stripped.startswith("["):
            break
        if in_project and stripped.startswith("version"):
            _, value = stripped.split("=", 1)
            expected = value.strip().strip("\"'")
            break

    assert bt.__version__ == expected


def test_package_version_uses_installed_metadata_when_pyproject_is_unavailable(
    monkeypatch,
):
    monkeypatch.setattr(_metadata, "_pyproject_version", lambda: None)

    assert _metadata.package_version("pip") != "0+unknown"
    assert _metadata.package_version("definitely_missing_bonesistools_dist") == (
        "0+unknown"
    )


def test_installed_version_reraises_unexpected_metadata_errors(monkeypatch):
    class BrokenMetadata:
        @staticmethod
        def version(distribution):
            raise RuntimeError(f"broken metadata for {distribution}")

    monkeypatch.setattr(
        _metadata.importlib,
        "import_module",
        lambda name: BrokenMetadata,
    )

    with pytest.raises(RuntimeError, match="broken metadata"):
        _metadata._installed_version("bonesistools")


def test_pyproject_version_handles_missing_project_version(monkeypatch, tmp_path):
    pyproject = tmp_path / "pyproject.toml"

    monkeypatch.setattr(_metadata, "_find_pyproject", lambda: None)
    assert _metadata._pyproject_version() is None

    pyproject.write_text("[build-system]\nrequires = []\n", encoding="utf-8")
    monkeypatch.setattr(_metadata, "_find_pyproject", lambda: pyproject)
    assert _metadata._pyproject_version() is None

    pyproject.write_text("[project]\nname = 'pkg'\n[tool.ruff]\n", encoding="utf-8")
    assert _metadata._pyproject_version() is None


def test_find_pyproject_returns_none_outside_source_tree(monkeypatch, tmp_path):
    module_path = tmp_path / "pkg" / "_metadata.py"
    module_path.parent.mkdir()
    module_path.write_text("", encoding="utf-8")
    monkeypatch.setattr(_metadata, "__file__", str(module_path))

    assert _metadata._find_pyproject() is None
