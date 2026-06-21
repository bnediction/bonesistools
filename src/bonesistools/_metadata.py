#!/usr/bin/env python

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, Optional, cast


def package_version(distribution: str = "bonesistools") -> str:
    """
    Return the installed package version.

    When the package is imported directly from a source checkout, distribution
    metadata can be unavailable. In that case, read the version from
    `pyproject.toml`.
    """

    source_version = _pyproject_version()
    if source_version is not None:
        return source_version

    installed_version = _installed_version(distribution)
    return installed_version if installed_version is not None else "0+unknown"


def _installed_version(distribution: str) -> Optional[str]:

    try:
        metadata = importlib.import_module("importlib.metadata")
    except ImportError:  # pragma: no cover - Python 3.7 compatibility
        metadata = importlib.import_module("importlib_metadata")

    try:
        return cast(str, cast(Any, metadata).version(distribution))
    except Exception as error:
        if error.__class__.__name__ == "PackageNotFoundError":
            return None
        raise


def _pyproject_version() -> Optional[str]:

    pyproject = _find_pyproject()
    if pyproject is None:
        return None

    in_project = False
    for line in pyproject.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped == "[project]":
            in_project = True
            continue
        if in_project and stripped.startswith("["):
            return None
        if in_project and stripped.startswith("version"):
            _, value = stripped.split("=", 1)
            return value.strip().strip("\"'")

    return None


def _find_pyproject() -> Optional[Path]:

    for parent in Path(__file__).resolve().parents:
        pyproject = parent / "pyproject.toml"
        if pyproject.exists():
            return pyproject

    return None
