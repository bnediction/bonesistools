#!/usr/bin/env python

from pathlib import Path

import bonesistools as bt


def test_version_matches_pyproject():
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
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
