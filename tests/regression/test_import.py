#!/usr/bin/env python

import importlib
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


def test_root_namespace_exposes_short_aliases_only():
    assert {"omics", "logic", "resources"} <= set(dir(bt))

    for name in ["sct", "bpy", "dbs", "sctools", "boolpy", "databases"]:
        assert name not in dir(bt)


def test_deprecated_root_aliases_remain_accessible():
    aliases = {
        "sct": bt.omics,
        "bpy": bt.logic,
        "dbs": bt.resources,
    }

    for name, module in aliases.items():
        with pytest.warns(FutureWarning, match=f"`bt.{name}`"):
            assert getattr(bt, name) is module


def test_legacy_package_import_paths_remain_compatible():
    aliases = {
        "bonesistools.sctools": bt.omics,
        "bonesistools.boolpy": bt.logic,
        "bonesistools.databases": bt.resources,
    }

    for name, module in aliases.items():
        assert importlib.import_module(name) is module


def test_deprecated_public_names_are_hidden_from_tab_completion():
    datasets_module = importlib.import_module("bonesistools.omics.datasets")
    deprecated_names = {
        bt.logic.ba: ["read_hypercube", "read_hypercubes"],
        bt.logic.bn: ["bn_to_pydot", "read_bnet", "read_bnet_directory"],
        bt.logic.ig: ["read_influence_graph"],
        bt.omics.pp: [
            "regress_out",
            "sort_anndata",
            "standardize_gene_identifiers",
            "transfer_obs_its",
            "transfer_obs_sti",
            "var_names_merge_duplicates",
        ],
        bt.omics.tl: [
            "Knnbs",
            "anndata_to_dataframe",
            "calculate_logfoldchanges",
            "kneighbors_graph",
            "mitochondrial_genes",
            "ribosomal_genes",
            "to_csv",
            "to_mtx",
            "to_npz",
        ],
        bt.omics.pl: [
            "add_graph",
            "boxplot",
            "draw_paga",
            "ecdf_plot",
            "embedding_plot",
            "kde_plot",
        ],
        bt.omics.io: ["nestorowa"],
        datasets_module: [
            "available",
            "clear",
            "from_geo",
            "info",
            "load",
            "nestorowa",
        ],
        bt.resources.ncbi: ["GeneSynonyms"],
        bt.resources.omnipath: ["load_collectri_grn", "load_dorothea_grn"],
    }

    for module, names in deprecated_names.items():
        module_dir = dir(module)
        for name in names:
            assert name not in module_dir
            assert hasattr(module, name)


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
