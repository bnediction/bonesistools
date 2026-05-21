#!/usr/bin/env python

import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from bonesistools.databases.omnipath import _collectri, _dorothea


def _install_fake_decoupler(monkeypatch, **members):
    fake_decoupler = SimpleNamespace(**members)
    monkeypatch.setitem(sys.modules, "decoupler", fake_decoupler)
    return fake_decoupler


def test_load_dorothea_grn_uses_op_wrapper_and_builds_signed_graph(
    monkeypatch, tmp_path
):
    calls = []

    def dorothea(**kwargs):
        calls.append(kwargs)
        return pd.DataFrame(
            {
                "source": ["Tf1", "Tf2"],
                "target": ["Gene1", "Gene2"],
                "weight": [-0.2, 0.7],
                "confidence": ["A", "B"],
            }
        )

    _install_fake_decoupler(monkeypatch, op=SimpleNamespace(dorothea=dorothea))
    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))

    grn = _dorothea.load_dorothea_grn(
        organism="mouse",
        levels=["A", "B"],
        wrapper="op",
        reload=True,
        license="academic",
    )

    assert calls == [
        {
            "organism": "mouse",
            "levels": ["A", "B"],
            "license": "academic",
        }
    ]
    assert sorted(grn.edges(data=True)) == [
        ("Tf1", "Gene1", {"confidence": "A", "sign": -1}),
        ("Tf2", "Gene2", {"confidence": "B", "sign": 1}),
    ]


def test_load_dorothea_grn_reads_cache_without_importing_decoupler(
    monkeypatch, tmp_path
):
    cache_dir = tmp_path / ".cache"
    cache_dir.mkdir()
    cache_file = cache_dir / "dorothea_op_mouse_A.csv"
    pd.DataFrame(
        {
            "source": ["Tf"],
            "target": ["Gene"],
            "weight": [1.0],
            "confidence": ["A"],
        }
    ).to_csv(cache_file, index=False)

    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))
    monkeypatch.delitem(sys.modules, "decoupler", raising=False)

    grn = _dorothea.load_dorothea_grn(
        organism="mouse",
        levels=["A"],
        wrapper="op",
    )

    assert list(grn.edges(data=True)) == [
        ("Tf", "Gene", {"confidence": "A", "sign": 1})
    ]


def test_load_dorothea_grn_rejects_invalid_wrapper():
    with pytest.raises(ValueError, match="invalid argument value for 'wrapper'"):
        _dorothea.load_dorothea_grn(wrapper="legacy")


def test_load_collectri_grn_uses_op_wrapper_when_legacy_wrapper_is_missing(
    monkeypatch,
):
    calls = []

    def get_collectri(**kwargs):
        raise AttributeError("legacy wrapper is not available")

    def collectri(**kwargs):
        calls.append(kwargs)
        return pd.DataFrame(
            {
                "source": ["Tf1", "Tf2"],
                "target": ["Gene1", "Gene2"],
                "weight": [-1.0, 1.0],
                "PMID": ["1", "2"],
            }
        )

    _install_fake_decoupler(
        monkeypatch,
        get_collectri=get_collectri,
        op=SimpleNamespace(collectri=collectri),
    )

    grn = _collectri.load_collectri_grn(
        organism="mouse",
        split_complexes=True,
        remove_pmid=True,
        license="academic",
    )

    assert calls == [
        {
            "organism": "mouse",
            "remove_complexes": True,
            "license": "academic",
        }
    ]
    assert sorted(grn.edges(data=True)) == [
        ("Tf1", "Gene1", {"sign": -1}),
        ("Tf2", "Gene2", {"sign": 1}),
    ]


def test_load_collectri_grn_rejects_invalid_boolean_arguments():
    with pytest.raises(
        TypeError, match="unsupported argument type for 'split_complexes'"
    ):
        _collectri.load_collectri_grn(split_complexes="yes")

    with pytest.raises(TypeError, match="unsupported argument type for 'remove_pmid'"):
        _collectri.load_collectri_grn(remove_pmid="yes")
