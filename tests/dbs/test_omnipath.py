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
                "source": ["Tf1", "Tf2", "Tf3"],
                "target": ["Gene1", "Gene2", "Gene3"],
                "weight": [-0.2, 0.7, 1.0],
                "confidence": ["A", "B", "D"],
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
            "levels": ["A", "B", "C", "D"],
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
    cache_file = cache_dir / "dorothea_op_mouse_complete.csv"
    pd.DataFrame(
        {
            "source": ["TfA", "TfB"],
            "target": ["GeneA", "GeneB"],
            "weight": [1.0, -1.0],
            "confidence": ["A", "B"],
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
        ("TfA", "GeneA", {"confidence": "A", "sign": 1})
    ]


def test_load_dorothea_grn_get_wrapper_warns_and_filters_locally(monkeypatch, tmp_path):
    calls = []

    def get_dorothea(**kwargs):
        calls.append(kwargs)
        return pd.DataFrame(
            {
                "source": ["TfA", "TfB"],
                "target": ["GeneA", "GeneB"],
                "weight": [1.0, -1.0],
                "confidence": ["A", "B"],
            }
        )

    _install_fake_decoupler(monkeypatch, get_dorothea=get_dorothea)
    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))

    with pytest.warns(UserWarning, match="legacy `decoupler.get_dorothea`"):
        grn = _dorothea.load_dorothea_grn(
            organism="mouse",
            levels=["B"],
            wrapper="get",
            reload=True,
        )

    assert calls == [{"organism": "mouse"}]
    assert list(grn.edges(data=True)) == [
        ("TfB", "GeneB", {"confidence": "B", "sign": -1})
    ]


def test_load_dorothea_grn_get_wrapper_explains_missing_legacy_api(
    monkeypatch, tmp_path
):
    _install_fake_decoupler(monkeypatch, __version__="2.1.6", op=SimpleNamespace())
    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))

    with pytest.warns(UserWarning, match="legacy `decoupler.get_dorothea`"):
        with pytest.raises(
            AttributeError,
            match=r"2\.1\.6.*requires decoupler<2\.0\.0",
        ):
            _dorothea.load_dorothea_grn(
                organism="mouse",
                wrapper="get",
                reload=True,
            )


def test_load_dorothea_grn_op_wrapper_explains_missing_current_api(
    monkeypatch, tmp_path
):
    _install_fake_decoupler(monkeypatch, __version__="1.9.2")
    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))

    with pytest.raises(
        AttributeError,
        match=r"1\.9\.2.*requires decoupler>=2\.0\.0",
    ):
        _dorothea.load_dorothea_grn(
            organism="mouse",
            wrapper="op",
            reload=True,
        )


def test_load_dorothea_grn_default_levels_and_genesyn(monkeypatch, tmp_path):
    class FakeGeneSynonyms:
        def __init__(self):
            self.calls = []

        def __call__(self, graph, **kwargs):
            self.calls.append((graph, kwargs))

    calls = []

    def dorothea(**kwargs):
        calls.append(kwargs)
        return pd.DataFrame(
            {
                "source": ["TfA", "TfB", "TfC", "TfD"],
                "target": ["GeneA", "GeneB", "GeneC", "GeneD"],
                "weight": [1.0, -1.0, 1.0, -1.0],
                "confidence": ["A", "B", "C", "D"],
            }
        )

    genesyn = FakeGeneSynonyms()

    _install_fake_decoupler(monkeypatch, op=SimpleNamespace(dorothea=dorothea))
    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))
    monkeypatch.setattr(_dorothea, "GeneSynonyms", FakeGeneSynonyms)

    grn = _dorothea.load_dorothea_grn(
        organism=9606,
        genesyn=genesyn,
        gene_identifier_type="gene_id",
        wrapper="op",
        reload=True,
    )

    assert calls == [{"organism": 9606, "levels": ["A", "B", "C", "D"]}]
    assert sorted(grn.edges(data=True)) == [
        ("TfA", "GeneA", {"confidence": "A", "sign": 1}),
        ("TfB", "GeneB", {"confidence": "B", "sign": -1}),
        ("TfC", "GeneC", {"confidence": "C", "sign": 1}),
    ]
    assert genesyn.calls == [
        (
            grn,
            {
                "input_identifier_type": "name",
                "output_identifier_type": "gene_id",
                "copy": False,
            },
        )
    ]


def test_load_dorothea_grn_rejects_invalid_wrapper():
    with pytest.raises(ValueError, match="invalid argument value for 'wrapper'"):
        _dorothea.load_dorothea_grn(wrapper="legacy")


def test_load_dorothea_grn_rejects_invalid_organism_and_genesyn(
    monkeypatch,
    tmp_path,
):
    with pytest.raises(TypeError, match="unsupported argument type for 'organism'"):
        _dorothea.load_dorothea_grn(organism=object())

    cache_dir = tmp_path / ".cache"
    cache_dir.mkdir()
    pd.DataFrame(
        {
            "source": ["TfA"],
            "target": ["GeneA"],
            "weight": [1.0],
            "confidence": ["A"],
        }
    ).to_csv(cache_dir / "dorothea_op_mouse_complete.csv", index=False)

    monkeypatch.setattr(_dorothea, "__file__", str(tmp_path / "_dorothea.py"))

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        _dorothea.load_dorothea_grn(genesyn=object())


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

    with pytest.raises(TypeError, match="unsupported argument type for 'organism'"):
        _collectri.load_collectri_grn(organism=object())


def test_load_collectri_grn_uses_legacy_wrapper_and_genesyn(monkeypatch):
    class FakeGeneSynonyms:
        def __init__(self):
            self.calls = []

        def __call__(self, graph, **kwargs):
            self.calls.append((graph, kwargs))

    calls = []

    def get_collectri(**kwargs):
        calls.append(kwargs)
        return pd.DataFrame(
            {
                "source": ["Tf1", "Tf2"],
                "target": ["Gene1", "Gene2"],
                "weight": [1.0, -1.0],
                "references": ["ref1", "ref2"],
            }
        )

    genesyn = FakeGeneSynonyms()

    _install_fake_decoupler(monkeypatch, get_collectri=get_collectri)
    monkeypatch.setattr(_collectri, "GeneSynonyms", FakeGeneSynonyms)

    grn = _collectri.load_collectri_grn(
        organism=9606,
        genesyn=genesyn,
        gene_identifier_type="ensembl_id",
        license="academic",
    )

    assert calls == [
        {
            "organism": 9606,
            "split_complexes": False,
            "license": "academic",
        }
    ]
    assert sorted(grn.edges(data=True)) == [
        ("Tf1", "Gene1", {"references": "ref1", "sign": 1}),
        ("Tf2", "Gene2", {"references": "ref2", "sign": -1}),
    ]
    assert genesyn.calls == [
        (
            grn,
            {
                "input_identifier_type": "name",
                "output_identifier_type": "ensembl_id",
                "copy": False,
            },
        )
    ]


def test_load_collectri_grn_rejects_invalid_genesyn(monkeypatch):
    def get_collectri(**kwargs):
        return pd.DataFrame(
            {
                "source": ["Tf1"],
                "target": ["Gene1"],
                "weight": [1.0],
            }
        )

    _install_fake_decoupler(monkeypatch, get_collectri=get_collectri)

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        _collectri.load_collectri_grn(genesyn=object())
