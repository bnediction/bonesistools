#!/usr/bin/env python

from typing import Any, cast

import pytest

import bonesistools as bt
from bonesistools.omics.preprocessing import _classification


def test_gene_classification_marks_variable_names(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_identifiers", fake_gene_synonyms_cls)
    mini_adata.var_names = ["mt-Co1", "Rps1", "Other"]

    bt.omics.pp.mitochondrial_genes(mini_adata, key="mt")
    bt.omics.pp.ribosomal_genes(mini_adata, key="rps")

    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]


def test_gene_classification_uses_requested_organism(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_identifiers", fake_gene_synonyms_cls)
    mini_adata.var_names = ["MT-ND1", "RPS1", "Other"]

    bt.omics.pp.mitochondrial_genes(mini_adata, key="mt", organism="human")
    bt.omics.pp.ribosomal_genes(mini_adata, key="rps", organism="human")

    assert fake_gene_synonyms_cls.organism == "human"
    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]


def test_gene_classification_reuses_identifiers(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    def raise_if_called(**_):
        raise AssertionError("gene identifier converter should be reused")

    monkeypatch.setattr(_classification, "create_identifiers", raise_if_called)
    identifiers = fake_gene_synonyms_cls(organism="human")
    mini_adata.var_names = ["MT-ND1", "RPS1", "Other"]

    bt.omics.pp.mitochondrial_genes(
        mini_adata,
        key="mt",
        organism="mouse",
        identifiers=identifiers,
    )
    bt.omics.pp.ribosomal_genes(
        mini_adata,
        key="rps",
        organism="mouse",
        identifiers=identifiers,
    )

    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]


@pytest.mark.parametrize("function_name", ["mitochondrial_genes", "ribosomal_genes"])
def test_gene_classification_accepts_deprecated_genesyn(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
    function_name,
):
    def raise_if_called(**_):
        raise AssertionError("gene identifier converter should be reused")

    monkeypatch.setattr(_classification, "create_identifiers", raise_if_called)
    identifiers = fake_gene_synonyms_cls(organism="human")
    mini_adata.var_names = ["MT-ND1", "RPS1", "Other"]
    function = getattr(bt.omics.pp, function_name)

    with pytest.warns(FutureWarning, match="genesyn.*identifiers"):
        function(mini_adata, genesyn=identifiers)


def test_ribosomal_gene_classification_uses_official_symbols(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_identifiers", fake_gene_synonyms_cls)
    mini_adata.var_names = ["Rps1", "Mrpl56", "Lactb"]

    bt.omics.pp.ribosomal_genes(mini_adata, key="rps")

    assert mini_adata.var["rps"].tolist() == [True, False, False]


def test_gene_classification_marks_observation_names(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_identifiers", fake_gene_synonyms_cls)
    obs_adata = mini_adata.copy()
    obs_adata.obs_names = ["mt-Co1", "Rps1", "Other1", "Other2"]

    bt.omics.pp.mitochondrial_genes(obs_adata, axis="obs", key="mt_obs")
    bt.omics.pp.ribosomal_genes(obs_adata, axis="obs", key="rps_obs")

    assert obs_adata.obs["mt_obs"].tolist() == [True, False, False, False]
    assert obs_adata.obs["rps_obs"].tolist() == [False, True, False, False]


def test_gene_classification_warns_on_deprecated_integer_axis(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_identifiers", fake_gene_synonyms_cls)
    obs_adata = mini_adata.copy()
    obs_adata.obs_names = ["mt-Co1", "Rps1", "Other1", "Other2"]

    with pytest.warns(FutureWarning):
        bt.omics.pp.ribosomal_genes(obs_adata, axis=0, key="rps_obs")

    assert obs_adata.obs["rps_obs"].tolist() == [False, True, False, False]


def test_gene_classification_rejects_invalid_axis(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.omics.pp.mitochondrial_genes(mini_adata, axis=cast(Any, "bad"))

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.omics.pp.ribosomal_genes(mini_adata, axis=cast(Any, "bad"))


def test_gene_classification_keeps_deprecated_tools_aliases(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    from bonesistools.omics.tools import _classification as tools_classification

    monkeypatch.setattr(
        tools_classification,
        "create_identifiers",
        fake_gene_synonyms_cls,
    )
    mini_adata.var_names = ["mt-Co1", "Rps1", "Other"]

    with pytest.warns(FutureWarning):
        getattr(bt.omics.tl, "mitochondrial_genes")(mini_adata, key="mt")
    with pytest.warns(FutureWarning):
        getattr(bt.omics.tl, "ribosomal_genes")(mini_adata, key="rps")

    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]
