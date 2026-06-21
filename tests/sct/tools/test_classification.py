#!/usr/bin/env python

from typing import Any, cast

import pytest

import bonesistools as bt
from bonesistools.sctools.tools import _classification


def test_gene_classification_marks_variable_names(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_gene_synonyms", fake_gene_synonyms_cls)
    mini_adata.var_names = ["mt-Co1", "Rps1", "Other"]

    bt.sct.tl.mitochondrial_genes(mini_adata, key="mt")
    bt.sct.tl.ribosomal_genes(mini_adata, key="rps")

    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]


def test_gene_classification_marks_observation_names(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_gene_synonyms", fake_gene_synonyms_cls)
    obs_adata = mini_adata.copy()
    obs_adata.obs_names = ["mt-Co1", "Rps1", "Other1", "Other2"]

    bt.sct.tl.mitochondrial_genes(obs_adata, axis="obs", key="mt_obs")
    bt.sct.tl.ribosomal_genes(obs_adata, axis="obs", key="rps_obs")

    assert obs_adata.obs["mt_obs"].tolist() == [True, False, False, False]
    assert obs_adata.obs["rps_obs"].tolist() == [False, True, False, False]


def test_gene_classification_warns_on_deprecated_integer_axis(
    monkeypatch,
    mini_adata,
    fake_gene_synonyms_cls,
):
    monkeypatch.setattr(_classification, "create_gene_synonyms", fake_gene_synonyms_cls)
    obs_adata = mini_adata.copy()
    obs_adata.obs_names = ["mt-Co1", "Rps1", "Other1", "Other2"]

    with pytest.warns(FutureWarning):
        bt.sct.tl.ribosomal_genes(obs_adata, axis=0, key="rps_obs")

    assert obs_adata.obs["rps_obs"].tolist() == [False, True, False, False]


def test_gene_classification_rejects_invalid_axis(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.mitochondrial_genes(mini_adata, axis=cast(Any, "bad"))

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.ribosomal_genes(mini_adata, axis=cast(Any, "bad"))
