#!/usr/bin/env python

from types import SimpleNamespace

import networkx as nx
import numpy as np
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _classification


def test_kneighbors_graph_can_use_observation_names(mini_adata):
    graph = bt.sct.tl.kneighbors_graph(
        mini_adata,
        n_neighbors=1,
        use_rep="X_pca",
        create_using=nx.DiGraph,
        index_or_name="name",
    )

    assert set(graph.nodes) == set(mini_adata.obs_names)
    assert graph.number_of_edges() == mini_adata.n_obs
    assert all("distance" in data for _, _, data in graph.edges(data=True))


def test_kneighbors_graph_rejects_invalid_node_label_mode(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'index_or_name'"):
        bt.sct.tl.kneighbors_graph(
            mini_adata,
            n_neighbors=1,
            use_rep="X_pca",
            index_or_name="bad",
        )


def test_shared_neighbors_stores_distances_connectivities_and_metadata(mini_adata):
    result = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=0,
        snn_key="snn",
        copy=True,
    )

    assert "snn" not in mini_adata.uns
    assert result.uns["snn"]["distances_key"] == "snn_distances"
    assert result.uns["snn"]["connectivities_key"] == "snn_connectivities"
    assert result.obsp["snn_distances"].shape == (4, 4)
    assert result.obsp["snn_connectivities"].shape == (4, 4)


def test_shared_neighbors_validates_source_graph_and_pruning(mini_adata):
    with pytest.raises(KeyError, match="neighborhood graph not found"):
        bt.sct.tl.shared_neighbors(mini_adata, knn_key="missing")

    with pytest.raises(ValueError, match="invalid argument value for 'prune_snn'"):
        bt.sct.tl.shared_neighbors(mini_adata, prune_snn=-1)


def test_get_paga_graph_builds_directed_cluster_graph(mini_adata):
    mini_adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])

    graph = bt.sct.tl.get_paga_graph(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        edges="paga_edges",
        threshold=0.1,
    )

    assert set(graph.nodes) == {"A", "B"}
    assert list(graph.edges) == [("B", "A")]
    assert np.allclose(graph.nodes["A"]["pos"], np.array([0.1, 0.05, 1.05]))


def test_get_paga_graph_validates_threshold(mini_adata):
    mini_adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])

    with pytest.raises(TypeError, match="unsupported argument type for 'threshold'"):
        bt.sct.tl.get_paga_graph(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            edges="paga_edges",
            threshold=1,
        )

    with pytest.raises(ValueError, match="invalid argument value for 'threshold'"):
        bt.sct.tl.get_paga_graph(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            edges="paga_edges",
            threshold=-0.1,
        )


class _FakeGeneSynonyms:
    gene_aliases_mapping = {
        "name": {
            "mt-Co1": SimpleNamespace(value=b"mt_gene"),
            "Rps1": SimpleNamespace(value=b"rps_gene"),
        }
    }

    def get_gene_id(self, gene, input_identifier_type):
        return {
            "mt-Co1": "mt_gene",
            "Rps1": "rps_gene",
            "Other": "other_gene",
        }[gene]


def test_mitochondrial_and_ribosomal_gene_classification(monkeypatch, mini_adata):
    monkeypatch.setattr(_classification, "GeneSynonyms", _FakeGeneSynonyms)
    mini_adata.var_names = ["mt-Co1", "Rps1", "Other"]

    bt.sct.tl.mitochondrial_genes(mini_adata, key="mt")
    bt.sct.tl.ribosomal_genes(mini_adata, key="rps")

    assert mini_adata.var["mt"].tolist() == [True, False, False]
    assert mini_adata.var["rps"].tolist() == [False, True, False]


def test_gene_classification_rejects_invalid_axis(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.mitochondrial_genes(mini_adata, axis="bad")
