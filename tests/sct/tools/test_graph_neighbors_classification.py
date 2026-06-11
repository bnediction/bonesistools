#!/usr/bin/env python

import networkx as nx
import numpy as np
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt
from bonesistools.sctools.tools import _classification


def test_knn_graph_can_use_observation_names(mini_adata):
    graph = bt.sct.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        use_rep="X_pca",
        create_using=nx.DiGraph,
        index_or_name="name",
    )

    assert set(graph.nodes) == set(mini_adata.obs_names)
    assert graph.number_of_edges() == mini_adata.n_obs
    assert {
        (source, target): data["distance"]
        for source, target, data in graph.edges(data=True)
    } == pytest.approx(
        {
            ("c1", "c2"): np.sqrt(0.06),
            ("c2", "c1"): np.sqrt(0.06),
            ("c3", "c4"): np.sqrt(0.06),
            ("c4", "c3"): np.sqrt(0.06),
        }
    )

    index_graph = bt.sct.tl.knn_graph(
        mini_adata,
        n_neighbors=1,
        use_rep="X_pca",
        create_using=nx.DiGraph,
        index_or_name="index",
    )

    assert set(index_graph.nodes) == set(range(mini_adata.n_obs))


def test_knn_graph_rejects_invalid_node_label_mode(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'index_or_name'"):
        bt.sct.tl.knn_graph(
            mini_adata,
            n_neighbors=1,
            use_rep="X_pca",
            index_or_name="bad",
        )


def test_shared_neighbors_stores_distances_connectivities_and_metadata(
    mini_adata,
    expected_mini_pca2_distances,
    expected_mini_snn_connectivities,
):
    result = bt.sct.tl.shared_neighbors(
        mini_adata,
        prune_snn=0,
        snn_key="snn",
        copy=True,
    )

    assert "snn" not in mini_adata.uns
    assert result.uns["snn"]["distances_key"] == "snn_distances"
    assert result.uns["snn"]["connectivities_key"] == "snn_connectivities"

    expected_distances = np.zeros_like(expected_mini_pca2_distances)
    expected_distances[0, 3] = expected_mini_pca2_distances[0, 3]
    expected_distances[1, 2] = expected_mini_pca2_distances[1, 2]
    expected_distances[2, 1] = expected_mini_pca2_distances[2, 1]
    expected_distances[3, 0] = expected_mini_pca2_distances[3, 0]

    assert np.allclose(
        result.obsp["snn_distances"].toarray(),
        expected_distances,
    )
    assert np.allclose(
        result.obsp["snn_connectivities"].toarray(),
        expected_mini_snn_connectivities / 3,
    )


def test_shared_neighbors_validates_source_graph_and_pruning(mini_adata):
    with pytest.raises(KeyError, match="neighborhood graph not found"):
        bt.sct.tl.shared_neighbors(mini_adata, knn_key="missing")

    with pytest.raises(ValueError, match="invalid argument value for 'prune_snn'"):
        bt.sct.tl.shared_neighbors(mini_adata, prune_snn=-1)


def test_extract_paga_graph_builds_directed_cluster_graph(
    mini_adata,
    expected_mini_cluster_barycenters,
):
    mini_adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])

    graph = bt.sct.tl.extract_paga_graph(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        edges="paga_edges",
        threshold=0.1,
    )

    assert set(graph.nodes) == {"A", "B"}
    assert list(graph.edges) == [("B", "A")]
    assert np.allclose(graph.nodes["A"]["pos"], expected_mini_cluster_barycenters["A"])
    assert np.allclose(graph.nodes["B"]["pos"], expected_mini_cluster_barycenters["B"])

    mini_adata.uns["connectivities"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])
    fallback_graph = bt.sct.tl.extract_paga_graph(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        edges="missing",
        threshold=0.1,
    )

    assert list(fallback_graph.edges) == [("B", "A")]


def test_extract_paga_graph_validates_threshold(mini_adata):
    mini_adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.0, 0.0]])

    with pytest.raises(KeyError, match="key 'missing' not found in adata.uns"):
        bt.sct.tl.extract_paga_graph(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            edges="missing",
        )

    with pytest.raises(TypeError, match="unsupported argument type for 'threshold'"):
        bt.sct.tl.extract_paga_graph(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            edges="paga_edges",
            threshold=1,
        )

    with pytest.raises(ValueError, match="invalid argument value for 'threshold'"):
        bt.sct.tl.extract_paga_graph(
            mini_adata,
            obs="cluster",
            use_rep="X_pca",
            edges="paga_edges",
            threshold=-0.1,
        )


def test_mitochondrial_and_ribosomal_gene_classification(
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

    obs_adata = mini_adata.copy()
    obs_adata.obs_names = ["mt-Co1", "Rps1", "Other1", "Other2"]

    bt.sct.tl.mitochondrial_genes(obs_adata, axis="obs", key="mt_obs")
    bt.sct.tl.ribosomal_genes(obs_adata, axis=0, key="rps_obs")

    assert obs_adata.obs["mt_obs"].tolist() == [True, False, False, False]
    assert obs_adata.obs["rps_obs"].tolist() == [False, True, False, False]


def test_gene_classification_rejects_invalid_axis(mini_adata):
    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.mitochondrial_genes(mini_adata, axis="bad")

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        bt.sct.tl.ribosomal_genes(mini_adata, axis="bad")
