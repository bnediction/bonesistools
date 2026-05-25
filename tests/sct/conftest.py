#!/usr/bin/env python

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


@pytest.fixture
def mini_adata():
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 0.0, 3.0],
                [2.0, 1.0, 0.0],
                [0.0, 3.0, 1.0],
                [4.0, 0.0, 2.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(["A", "A", "B", "B"]),
                "batch": ["b1", "b2", "b1", "b2"],
                "score": [0.1, 0.7, 0.2, 0.8],
            },
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(
            {"kind": ["keep", "drop", "keep"]},
            index=["g1", "g2", "g3"],
        ),
    )
    adata.layers["counts"] = adata.X.copy() + 1.0
    adata.obsm["X_pca"] = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.2, 0.1, 1.1],
            [2.0, 2.0, 0.0],
            [2.2, 2.1, 0.1],
        ]
    )
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
        "params": {
            "n_neighbors": 3,
            "n_pcs": 2,
            "use_rep": "X_pca",
        },
    }
    adata.obsp["distances"] = csr_matrix(
        [
            [0.0, 1.0, 2.0, 0.0],
            [1.0, 0.0, 0.0, 2.0],
            [2.0, 0.0, 0.0, 1.0],
            [0.0, 2.0, 1.0, 0.0],
        ]
    )
    adata.obsp["connectivities"] = csr_matrix(
        [
            [1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 1.0],
            [0.0, 1.0, 1.0, 1.0],
        ]
    )
    return adata


@pytest.fixture
def expected_mini_cluster_barycenters():
    return {
        "A": np.array([0.1, 0.05, 1.05]),
        "B": np.array([2.1, 2.05, 0.05]),
    }


@pytest.fixture
def expected_mini_pca2_distances():
    return np.array(
        [
            [0.0, np.sqrt(0.05), np.sqrt(8.0), np.sqrt(9.25)],
            [np.sqrt(0.05), 0.0, np.sqrt(6.85), np.sqrt(8.0)],
            [np.sqrt(8.0), np.sqrt(6.85), 0.0, np.sqrt(0.05)],
            [np.sqrt(9.25), np.sqrt(8.0), np.sqrt(0.05), 0.0],
        ]
    )


@pytest.fixture
def expected_mini_snn_connectivities():
    return np.array(
        [
            [0.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 2.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [2.0, 0.0, 0.0, 0.0],
        ]
    )


@pytest.fixture
def fake_gene_synonyms_cls():
    class FakeGeneSynonyms:
        gene_aliases_mapping = {
            "name": {
                "mt-Co1": SimpleNamespace(value=b"mt_gene"),
                "Rps1": SimpleNamespace(value=b"rps_gene"),
            }
        }

        _official_names = {
            "Tp53": "Trp53",
            "Myc": "Myc",
            "NF-kappaB": "Nfkb1",
            "unknown": "unknown",
            "mt-Co1": "mt-Co1",
            "Rps1": "Rps1",
            "Other": "Other",
        }
        _gene_ids = {
            "mt-Co1": "mt_gene",
            "Rps1": "rps_gene",
            "Other": "other_gene",
        }

        def __call__(
            self,
            data,
            axis="index",
            output_identifier_type="official_name",
            copy=True,
            **_,
        ):
            data = data.copy() if copy else data

            if output_identifier_type == "official_name":
                convert = self.get_official_name
            elif output_identifier_type == "gene_id":
                convert = self.get_gene_id
            else:

                def convert(gene, **__):
                    return gene

            if axis in [0, "index"]:
                data.index = [convert(gene) for gene in data.index]
            elif axis in [1, "columns"]:
                data.columns = [convert(gene) for gene in data.columns]

            return data if copy else None

        def get_official_name(self, gene, **_):
            return self._official_names.get(gene, gene)

        def get_gene_id(self, gene, input_identifier_type="name"):
            return self._gene_ids.get(gene, gene)

    return FakeGeneSynonyms
