#!/usr/bin/env python

"""
Analysis and utility tools for AnnData and MuData objects.

The `tl` namespace contains matrix extraction helpers, graph construction,
distance calculations, marker statistics, lightweight export functions and
feature classification utilities.
"""

from typing import List as _List

from ._binarize import DEABinarizer
from ._classification import mitochondrial_genes as mitochondrial_genes
from ._classification import ribosomal_genes as ribosomal_genes
from ._clustering import kmeans, leiden, louvain
from ._conversion import anndata_to_dataframe as anndata_to_dataframe
from ._conversion import to_dataframe
from ._embedding import pca, spectral, tsne, umap
from ._graph import paga
from ._markers import calculate_logfoldchanges as calculate_logfoldchanges
from ._markers import (
    dea,
    logfoldchanges,
    ora,
    smirnov_tests,
)
from ._maths import barycenters, pairwise_distances
from ._neighbors import (
    KNNSC,
    knn_graph,
    neighbors,
    shared_neighbors,
)
from ._neighbors import Knnbs as Knnbs
from ._neighbors import kneighbors_graph as kneighbors_graph
from ._regress import regress_out
from ._stats import welch_tests, wilcoxon_tests
from ._utils import get_expression, get_pairwise, get_representation
from ._write import to_csv as to_csv
from ._write import to_mtx as to_mtx
from ._write import to_npz as to_npz

__all__ = [
    "get_expression",
    "get_representation",
    "get_pairwise",
    "DEABinarizer",
    "KNNSC",
    "knn_graph",
    "neighbors",
    "shared_neighbors",
    "kmeans",
    "leiden",
    "louvain",
    "pca",
    "spectral",
    "tsne",
    "umap",
    "to_dataframe",
    "paga",
    "pairwise_distances",
    "barycenters",
    "regress_out",
    "welch_tests",
    "wilcoxon_tests",
    "logfoldchanges",
    "dea",
    "ora",
    "smirnov_tests",
]


def __dir__() -> _List[str]:
    hidden = {
        "Knnbs",
        "anndata_to_dataframe",
        "calculate_logfoldchanges",
        "kneighbors_graph",
        "mitochondrial_genes",
        "ribosomal_genes",
        "to_csv",
        "to_mtx",
        "to_npz",
    }
    return sorted((set(globals()) | set(__all__)) - hidden)
