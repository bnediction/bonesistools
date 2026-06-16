#!/usr/bin/env python

"""
Analysis and utility tools for AnnData and MuData objects.

The `tl` namespace contains matrix extraction helpers, graph construction,
distance calculations, marker statistics, lightweight export functions and
feature classification utilities.
"""

from typing import List as _List

from ._classification import mitochondrial_genes, ribosomal_genes
from ._clustering import leiden
from ._conversion import anndata_to_dataframe
from ._embedding import pca, spectral, tsne, umap
from ._graph import extract_paga_graph, get_paga_graph
from ._markers import (
    calculate_logfoldchanges,
    logfoldchanges,
    ora,
    smirnov_tests,
)
from ._maths import barycenters, pairwise_distances
from ._neighbors import (
    KNNSC,
    Knnbs,
    kneighbors_graph,
    knn_graph,
    neighbors,
    shared_neighbors,
)
from ._regress import regress_out
from ._stats import wilcoxon_tests
from ._utils import get_expression, get_pairwise, get_representation
from ._write import (
    to_csv,
    to_csv_or_mtx,
    to_csv_or_npz,
    to_mtx,
    to_npz,
)

__all__ = [
    "get_expression",
    "get_representation",
    "get_pairwise",
    "KNNSC",
    "Knnbs",
    "knn_graph",
    "kneighbors_graph",
    "neighbors",
    "shared_neighbors",
    "leiden",
    "pca",
    "spectral",
    "tsne",
    "umap",
    "anndata_to_dataframe",
    "extract_paga_graph",
    "get_paga_graph",
    "pairwise_distances",
    "barycenters",
    "regress_out",
    "wilcoxon_tests",
    "to_csv",
    "to_csv_or_mtx",
    "to_csv_or_npz",
    "to_mtx",
    "to_npz",
    "logfoldchanges",
    "calculate_logfoldchanges",
    "ora",
    "smirnov_tests",
    "mitochondrial_genes",
    "ribosomal_genes",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
