#!/usr/bin/env python

"""
Analysis and utility tools for AnnData and MuData objects.

The `tl` namespace contains matrix extraction helpers, graph construction,
distance calculations, marker statistics, lightweight export functions and
feature classification utilities.
"""

from typing import List

from ._utils import choose_mtx_representation, choose_representation

from ._neighbors import (
    Knnbs,
    kneighbors_graph,
    shared_neighbors,
)
from ._conversion import anndata_to_dataframe
from ._graph import get_paga_graph
from ._maths import pairwise_distances, barycenters
from ._regress import regress_out

from ._write import (
    to_csv,
    to_csv_or_mtx,
    to_csv_or_npz,
    to_mtx,
    to_npz,
)

from ._markers import (
    calculate_logfoldchanges,
    hypergeometric_test,
    smirnov_tests,
)

from ._classification import mitochondrial_genes, ribosomal_genes

__all__ = [
    "choose_mtx_representation",
    "choose_representation",
    "Knnbs",
    "kneighbors_graph",
    "shared_neighbors",
    "anndata_to_dataframe",
    "get_paga_graph",
    "pairwise_distances",
    "barycenters",
    "regress_out",
    "to_csv",
    "to_csv_or_mtx",
    "to_csv_or_npz",
    "to_mtx",
    "to_npz",
    "calculate_logfoldchanges",
    "hypergeometric_test",
    "smirnov_tests",
    "mitochondrial_genes",
    "ribosomal_genes",
]


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
