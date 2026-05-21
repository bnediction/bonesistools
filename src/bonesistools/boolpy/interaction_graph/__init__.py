#!/usr/bin/env python

"""
Utilities for biological interaction graphs and regulatory domains.

The `ig` sub-package provides parsers, analyses and graph-related
utilities for interaction networks.
"""

from ._parser import read_interaction_graph
from ._algorithms import depth_first_extraction
from ._graphinfo import get_edge_sign, get_path_sign, path_to_string, statistics

__all__ = [
    "read_interaction_graph",
    "depth_first_extraction",
    "get_edge_sign",
    "get_path_sign",
    "path_to_string",
    "statistics",
]


def __dir__():
    return sorted(set(globals()) | set(__all__))
