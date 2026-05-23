#!/usr/bin/env python

"""
Utilities for biological interaction graphs and regulatory domains.

The `ig` sub-package provides parsers, analyses and graph-related
utilities for interaction networks.
"""

from ._parser import read_interaction_graph
from ._algorithms import walks_from
from ._influence_graph import InfluenceGraph

__all__ = [
    "read_interaction_graph",
    "walks_from",
    "InfluenceGraph",
]


def __dir__():
    return sorted(set(globals()) | set(__all__))
