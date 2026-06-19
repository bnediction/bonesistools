#!/usr/bin/env python

"""
Utilities for biological influence graphs and regulatory domains.

The `ig` sub-package provides parsers, analyses and graph-related
utilities for signed regulatory networks.
"""

from typing import List as _List

from ._algorithms import walks_from
from ._distances import distance, similarity
from ._influence_graph import AggregatedInfluenceGraph, InfluenceGraph
from ._parser import read_influence_graph
from ._scoring import (
    InteractionScore,
    infer_signed_interactions,
    infer_signed_interactions_from_walks,
    interaction_scores_from_walks,
)

__all__ = [
    "read_influence_graph",
    "walks_from",
    "distance",
    "similarity",
    "AggregatedInfluenceGraph",
    "InfluenceGraph",
    "InteractionScore",
    "infer_signed_interactions",
    "infer_signed_interactions_from_walks",
    "interaction_scores_from_walks",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
