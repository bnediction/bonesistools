#!/usr/bin/env python

"""
Deprecated GRN namespace.

This module re-exports influence-graph utilities for backward
compatibility. Prefer `bonesistools.bpy.ig` for new code.
"""

from typing import List as _List

from ..boolpy.influence_graph import (
    InfluenceGraph,
    InteractionScore,
    infer_signed_interactions,
    infer_signed_interactions_from_walks,
    interaction_scores_from_walks,
    read_influence_graph,
    walks_from,
)

__all__ = [
    "InfluenceGraph",
    "InteractionScore",
    "infer_signed_interactions",
    "infer_signed_interactions_from_walks",
    "interaction_scores_from_walks",
    "read_influence_graph",
    "walks_from",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
