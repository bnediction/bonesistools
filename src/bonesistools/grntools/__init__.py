#!/usr/bin/env python

"""
Deprecated GRN namespace.

This module re-exports interaction-graph utilities for backward
compatibility. Prefer `bonesistools.bpy.ig` for new code.
"""

from typing import List

from ..boolpy.interaction_graph import (
    InfluenceGraph,
    InteractionScore,
    infer_signed_interactions,
    infer_signed_interactions_from_walks,
    interaction_scores_from_walks,
    read_interaction_graph,
    walks_from,
)

__all__ = [
    "InfluenceGraph",
    "InteractionScore",
    "infer_signed_interactions",
    "infer_signed_interactions_from_walks",
    "interaction_scores_from_walks",
    "read_interaction_graph",
    "walks_from",
]


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
