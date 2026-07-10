#!/usr/bin/env python

"""
Utilities for biological influence graphs and regulatory domains.

The `ig` sub-package provides parsers, analyses and graph-related
utilities for signed regulatory networks.
"""

from typing import Any as _Any
from typing import List as _List

from ..._warnings import _warn_deprecated
from ._algorithms import walks_from
from ._distances import distance, similarity
from ._influence_graph import AggregatedInfluenceGraph, InfluenceGraph
from ._scoring import (
    infer_signed_interactions,
    infer_signed_interactions_from_walks,
    interaction_scores_from_walks,
)

__all__ = [
    "walks_from",
    "distance",
    "similarity",
    "AggregatedInfluenceGraph",
    "InfluenceGraph",
    "infer_signed_interactions",
    "infer_signed_interactions_from_walks",
    "interaction_scores_from_walks",
]


def read_influence_graph(*args: _Any, **kwargs: _Any) -> object:
    """
    Deprecated. Read an influence graph from a tabular file.

    Use `bt.logic.io.read_influence_graph(...)` instead.
    """

    _warn_deprecated(
        "`bt.logic.ig.read_influence_graph`",
        replacement="`bt.logic.io.read_influence_graph`",
        stacklevel=2,
    )

    from ._parser import read_influence_graph as _read_influence_graph

    return _read_influence_graph(*args, **kwargs)


def __dir__() -> _List[str]:
    hidden = {"read_influence_graph"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
