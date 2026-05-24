#!/usr/bin/env python

"""
Deprecated GRN namespace.

This module re-exports interaction-graph utilities for backward
compatibility. Prefer `bonesistools.bpy.ig` for new code.
"""

from typing import List

from ..boolpy.interaction_graph import __all__ as _INTERACTION_GRAPH_ALL
from ..boolpy.interaction_graph import *  # type: ignore  # noqa: F401,F403

__all__ = list(_INTERACTION_GRAPH_ALL)


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
