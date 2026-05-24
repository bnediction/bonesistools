#!/usr/bin/env python

"""
Utilities for Boolean networks and logical modelling.

The `bn` sub-package provides data structures, analyses and helper functions
for BooleanNetworkLike objects used in logical modelling frameworks.
"""

from typing import Any, List

from . import _typing as typing

from ._parser import read_bnet, read_bnet_directory
from ._network import (
    BooleanNetwork,
    BooleanNetworkEnsemble,
    ratio_edge_style,
)

__all__ = [
    "read_bnet",
    "read_bnet_directory",
    "ratio_edge_style",
    "BooleanNetwork",
    "BooleanNetworkEnsemble",
    "bn_to_pydot",
    "typing",
]


def bn_to_pydot(bn: Any, **kwargs: Any) -> Any:
    """
    Deprecated. Convert a BooleanNetworkLike object into a pydot graph.

    Use `BooleanNetwork(bn).to_pydot(**kwargs)` instead.
    """

    import warnings

    warnings.warn(
        "`bt.bpy.bn.bn_to_pydot` is deprecated; use "
        "`bt.bpy.bn.BooleanNetwork(bn).to_pydot()` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return BooleanNetwork(bn).to_pydot(**kwargs)


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
