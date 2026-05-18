#!/usr/bin/env python

"""
Utilities for Boolean networks and logical modelling.

The `bn` sub-package provides data structures, analyses and helper functions
for Boolean network-like objects used in logical modelling frameworks.
"""

from . import _typing as typing

from ._network import BooleanNetwork, BooleanNetworkEnsemble
from ._hypercube import Hypercube, HypercubeCollection

__all__ = [
    "BooleanNetwork",
    "BooleanNetworkEnsemble",
    "Hypercube",
    "HypercubeCollection",
    "bn_to_pydot",
    "typing",
]


def bn_to_pydot(bn, **kwargs):
    """
    Deprecated. Convert a Boolean-network-like object into a pydot graph.

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
