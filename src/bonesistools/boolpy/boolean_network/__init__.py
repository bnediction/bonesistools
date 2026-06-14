#!/usr/bin/env python

"""
Utilities for Boolean networks and logical modelling.

The `bn` sub-package provides data structures, analyses and helper functions
for BooleanNetworkLike objects used in logical modelling frameworks.
"""

from typing import Any as _Any
from typing import List as _List

from ..._warnings import _warn_deprecated
from . import _typing as typing
from ._network import (
    BooleanNetwork,
    BooleanNetworkEnsemble,
    ratio_edge_style,
)
from ._parser import read_bnet, read_bnet_directory

__all__ = [
    "BooleanNetwork",
    "BooleanNetworkEnsemble",
    "ratio_edge_style",
    "read_bnet",
    "read_bnet_directory",
    "bn_to_pydot",
    "typing",
]


def bn_to_pydot(bn: _Any, **kwargs: _Any) -> _Any:
    """
    Deprecated. Convert a BooleanNetworkLike object into a pydot graph.

    Use `BooleanNetwork(bn).to_pydot(**kwargs)` instead.
    """

    _warn_deprecated(
        "`bt.bpy.bn.bn_to_pydot`",
        replacement="`bt.bpy.bn.BooleanNetwork(bn).to_pydot()`",
        stacklevel=2,
    )

    return BooleanNetwork(bn).to_pydot(**kwargs)


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
