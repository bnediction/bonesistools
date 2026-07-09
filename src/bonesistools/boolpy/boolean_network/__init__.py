#!/usr/bin/env python

"""
Utilities for Boolean networks and logical modelling.

The `bn` sub-package provides data structures, analyses and helper functions
for BooleanNetworkLike objects used in logical modelling frameworks.
"""

from typing import Any as _Any
from typing import List as _List

from ..._warnings import _warn_deprecated
from ._network import (
    BooleanNetwork,
    BooleanNetworkEnsemble,
)

__all__ = [
    "BooleanNetwork",
    "BooleanNetworkEnsemble",
]


def read_bnet(*args: _Any, **kwargs: _Any) -> BooleanNetwork:
    """
    Deprecated. Read a Boolean network from a `.bnet` file.

    Use `bt.bpy.io.read_bnet(...)` instead.
    """

    _warn_deprecated(
        "`bt.bpy.bn.read_bnet`",
        replacement="`bt.bpy.io.read_bnet`",
        stacklevel=2,
    )

    from ._parser import read_bnet as _read_bnet

    return _read_bnet(*args, **kwargs)


def read_bnet_directory(*args: _Any, **kwargs: _Any) -> BooleanNetworkEnsemble:
    """
    Deprecated. Read `.bnet` files from a directory.

    Use `bt.bpy.io.read_bnet_directory(...)` instead.
    """

    _warn_deprecated(
        "`bt.bpy.bn.read_bnet_directory`",
        replacement="`bt.bpy.io.read_bnet_directory`",
        stacklevel=2,
    )

    from ._parser import read_bnet_directory as _read_bnet_directory

    return _read_bnet_directory(*args, **kwargs)


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
    hidden = {
        "bn_to_pydot",
        "read_bnet",
        "read_bnet_directory",
        "typing",
    }
    return sorted((set(globals()) | set(__all__)) - hidden)
