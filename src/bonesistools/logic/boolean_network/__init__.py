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
from ._symbolic import SymbolicConfigurationSet, SymbolicTransitionSystem

__all__ = [
    "BooleanNetwork",
    "BooleanNetworkEnsemble",
    "SymbolicTransitionSystem",
    "SymbolicConfigurationSet",
]


def read_bnet(*args: _Any, **kwargs: _Any) -> BooleanNetwork:
    """
    Deprecated. Read a Boolean network from a `.bnet` file.

    Use `bt.logic.io.read_bnet(...)` instead.
    """

    _warn_deprecated(
        "`bt.logic.bn.read_bnet`",
        replacement="`bt.logic.io.read_bnet`",
        stacklevel=2,
    )

    from ._parser import read_bnet as _read_bnet

    return _read_bnet(*args, **kwargs)


def read_bnet_directory(*args: _Any, **kwargs: _Any) -> BooleanNetworkEnsemble:
    """
    Deprecated. Read `.bnet` files from a directory.

    Use `bt.logic.io.read_bnet_directory(...)` instead.
    """

    _warn_deprecated(
        "`bt.logic.bn.read_bnet_directory`",
        replacement="`bt.logic.io.read_bnet_directory`",
        stacklevel=2,
    )

    from ._parser import read_bnet_directory as _read_bnet_directory

    return _read_bnet_directory(*args, **kwargs)


def __dir__() -> _List[str]:
    hidden = {
        "read_bnet",
        "read_bnet_directory",
        "typing",
    }
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
