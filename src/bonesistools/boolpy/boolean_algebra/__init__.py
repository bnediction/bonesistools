#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides structures and operations for Boolean algebra,
including partial Boolean values, Kleene truth values, hypercubes,
configuration sets, implicants and predecessor inference utilities.
"""

from typing import Any as _Any
from typing import Dict as _Dict
from typing import List as _List

from ..._warnings import _warn_deprecated
from ._algebra import (
    BooleanPredecessorInference,
)
from ._boolean import PartialBoolean
from ._configuration import ConfigurationSet
from ._hypercube import Hypercube, HypercubeCollection
from ._kleene import (
    KleeneValue,
    diff,
    join,
    meet,
)
from ._representation import rule_to_string
from ._structure import (
    dnf_implicants,
    expressions_equivalent,
    prime_implicants,
)

__all__ = [
    "BooleanPredecessorInference",
    "PartialBoolean",
    "KleeneValue",
    "ConfigurationSet",
    "Hypercube",
    "HypercubeCollection",
    "diff",
    "meet",
    "join",
    "expressions_equivalent",
    "dnf_implicants",
    "prime_implicants",
    "rule_to_string",
]


def read_hypercube(*args: _Any, **kwargs: _Any) -> Hypercube:
    """
    Deprecated. Read a hypercube from a JSON file.

    Use `bt.bpy.io.read_hypercube(...)` instead.
    """

    _warn_deprecated(
        "`bt.bpy.ba.read_hypercube`",
        replacement="`bt.bpy.io.read_hypercube`",
        stacklevel=2,
    )

    from ._parser import read_hypercube as _read_hypercube

    return _read_hypercube(*args, **kwargs)


def read_hypercubes(*args: _Any, **kwargs: _Any) -> _Dict[str, Hypercube]:
    """
    Deprecated. Read named hypercubes from a CSV, TSV or JSON file.

    Use `bt.bpy.io.read_hypercubes(...)` instead.
    """

    _warn_deprecated(
        "`bt.bpy.ba.read_hypercubes`",
        replacement="`bt.bpy.io.read_hypercubes`",
        stacklevel=2,
    )

    from ._parser import read_hypercubes as _read_hypercubes

    return _read_hypercubes(*args, **kwargs)


def __dir__() -> _List[str]:
    hidden = {"read_hypercube", "read_hypercubes"}
    return sorted((set(globals()) | set(__all__)) - hidden)
