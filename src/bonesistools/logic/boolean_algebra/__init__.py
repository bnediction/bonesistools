#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides structures and operations for Boolean algebra,
including partial Boolean values, Kleene truth values, hypercubes,
configuration sets, ROBDDs, implicants and predecessor inference utilities.
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
from ._hypercube import Hypercube, HypercubeChanges, HypercubeCollection
from ._kleene import (
    KleeneValue,
    diff,
    join,
    meet,
)
from ._representation import rule_to_string
from ._robdd import ROBDD
from ._structure import (
    dnf_implicants,
    equivalence,
    prime_implicants,
)

__all__ = [
    "BooleanPredecessorInference",
    "PartialBoolean",
    "KleeneValue",
    "ConfigurationSet",
    "Hypercube",
    "HypercubeChanges",
    "HypercubeCollection",
    "ROBDD",
    "diff",
    "meet",
    "join",
    "equivalence",
    "dnf_implicants",
    "prime_implicants",
    "rule_to_string",
]


def read_hypercubes(*args: _Any, **kwargs: _Any) -> _Dict[str, Hypercube]:
    """
    Deprecated. Read named hypercubes from a CSV, TSV or JSON file.

    Use `bt.logic.io.read_hypercubes(...)` instead.
    """

    _warn_deprecated(
        "`bt.logic.ba.read_hypercubes`",
        replacement="`bt.logic.io.read_hypercubes`",
        stacklevel=2,
    )

    from ._parser import read_hypercubes as _read_hypercubes

    return _read_hypercubes(*args, **kwargs)


def __dir__() -> _List[str]:
    hidden = {"read_hypercubes"}
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
