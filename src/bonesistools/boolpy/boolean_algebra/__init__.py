#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides helper structures and operations for
Boolean algebra, including partial Boolean values and Boolean
differential and predecessor inference utilities.
"""

from typing import Any as _Any
from typing import Dict as _Dict
from typing import List as _List

from ..._warnings import _warn_deprecated
from ._algebra import (
    BooleanPredecessorInference,
)
from ._boolean import PartialBoolean
from ._hypercube import Hypercube, HypercubeCollection
from ._kleene import (
    KleeneValue,
    KleeneValueLike,
    diff,
    join,
    meet,
)
from ._representation import rule_to_string
from ._structure import dnf_to_structure, expressions_equivalent
from ._typing import (
    BooleanRule,
    ConfigurationLike,
    HypercubeLike,
    PartialBooleanLike,
    is_boolean_expression_available,
    is_boolean_expression_like,
    is_boolean_rule_like,
    is_configuration_like,
    is_hypercube_like,
    is_kleene_value_like,
    is_partial_boolean_like,
)

__all__ = [
    "BooleanPredecessorInference",
    "PartialBoolean",
    "KleeneValue",
    "Hypercube",
    "HypercubeCollection",
    "diff",
    "meet",
    "join",
    "expressions_equivalent",
    "dnf_to_structure",
    "rule_to_string",
    "read_hypercube",
    "read_hypercubes",
    "BooleanRule",
    "ConfigurationLike",
    "KleeneValueLike",
    "PartialBooleanLike",
    "HypercubeLike",
    "is_boolean_expression_available",
    "is_boolean_expression_like",
    "is_boolean_rule_like",
    "is_configuration_like",
    "is_kleene_value_like",
    "is_partial_boolean_like",
    "is_hypercube_like",
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
    return sorted(set(globals()) | set(__all__))
