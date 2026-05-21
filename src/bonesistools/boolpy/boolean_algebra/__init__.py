#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides helper structures and operations for
Boolean algebra, including partial Boolean values and Boolean
differential calculus.
"""

from ._algebra import BooleanDifferentialCalculus
from ._boolean import PartialBoolean
from ._hypercube import Hypercube, HypercubeCollection
from ._structure import expressions_equivalent, dnf_to_structure
from ._representation import rule_to_string
from ._parser import read_hypercube, read_hypercubes
from ._typing import (
    BooleanRule,
    PartialBooleanLike,
    HypercubeLike,
    is_boolean_expression_available,
    is_boolean_expression_like,
    is_boolean_rule_like,
    is_partial_boolean_like,
    is_hypercube_like,
)

__all__ = [
    "BooleanDifferentialCalculus",
    "PartialBoolean",
    "Hypercube",
    "HypercubeCollection",
    "expressions_equivalent",
    "dnf_to_structure",
    "rule_to_string",
    "read_hypercube",
    "read_hypercubes",
    "BooleanRule",
    "PartialBooleanLike",
    "HypercubeLike",
    "is_boolean_expression_available",
    "is_boolean_expression_like",
    "is_boolean_rule_like",
    "is_partial_boolean_like",
    "is_hypercube_like",
]


def __dir__():
    return sorted(set(globals()) | set(__all__))
