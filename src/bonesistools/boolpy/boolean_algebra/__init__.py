#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides helper structures and operations for
Boolean algebra, including partial Boolean values and Boolean
differential and predecessor inference utilities.
"""

from typing import List

from ._algebra import BooleanPredecessorInference, PartialBooleanDifferential
from ._boolean import PartialBoolean
from ._hypercube import Hypercube, HypercubeCollection
from ._parser import read_hypercube, read_hypercubes
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
    is_partial_boolean_like,
)

__all__ = [
    "BooleanPredecessorInference",
    "PartialBooleanDifferential",
    "PartialBoolean",
    "Hypercube",
    "HypercubeCollection",
    "expressions_equivalent",
    "dnf_to_structure",
    "rule_to_string",
    "read_hypercube",
    "read_hypercubes",
    "BooleanRule",
    "ConfigurationLike",
    "PartialBooleanLike",
    "HypercubeLike",
    "is_boolean_expression_available",
    "is_boolean_expression_like",
    "is_boolean_rule_like",
    "is_configuration_like",
    "is_partial_boolean_like",
    "is_hypercube_like",
]


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
