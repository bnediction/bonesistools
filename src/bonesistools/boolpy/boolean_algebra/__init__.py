#!/usr/bin/env python

"""
Utilities for Boolean algebra and partial Boolean abstractions.

The `ba` sub-package provides helper structures and operations for
Boolean algebra, including partial Boolean values and Boolean
differential calculus.
"""

from ._algebra import BooleanDifferentialCalculus
from ._boolean import PartialBoolean
from ._structure import expressions_equivalent, dnf_to_structure
from ._representation import rule_to_string

__all__ = [
    BooleanDifferentialCalculus,
    PartialBoolean,
    expressions_equivalent,
    dnf_to_structure,
    rule_to_string,
]
