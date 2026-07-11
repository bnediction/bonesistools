#!/usr/bin/env python

"""
Input/output helpers for Boolean modelling objects.

The `io` sub-package centralizes readers and writers for Boolean networks,
influence graphs and Boolean algebra objects.
"""

from typing import List as _List

from ._boolean_network import read_bnet, read_bnet_directory
from ._ginml import read_ginml, read_zginml
from ._hypercube import read_hypercube, read_hypercubes
from ._influence_graph import read_influence_graph
from ._sbml import read_sbml

__all__ = [
    "read_bnet",
    "read_bnet_directory",
    "read_ginml",
    "read_hypercube",
    "read_hypercubes",
    "read_influence_graph",
    "read_sbml",
    "read_zginml",
]


def __dir__() -> _List[str]:
    return sorted(name for name in set(globals()) | set(__all__) if name[0] != "_")
