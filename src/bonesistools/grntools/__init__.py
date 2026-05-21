#!/usr/bin/env python

"""
Deprecated GRN namespace.

This module re-exports interaction-graph utilities for backward
compatibility. Prefer `bonesistools.bpy.ig` for new code.
"""

from ..boolpy.ig import __all__ as __all__
from ..boolpy.ig import *  # type: ignore  # noqa: F401,F403

__all__ = list(__all__)


def __dir__():
    return sorted(set(globals()) | set(__all__))
