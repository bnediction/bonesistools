#!/usr/bin/env python

"""
Compatibility imports shared across the package.
"""

try:
    from typing import Literal, get_args
except ImportError:
    from typing_extensions import Literal, get_args

__all__ = [
    "Literal",
    "get_args",
]
