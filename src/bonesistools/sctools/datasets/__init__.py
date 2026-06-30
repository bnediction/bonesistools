#!/usr/bin/env python

"""
Bundled example datasets for single-cell workflows.

The `datasets` namespace exposes AnnData-oriented dataset loaders and
lightweight objects used in examples, tests and tutorials.
"""

from __future__ import annotations

from typing import List as _List

from ._geo import from_geo
from ._nestorowa import nestorowa

del annotations

__all__ = [
    "from_geo",
    "nestorowa",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
