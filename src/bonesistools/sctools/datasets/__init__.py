#!/usr/bin/env python

"""
Bundled example datasets for single-cell workflows.

The `datasets` namespace exposes lightweight AnnData objects used in examples,
tests and tutorials.
"""

from __future__ import annotations

import sys as _sys
from typing import List as _List

if _sys.version_info >= (3, 9):
    from importlib.resources import files as _files
else:
    from importlib_resources import files as _files

import anndata as _ad
from anndata import AnnData as _AnnData

from ._geo import from_geo

del annotations

__all__ = [
    "from_geo",
    "nestorowa",
]


def nestorowa() -> _AnnData:
    """
    Load the Nestorowa et al. hematopoietic stem and progenitor cell dataset.
    The dataset is a lightweight version intended for examples, testing and tutorials.

    Returns
    -------
    AnnData
        Preprocessed AnnData object containing highly variable genes.

    References
    ----------
    Nestorowa et al. (2016). A single-cell resolution map of mouse
    hematopoietic stem and progenitor cell differentiation. Blood, 128(8),
    e20-e31.
    """

    path = _files(__package__) / "nestorowa_hvg.h5ad"

    return _ad.read_h5ad(str(path))


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
