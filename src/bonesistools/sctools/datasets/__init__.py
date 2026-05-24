#!/usr/bin/env python

"""
Bundled example datasets for single-cell workflows.

The `datasets` namespace exposes lightweight AnnData objects used in examples,
tests and tutorials.
"""

from __future__ import annotations

import sys
from typing import List

if sys.version_info >= (3, 9):
    from importlib.resources import files
else:
    from importlib_resources import files

import anndata as ad
from anndata import AnnData

__all__ = [
    "nestorowa",
]


def nestorowa() -> AnnData:
    """
    Load the Nestorowa et al. hematopoietic stem and progenitor cell dataset.
    The dataset is a lightweight version intended for examples, testing and tutorials.

    Returns
    -------
    AnnData
        Preprocessed AnnData object containing highly variable genes.

    Notes
    -----
    Dataset derived from:

        Nestorowa et al. (2016)
        "A single-cell resolution map of mouse hematopoietic stem and "
        "progenitor cell differentiation"
        Blood 128(8): e20-e31
    """

    path = files(__package__) / "nestorowa_hvg.h5ad"

    return ad.read_h5ad(str(path))


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
