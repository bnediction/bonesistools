#!/usr/bin/env python

"""
Bundled example datasets for single-cell workflows.

The `datasets` namespace exposes lightweight AnnData objects used in examples,
tests and tutorials.
"""

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

import anndata as ad


def nestorowa():
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

    return ad.read_h5ad(path)


__all__ = [
    "nestorowa",
]


def __dir__():
    return sorted(set(globals()) | set(__all__))
