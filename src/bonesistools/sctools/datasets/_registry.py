#!/usr/bin/env python

"""
Dataset registry and cache management.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Callable, Dict, List

import pandas as pd
from anndata import AnnData

from ..._validation import _as_boolean
from ._metadata import _DATASETS
from ._nestorowa import _load_nestorowa
from ._pbmc3k import _load_pbmc3k

_DatasetLoader = Callable[[Path, bool], AnnData]

_DATASET_LOADERS: Dict[str, _DatasetLoader] = {
    "nestorowa": _load_nestorowa,
    "pbmc3k": _load_pbmc3k,
}


def load(
    name: str,
    quiet: bool = False,
) -> AnnData:
    """
    Load a built-in single-cell dataset.

    Use `bt.sct.datasets.info(name)` for dataset metadata, including source,
    licensing, and citation information.

    Parameters
    ----------
    name: str
        Registered dataset name.
    quiet: bool (default: False)
        Whether to suppress download progress messages.

    Examples
    --------
    >>> adata = bt.sct.datasets.load("nestorowa")
    >>> adata = bt.sct.datasets.load("pbmc3k")

    Returns
    -------
    AnnData
        Requested dataset.
    """

    dataset_name = _as_dataset_name(name)
    quiet = _as_boolean(quiet, "quiet")
    loader = _DATASET_LOADERS[dataset_name]

    return loader(
        _dataset_cache_dir(dataset_name),
        quiet,
    )


def info(name: str) -> Dict[str, object]:
    """
    Return metadata for a registered dataset.

    Parameters
    ----------
    name: str
        Dataset name.

    Examples
    --------
    >>> metadata = bt.sct.datasets.info("pbmc3k")

    Returns
    -------
    dict
        Dataset metadata, including source, license, URL and citation.
    """

    dataset_name = _as_dataset_name(name)

    metadata: Dict[str, object] = {"name": dataset_name}
    metadata.update(_DATASETS[dataset_name])

    return metadata


def available() -> pd.DataFrame:
    """
    List registered datasets.

    Examples
    --------
    >>> bt.sct.datasets.available()

    Returns
    -------
    pandas.DataFrame
        Dataset table indexed by dataset name.
    """

    records: List[Dict[str, object]] = []
    for name, entry in _DATASETS.items():
        record: Dict[str, object] = {"name": name}
        record.update(entry)
        records.append(record)

    return pd.DataFrame.from_records(records, index="name")


def clear(*names: str) -> None:
    """
    Remove cached dataset files.

    Parameters
    ----------
    *names: str
        Dataset names to clear. If no names are provided, remove the whole
        dataset cache.

    Examples
    --------
    >>> bt.sct.datasets.clear("pbmc3k")
    >>> bt.sct.datasets.clear()
    """

    if len(names) == 0:
        shutil.rmtree(_dataset_cache_root(), ignore_errors=True)
        return

    for name in names:
        dataset_name = _as_dataset_name(name)
        shutil.rmtree(_dataset_cache_dir(dataset_name), ignore_errors=True)


def _as_dataset_name(name: str) -> str:

    if not isinstance(name, str):
        raise TypeError(
            f"unsupported argument type for 'name': "
            f"expected {str} but received {type(name)}"
    )

    dataset_name = name.strip().lower()
    if dataset_name not in _DATASETS:
        available = ", ".join(repr(name) for name in sorted(_DATASETS))
        raise ValueError(
            f"unknown dataset {name!r}; available datasets are "
            f"{available}"
        )

    return dataset_name


def _dataset_cache_root() -> Path:

    root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))

    return root / "bonesistools" / "datasets"


def _dataset_cache_dir(name: str) -> Path:

    return _dataset_cache_root() / name
