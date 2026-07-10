#!/usr/bin/env python

"""
Dataset registry and cache management.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional

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
_OBSERVATION_COUNT_FIELDS = ("cells", "barcodes", "observations")
_FEATURE_COUNT_FIELDS = ("genes", "peaks")


def load(
    name: str,
    quiet: bool = False,
) -> AnnData:
    """
    Load a built-in single-cell dataset.

    Use `bt.omics.io.info(name)` for dataset metadata, including source,
    licensing, and citation information.

    Parameters
    ----------
    name: str
        Registered dataset name.
    quiet: bool (default: False)
        Whether to suppress download progress messages.

    Examples
    --------
    >>> adata = bt.omics.io.load("nestorowa")
    >>> adata = bt.omics.io.load("pbmc3k")

    Returns
    -------
    AnnData
        Requested dataset.
    """

    dataset_name = _as_dataset_name(name)
    quiet = _as_boolean(quiet, "quiet")
    loader = _DATASET_LOADERS[dataset_name]
    adata = loader(
        _dataset_cache_dir(dataset_name),
        quiet,
    )
    _validate_dataset_shape(dataset_name, adata)

    return adata


def info(name: str) -> Dict[str, object]:
    """
    Return metadata for a registered dataset.

    Parameters
    ----------
    name: str
        Dataset name.

    Examples
    --------
    >>> metadata = bt.omics.io.info("pbmc3k")

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
    >>> bt.omics.io.available()

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
    >>> bt.omics.io.clear("pbmc3k")
    >>> bt.omics.io.clear()
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
            f"unknown dataset {name!r}; available datasets are " f"{available}"
        )

    return dataset_name


def _validate_dataset_shape(
    name: str,
    adata: AnnData,
) -> None:

    metadata = _DATASETS[name]
    expected_observations = _expected_count(
        metadata,
        _OBSERVATION_COUNT_FIELDS,
    )
    expected_features = _expected_feature_count(metadata)
    errors = []
    if expected_observations is not None and adata.n_obs != expected_observations:
        errors.append(
            f"expected {expected_observations} observations but found " f"{adata.n_obs}"
        )
    if expected_features is not None and adata.n_vars != expected_features:
        errors.append(f"expected {expected_features} features but found {adata.n_vars}")
    if len(errors) == 0:
        return

    details = "; ".join(errors)
    raise ValueError(
        f"invalid cached dataset {name!r}: {details}. "
        f"Clear cached files with bt.omics.io.clear({name!r}) and retry."
    )


def _expected_count(
    metadata: Dict[str, Any],
    fields: Iterable[str],
) -> Optional[int]:

    for field in fields:
        if field in metadata:
            return _metadata_count(metadata[field], field)

    return None


def _expected_feature_count(metadata: Dict[str, Any]) -> Optional[int]:

    if "features" in metadata:
        return _metadata_count(metadata["features"], "features")

    counts = [
        _metadata_count(metadata[field], field)
        for field in _FEATURE_COUNT_FIELDS
        if field in metadata
    ]
    if len(counts) == 0:
        return None

    return sum(counts)


def _metadata_count(
    value: object,
    field: str,
) -> int:

    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(
            f"invalid dataset metadata for '{field}': expected {int} "
            f"but received {type(value)}"
        )
    if value < 0:
        raise ValueError(
            f"invalid dataset metadata for '{field}': expected non-negative "
            f"value but received {value!r}"
        )

    return value


def _dataset_cache_root() -> Path:

    root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))

    return root / "bonesistools" / "datasets"


def _dataset_cache_dir(name: str) -> Path:

    return _dataset_cache_root() / name
