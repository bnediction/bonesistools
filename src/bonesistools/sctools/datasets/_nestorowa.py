#!/usr/bin/env python

"""
Utilities for loading the Nestorowa hematopoiesis dataset.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, cast

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from ..._validation import _as_boolean
from ._geo import _download

_NESTOROWA_URLS = {
    "read_counts": (
        "GSE81682_HTSeq_counts.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/"
        "GSE81682_HTSeq_counts.txt.gz",
    ),
    "cell_types": (
        "all_cell_types.txt",
        "http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt",
    ),
    "cluster_ids": (
        "cluster_ids.txt",
        "http://blood.stemcells.cam.ac.uk/data/cluster_ids.txt",
    ),
}

_LABEL_GROUPS = {
    "HSC": ["LTHSC_broad", "STHSC_broad", "LTHSC", "STHSC", "ESLAM", "HSC1"],
    "LMPP": ["LMPP_broad", "LMPP"],
    "MPP": [
        "MPP_broad",
        "MPP1_broad",
        "MPP2_broad",
        "MPP3_broad",
        "MPP",
        "MPP1",
        "MPP2",
        "MPP3",
    ],
    "CMP": ["CMP_broad", "CMP"],
    "MEP": ["MEP_broad", "MEP"],
    "GMP": ["GMP_broad", "GMP"],
}

_CLUSTER_LABEL_FALLBACKS = {
    "purple": ["HSC"],
    "deeppink": ["MEP"],
    "gold": ["CMP", "GMP"],
    "darkturquoise": ["CMP", "MPP", "LMPP"],
}


def _load_nestorowa(
    cache_dir: Optional[Path] = None,
    quiet: bool = False,
) -> AnnData:
    """
    Load the Nestorowa et al. hematopoietic stem and progenitor cell dataset.

    Parameters
    ----------
    quiet: bool (default: False)
        Whether to suppress download progress messages.

    Returns
    -------
    AnnData
        Annotated Nestorowa dataset.

        - `adata.X`: raw read counts;
        - `adata.obs["label"]`: broad hematopoietic cell label;
        - `adata.obs["clusters"]`: original color-cluster labels encoded as
          categorical integer labels;
        - `adata.var["ensembl"]`: original Ensembl gene identifier;
        - `adata.var["symbol"]`: converted official gene symbol when available;
        - `adata.uns["nestorowa"]`: source and labeling metadata.

    References
    ----------
    Nestorowa et al. (2016). A single-cell resolution map of mouse
    hematopoietic stem and progenitor cell differentiation. Blood, 128(8),
    e20-e31.
    """

    quiet = _as_boolean(quiet, "quiet")
    if cache_dir is None:
        raise ValueError("missing Nestorowa cache directory")

    paths = _nestorowa_files(cache_dir)
    _download_nestorowa_files(paths, quiet=quiet)

    return _read_full_nestorowa(paths)


def _download_nestorowa_files(paths: Dict[str, Path], quiet: bool) -> None:

    for key, path in paths.items():
        _, url = _NESTOROWA_URLS[key]
        if not path.exists():
            _download(url, path, quiet=quiet)


def _read_full_nestorowa(paths: Dict[str, Path]) -> AnnData:

    read_counts = pd.read_csv(paths["read_counts"], index_col=0, sep="\t").transpose()
    cell_types = pd.read_csv(paths["cell_types"], index_col=0, sep="\t")
    cluster_ids = cast(
        pd.Series,
        pd.read_csv(
            paths["cluster_ids"],
            header=None,
            index_col=0,
            names=["cluster"],
            sep=" ",
            dtype="category",
        )["cluster"],
    )

    read_counts = read_counts.loc[
        cluster_ids.index,
        read_counts.columns.str.startswith("ENS"),
    ]
    cell_types = cell_types.loc[cluster_ids.index, :]
    metadata, label_counts = _resolve_cell_metadata(cell_types, cluster_ids)

    adata = AnnData(X=csr_matrix(read_counts.to_numpy()))
    adata.obs_names = read_counts.index.astype(str).tolist()
    adata.var_names = read_counts.columns.astype(str).tolist()

    adata.obs["label"] = metadata.loc[adata.obs_names, "label"].astype("category")
    adata.obs["clusters"] = _numeric_cluster_categories(cluster_ids).loc[
        adata.obs_names
    ]

    adata.var["ensembl"] = adata.var_names.astype(str)
    _map_ensembl_to_official_names(adata)
    adata.uns["nestorowa"] = {
        "source": {key: url for key, (_, url) in _NESTOROWA_URLS.items()},
        "labels": label_counts,
    }

    return adata


def _resolve_cell_metadata(
    cell_types: pd.DataFrame,
    cluster_ids: pd.Series,
) -> Tuple[pd.DataFrame, Dict[str, int]]:

    cell_labels = pd.DataFrame(index=cell_types.index)
    for label, columns in _LABEL_GROUPS.items():
        cell_labels[label] = cell_types[columns].sum(axis=1)
    cell_labels = (cell_labels != 0).astype(int)

    metadata = pd.DataFrame(
        index=cell_labels.index,
        columns=pd.Index(["label", "cluster"]),
    )
    label_unique: List[str] = []
    label_multi: List[str] = []
    label_missing: List[str] = []

    for cell in cell_labels.index:
        selected = cell_labels.loc[cell, :] > 0
        n_labels = int(selected.sum())
        if n_labels == 1:
            label = str(cell_labels.columns[selected][0])
            label_unique.append(cell)
        elif n_labels > 1:
            label = _seeded_choice(tuple(cell_labels.columns[selected]))
            label_multi.append(cell)
        else:
            cluster = str(cluster_ids[cell])
            label = _seeded_choice(tuple(_CLUSTER_LABEL_FALLBACKS[cluster]))
            label_missing.append(cell)

        metadata.loc[cell, "label"] = label
        metadata.loc[cell, "cluster"] = cluster_ids[cell]

    return metadata, {
        "unique": len(label_unique),
        "multiple": len(label_multi),
        "missing": len(label_missing),
    }


def _seeded_choice(values: Tuple[str, ...]) -> str:

    rng = np.random.RandomState(2020)

    return str(rng.choice(values, 1)[0])


def _numeric_cluster_categories(cluster_ids: pd.Series) -> pd.Series:

    categories = cluster_ids.cat.categories
    mapping = {category: index for index, category in enumerate(categories)}

    return cluster_ids.cat.rename_categories(mapping)


def _map_ensembl_to_official_names(adata: AnnData) -> None:

    from ..preprocessing import (
        convert_gene_identifiers,
        merge_duplicate_vars,
    )

    convert_gene_identifiers(
        adata,
        axis="var",
        input_identifier_type="ensembl_id",
        output_identifier_type="official_name",
        copy=False,
    )
    convert_gene_identifiers(
        adata,
        axis="var",
        copy=False,
    )
    merge_duplicate_vars(adata, copy=False)
    adata.var["symbol"] = adata.var_names.astype(str)


def _nestorowa_files(directory: Path) -> Dict[str, Path]:
    return {key: directory / filename for key, (filename, _) in _NESTOROWA_URLS.items()}
