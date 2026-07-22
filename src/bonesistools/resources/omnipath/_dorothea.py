#!/usr/bin/env python

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, List, Optional, Union, cast

import numpy as np
import pandas as pd

from ..._compat import Literal
from ..._validation import _as_boolean
from ..._warnings import _rename_deprecated_arguments, _warn_deprecated
from ...logic.influence_graph import InfluenceGraph
from ..ncbi._identifiers import GeneIdentifiers
from ..ncbi._typing import GeneIdentifiersLike, OutputIdentifierType
from ._archive import (
    DorotheaFlavor,
    HcopVersion,
    OmnipathVersion,
    _deduplicate_dorothea,
    _filter_dataset_evidences,
    _format_dorothea,
    _normalize_organism,
    _normalize_signed_weights,
    _normalize_version,
    _read_omnipath_query,
    _tax_id,
    _translate_hcop,
    _truthy,
    list_interactions_versions,
    load_interactions_version,
)
from ._graph_cache import (
    _cache_source_identity,
    _cached_graph,
    _GraphBuild,
    _snapshot_from_dataframe,
)

OMNIPATH_INTERACTIONS_URL = "https://omnipathdb.org/interactions/?genesymbols=1&"
DorotheaWrapper = Literal["op", "get"]


@_rename_deprecated_arguments(
    genesyn="identifiers",
    gene_identifier_type="identifier_type",
)
def dorothea(
    organism: Union[str, int] = "mouse",
    *,
    levels: Optional[List[str]] = None,
    identifiers: Optional[GeneIdentifiersLike] = None,
    identifier_type: OutputIdentifierType = "symbol",
    version: OmnipathVersion = "latest",
    hcop_version: HcopVersion = "latest",
    compatibility: bool = False,
    flavor: Optional[DorotheaFlavor] = None,
    wrapper: Optional[DorotheaWrapper] = None,
    reload: Optional[bool] = None,
    **kwargs: Any,
) -> InfluenceGraph:
    """
    Load a DoRothEA signed regulatory network from OmniPath.

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Dated OmniPath archives support human, mouse and
        rat. If a dated archive lacks records for the requested organism,
        non-human networks are translated from human interactions with the
        bonesistools HCOP translator using decoupler-compatible expansion.
    levels: list of str, optional
        DoRothEA confidence levels to keep. If `None`, use `["A", "B", "C"]`.
    identifiers: GeneIdentifiers, optional
        GeneIdentifiers object used to convert graph node identifiers.
    identifier_type: OutputIdentifierType (default: "symbol")
        Output gene identifier type used when `identifiers` is provided.
    version: str or date (default: "latest")
        OmniPath resource version to load. `"latest"` uses the current OmniPath
        interactions endpoint with decoupler-compatible post-processing; dates
        load archived OmniPath interaction dumps. Use `dorothea.versions()` to
        inspect available version labels.
    hcop_version: "latest", "bundled", str path or Path (default: "latest")
        HCOP version used to translate human interactions to non-human
        organisms when required. `"latest"` downloads the current HGNC HCOP
        file. `"bundled"` uses the HCOP snapshots distributed with
        bonesistools. Paths resolve to custom HCOP-like files.
    compatibility: bool (default: False)
        If `True`, reproduce decoupler's source-target deduplication after
        orthology translation. If `False`, sign-conflicting interactions are
        preserved as distinct signed influences.
    flavor: {"modern", "legacy"}, optional
        DoRothEA construction flavor. If `None`, use `"modern"`.
        If `"modern"`, reproduce the current `decoupler.op.dorothea`
        construction.
        If `"legacy"`, reproduce the historical `decoupler.get_dorothea`
        construction.
    wrapper: {"op", "get"}, optional
        Deprecated alias for `flavor`. `wrapper="op"` maps to
        `flavor="modern"` and `wrapper="get"` maps to `flavor="legacy"`.
    reload: bool, optional
        Deprecated and ignored. Kept for compatibility with older calls.
    **kwargs: Any
        Deprecated and ignored. Kept for compatibility with older decoupler
        wrapper arguments.

    Returns
    -------
    InfluenceGraph
        Signed regulatory network.

    Notes
    -----
    Loaded networks remain cached until explicitly removed. Non-human
    translations may also use cached HCOP files. Use
    `bt.resources.cache.clear("omnipath", "hcop")` to remove both caches.

    References
    ----------
    Garcia-Alonso et al. (2019). Benchmark and integration of resources for
    the estimation of human transcription factor activities. Genome Research,
    29(8), 1363.
    """

    if levels is None:
        levels = ["A", "B", "C"]

    if not isinstance(organism, (str, int)):
        raise TypeError(
            f"unsupported argument type for 'organism': "
            f"expected {str} or {int} but received {type(organism)}"
        )

    compatibility = _as_boolean(compatibility, "compatibility")
    if identifiers is not None and not callable(identifiers):
        raise TypeError(
            f"unsupported argument type for 'identifiers': "
            f"expected {GeneIdentifiers} or compatible callable "
            f"but received {type(identifiers)}"
        )

    if wrapper is not None:
        if wrapper == "op":
            wrapper_flavor: DorotheaFlavor = "modern"
        elif wrapper == "get":
            wrapper_flavor = "legacy"
        else:
            raise ValueError(
                "invalid argument value for 'wrapper': expected 'op' or "
                f"'get', but received {wrapper!r}"
            )

        if flavor is not None and flavor != wrapper_flavor:
            raise ValueError(
                "conflicting DoRothEA construction arguments: "
                f"flavor={flavor!r} and wrapper={wrapper!r}"
            )

        warnings.warn(
            "`wrapper` is deprecated and will be removed in 2.0.0; "
            "use `flavor='modern'` instead of `wrapper='op'` and "
            "`flavor='legacy'` instead of `wrapper='get'`.",
            FutureWarning,
            stacklevel=3,
        )
        flavor = wrapper_flavor

    if flavor is None:
        flavor = "modern"
    elif flavor not in ["modern", "legacy"]:
        raise ValueError(
            "invalid argument value for 'flavor': expected 'modern' or "
            f"'legacy', but received {flavor!r}"
        )

    if reload is not None:
        warnings.warn(
            "`reload` is deprecated and will be removed in 2.0.0; DoRothEA "
            "versions are controlled with `version`.",
            FutureWarning,
            stacklevel=3,
        )
    if kwargs:
        warnings.warn(
            "additional DoRothEA wrapper keyword arguments are deprecated and "
            f"will be removed in 2.0.0; ignored: {', '.join(sorted(kwargs))}",
            FutureWarning,
            stacklevel=3,
        )

    organism_name = _normalize_organism(organism)
    version_label = _normalize_version(version)
    hcop_identity = (
        None
        if organism_name == "human"
        or (version_label == "latest" and flavor == "legacy")
        else _cache_source_identity(hcop_version)
    )
    # Level selection precedes orthology translation and edge deduplication.
    snapshot = _cached_graph(
        {
            "resource": "dorothea",
            "organism": organism_name,
            "version": version_label,
            "flavor": flavor,
            "hcop": hcop_identity,
            "levels": levels,
            "compatibility": compatibility,
        },
        lambda: _build_dorothea_graph(
            organism=organism,
            levels=levels,
            version=version,
            version_label=version_label,
            hcop_version=hcop_version,
            flavor=flavor,
            compatibility=compatibility,
        ),
    )
    grn = snapshot.graph

    if identifiers is None:
        return grn

    identifiers(
        grn,
        input_type="name",
        output_type=identifier_type,
        copy=False,
    )
    return grn


def _dorothea_versions() -> List[str]:
    """
    List available DoRothEA versions without loading the DoRothEA network.

    Returns
    -------
    list of str
        Accepted version labels. `"latest"` denotes the current OmniPath
        endpoint. Date ranges denote archived web-service snapshots; any date
        within one of these ranges can be passed as `version`.
    """

    return list_interactions_versions()


setattr(dorothea, "versions", _dorothea_versions)


def _load_latest_dorothea(
    organism: Union[str, int],
    levels: List[str],
    hcop_version: HcopVersion,
    flavor: DorotheaFlavor,
    compatibility: bool = False,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:

    if flavor == "legacy":
        return _load_latest_legacy_dorothea(
            organism=organism,
            levels=levels,
            compatibility=compatibility,
            _downloads=_downloads,
        )

    return _load_latest_modern_dorothea(
        organism=organism,
        levels=levels,
        hcop_version=hcop_version,
        compatibility=compatibility,
        _downloads=_downloads,
    )


def _load_latest_modern_dorothea(
    organism: Union[str, int],
    levels: List[str],
    hcop_version: HcopVersion,
    compatibility: bool = False,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:

    organism_name = _normalize_organism(organism)
    url = (
        OMNIPATH_INTERACTIONS_URL
        + f"datasets=dorothea&dorothea_levels={','.join(levels)}"
        + "&fields=dorothea_level&license=academic"
    )
    dorothea_db = _read_omnipath_query(url, _downloads=_downloads)
    dorothea_db = cast(
        pd.DataFrame,
        dorothea_db[
            [
                "source_genesymbol",
                "target_genesymbol",
                "is_stimulation",
                "is_inhibition",
                "consensus_stimulation",
                "dorothea_level",
            ]
        ],
    )
    dorothea_db = cast(
        pd.DataFrame,
        dorothea_db[
            ~dorothea_db.duplicated(
                subset=[
                    "source_genesymbol",
                    "dorothea_level",
                    "target_genesymbol",
                    "is_stimulation",
                    "is_inhibition",
                    "consensus_stimulation",
                ]
            )
        ],
    ).copy()
    dorothea_level = cast(pd.Series, dorothea_db["dorothea_level"])
    dorothea_db["dorothea_level"] = dorothea_level.astype(str).str.split(";").str[0]

    stimulation = _truthy(dorothea_db["is_stimulation"])
    inhibition = _truthy(dorothea_db["is_inhibition"])
    consensus_stimulation = _truthy(dorothea_db["consensus_stimulation"])

    signs = np.where(
        inhibition & ~(stimulation & consensus_stimulation),
        -1,
        1,
    )

    dorothea_db = pd.DataFrame(
        {
            "source": dorothea_db["source_genesymbol"],
            "target": dorothea_db["target_genesymbol"],
            "sign": signs,
            "confidence": dorothea_db["dorothea_level"],
        }
    ).sort_values("confidence")
    dorothea_db = cast(
        pd.DataFrame,
        dorothea_db[dorothea_db["confidence"].isin(levels)],
    )

    if organism_name != "human":
        dorothea_db = _translate_hcop(
            cast(pd.DataFrame, dorothea_db),
            organism_name,
            hcop_version=hcop_version,
        )

    return _deduplicate_dorothea(
        cast(pd.DataFrame, dorothea_db),
        compatibility=compatibility,
    )


def _load_latest_legacy_dorothea(
    organism: Union[str, int],
    levels: List[str],
    compatibility: bool = False,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:

    organism_name = _normalize_organism(organism)
    url = (
        "https://omnipathdb.org/interactions?"
        "datasets=dorothea"
        "&dorothea_levels=A,B,C,D"
        "&fields=curation_effort,dorothea_level,evidences,extra_attrs,"
        "references,sources"
        "&format=tsv"
        "&genesymbols=1"
        f"&organisms={_tax_id(organism_name)}"
    )
    dorothea_db = _read_omnipath_query(url, _downloads=_downloads)
    dorothea_db = _filter_dataset_evidences(dorothea_db, dataset="dorothea")
    dorothea_db = _format_dorothea(
        dorothea_db,
        source_column="source_genesymbol",
        target_column="target_genesymbol",
        version_label="latest",
        levels=levels,
    )
    return _deduplicate_dorothea(dorothea_db, compatibility=compatibility)


def _build_dorothea_graph(
    *,
    organism: Union[str, int],
    levels: List[str],
    version: OmnipathVersion,
    version_label: str,
    hcop_version: HcopVersion,
    flavor: DorotheaFlavor,
    compatibility: bool,
) -> _GraphBuild:

    downloads: List[Path] = []
    if version_label == "latest":
        dorothea_db = _load_latest_dorothea(
            organism=organism,
            levels=levels,
            hcop_version=hcop_version,
            flavor=flavor,
            compatibility=compatibility,
            _downloads=downloads,
        )
    else:
        dorothea_db = load_interactions_version(
            "dorothea",
            version=version,
            organism=organism,
            levels=levels,
            flavor=flavor,
            hcop_version=hcop_version,
            compatibility=compatibility,
            _downloads=downloads,
        )

    if "weight" in dorothea_db:
        dorothea_db = dorothea_db.rename(columns={"weight": "sign"})
    if "confidence" in dorothea_db:
        dorothea_db = cast(
            pd.DataFrame,
            dorothea_db[dorothea_db["confidence"].isin(levels)],
        )
    dorothea_db["sign"] = _normalize_signed_weights(dorothea_db["sign"].to_numpy())
    dorothea_db = _deduplicate_dorothea(
        dorothea_db,
        compatibility=compatibility,
    )
    return _snapshot_from_dataframe(dorothea_db, downloads=downloads)


def load_dorothea_grn(*args: Any, **kwargs: Any) -> InfluenceGraph:
    """
    Deprecated alias for `dorothea`.
    """

    _warn_deprecated(
        "`load_dorothea_grn`",
        replacement="`dorothea`",
        stacklevel=2,
    )

    return dorothea(*args, **kwargs)
