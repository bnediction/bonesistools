#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import Any, List, Optional, Union

from ..._validation import _as_boolean
from ..._warnings import _rename_deprecated_arguments, _warn_deprecated
from ...logic.influence_graph import InfluenceGraph
from ..ncbi._identifiers import GeneIdentifiers
from ..ncbi._typing import GeneIdentifiersLike, OutputIdentifierType
from ._archive import (
    OmnipathVersion,
    _normalize_organism,
    _normalize_signed_weights,
    _normalize_version,
    list_interactions_versions,
    load_interactions_version,
)
from ._graph_cache import (
    _cache_source_identity,
    _cached_graph,
    _GraphBuild,
    _snapshot_from_dataframe,
)


@_rename_deprecated_arguments(
    genesyn="identifiers",
    gene_identifier_type="identifier_type",
)
def collectri(
    organism: Union[str, int] = "mouse",
    *,
    split_complexes: bool = False,
    remove_pmid: bool = False,
    identifiers: Optional[GeneIdentifiersLike] = None,
    identifier_type: OutputIdentifierType = "symbol",
    version: OmnipathVersion = "latest",
) -> InfluenceGraph:
    """
    Load a CollecTRI signed regulatory network from OmniPath archives.

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Archived OmniPath versions support human, mouse
        and rat. If a dated archive lacks records for the requested organism,
        non-human networks are translated from human interactions with the
        bonesistools HCOP translator using decoupler-compatible expansion.
    split_complexes: bool (default: False)
        Whether to split regulatory complexes into subunits.
    remove_pmid: bool (default: False)
        Whether to remove publication-reference edge attributes.
    identifiers: GeneIdentifiers, optional
        GeneIdentifiers object used to convert graph node identifiers.
    identifier_type: OutputIdentifierType (default: "symbol")
        Output gene identifier type used when `identifiers` is provided.
    version: str or date (default: "latest")
        OmniPath resource version to load. `"latest"` uses the current OmniPath
        interactions endpoint; dates load archived OmniPath interaction dumps.
        Use `collectri.versions()` to inspect available version labels.

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
    Müller-Dott et al. (2023). Expanding the coverage of regulons from
    high-confidence prior knowledge for accurate estimation of transcription
    factor activities. Nucleic Acids Research, 51(20), 10934-10949.
    """

    if not isinstance(organism, (str, int)):
        raise TypeError(
            f"unsupported argument type for 'organism': "
            f"expected {str} or {int} but received {type(organism)}"
        )
    split_complexes = _as_boolean(split_complexes, "split_complexes")
    remove_pmid = _as_boolean(remove_pmid, "remove_pmid")
    if identifiers is not None and not callable(identifiers):
        raise TypeError(
            f"unsupported argument type for 'identifiers': "
            f"expected {GeneIdentifiers} or compatible callable "
            f"but received {type(identifiers)}"
        )

    organism_name = _normalize_organism(organism)
    version_label = _normalize_version(version)
    snapshot = _cached_graph(
        {
            "resource": "collectri",
            "organism": organism_name,
            "version": version_label,
            "hcop": (
                None
                if organism_name == "human"
                else _cache_source_identity("latest")
            ),
        },
        lambda: _build_canonical_collectri(
            organism=organism,
            version=version,
        ),
    )
    grn = snapshot.graph

    if remove_pmid:
        for _source, _target, data in grn.edges(data=True):
            data.pop("PMID", None)
            data.pop("references", None)

    if identifiers is None:
        return grn

    identifiers(
        grn,
        input_type="name",
        output_type=identifier_type,
        copy=False,
    )
    return grn


def _build_canonical_collectri(
    *,
    organism: Union[str, int],
    version: OmnipathVersion,
) -> _GraphBuild:

    downloads: List[Path] = []
    collectri_db = load_interactions_version(
        "collectri",
        version=version,
        organism=organism,
        _downloads=downloads,
    )
    if "weight" in collectri_db:
        collectri_db = collectri_db.rename(columns={"weight": "sign"})
    collectri_db["sign"] = _normalize_signed_weights(collectri_db["sign"].to_numpy())
    return _snapshot_from_dataframe(collectri_db, downloads=downloads)


def _collectri_versions() -> List[str]:
    """
    List available CollecTRI versions without loading the CollecTRI network.

    Returns
    -------
    list of str
        Accepted version labels. `"latest"` denotes the current OmniPath
        endpoint. Date ranges denote archived web-service snapshots; any date
        within one of these ranges can be passed as `version`.
    """

    return list_interactions_versions()


setattr(collectri, "versions", _collectri_versions)


def load_collectri_grn(*args: Any, **kwargs: Any) -> InfluenceGraph:
    """
    Deprecated alias for `collectri`.
    """

    _warn_deprecated(
        "`load_collectri_grn`",
        replacement="`collectri`",
        stacklevel=2,
    )

    return collectri(*args, **kwargs)
