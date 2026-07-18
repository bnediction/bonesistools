#!/usr/bin/env python

from __future__ import annotations

from typing import Any, List, Optional, Union

from ..._validation import _as_boolean
from ..._warnings import _warn_deprecated
from ...logic.influence_graph import InfluenceGraph
from ..ncbi._genesyn import GeneSynonyms
from ..ncbi._typing import GeneSynonymsLike, OutputIdentifierType
from ._archive import (
    OmnipathVersion,
    _normalize_signed_weights,
    list_interactions_versions,
    load_interactions_version,
)


def collectri(
    organism: Union[str, int] = "mouse",
    split_complexes: bool = False,
    remove_pmid: bool = False,
    genesyn: Optional[GeneSynonymsLike] = None,
    gene_identifier_type: OutputIdentifierType = "symbol",
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
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert graph node identifiers.
    gene_identifier_type: OutputIdentifierType (default: "symbol")
        Output gene identifier type used when `genesyn` is provided.
    version: str or date (default: "latest")
        OmniPath resource version to load. `"latest"` uses the current OmniPath
        interactions endpoint; dates load archived OmniPath interaction dumps.
        Use `collectri.versions()` to inspect available version labels.

    Returns
    -------
    InfluenceGraph
        Signed regulatory network.

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

    collectri_db = load_interactions_version(
        "collectri",
        version=version,
        organism=organism,
    )

    if "weight" in collectri_db:
        collectri_db = collectri_db.rename(columns={"weight": "sign"})
    if remove_pmid:
        reference_columns = [
            column
            for column in ("PMID", "references")
            if column in collectri_db.columns
        ]
        collectri_db = collectri_db.drop(columns=reference_columns)
    collectri_db["sign"] = _normalize_signed_weights(
        collectri_db["sign"].to_numpy()
    )

    grn = InfluenceGraph.from_dataframe(collectri_db)

    if genesyn is None:
        return grn

    if not callable(genesyn):
        raise TypeError(
            f"unsupported argument type for 'genesyn': "
            f"expected {GeneSynonyms} or compatible callable "
            f"but received {type(genesyn)}"
        )

    genesyn(
        grn,
        input_type="name",
        output_type=gene_identifier_type,
        copy=False,
    )
    return grn


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
