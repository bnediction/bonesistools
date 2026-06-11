#!/usr/bin/env python

from __future__ import annotations

import warnings
from typing import Any, List, Optional, Union, cast

import networkx as nx
import pandas as pd

from ..._compat import Literal
from ...boolpy.influence_graph import InfluenceGraph
from ..ncbi import GeneSynonyms, OutputIdentifierType
from ._archive import (
    DorotheaFlavor,
    HcopVersion,
    OmnipathVersion,
    _deduplicate_dorothea,
    _filter_dataset_evidences,
    _format_dorothea,
    _normalize_organism,
    _tax_id,
    _translate_hcop,
    _truthy,
    list_interactions_versions,
    load_interactions_version,
)

OMNIPATH_INTERACTIONS_URL = "https://omnipathdb.org/interactions/?genesymbols=1&"
DorotheaWrapper = Literal["op", "get"]


def dorothea(
    organism: Union[str, int] = "mouse",
    levels: Optional[List[str]] = None,
    genesyn: Optional[GeneSynonyms] = None,
    gene_identifier_type: OutputIdentifierType = "official_name",
    version: OmnipathVersion = "latest",
    hcop_version: HcopVersion = "latest",
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
        DoRothEA confidence levels to keep. If None, use `["A", "B", "C"]`.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert graph node identifiers.
    gene_identifier_type: OutputIdentifierType (default: "official_name")
        Output gene identifier type used when `genesyn` is provided.
    version: str or date (default: "latest")
        OmniPath resource version to load. `"latest"` uses the current OmniPath
        interactions endpoint with decoupler-compatible post-processing; dates
        load archived OmniPath interaction dumps. Use `dorothea.versions()` to
        inspect available version labels.
    hcop_version: "latest", "bundled" or pathlib.Path (default: "latest")
        HCOP version used to translate human interactions to non-human
        organisms when required. `"latest"` downloads the current HGNC HCOP
        file. `"bundled"` uses the HCOP snapshots distributed with
        bonesistools. Paths resolve to custom HCOP-like files.
    flavor: {"modern", "legacy"}, optional
        DoRothEA construction flavor. If None, use `"modern"`.
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

    Raises
    ------
    TypeError
        If `organism` or `genesyn` has an unsupported type.
    """

    if levels is None:
        levels = ["A", "B", "C"]

    if not isinstance(organism, (str, int)):
        raise TypeError(
            f"unsupported argument type for 'organism': "
            f"expected {str} or {int} but received {type(organism)}"
        )

    if genesyn is not None and not isinstance(genesyn, GeneSynonyms):
        raise TypeError(
            f"unsupported argument type for 'genesyn': "
            f"expected {GeneSynonyms} but received {type(genesyn)}"
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
            "`wrapper` is deprecated and will be removed in a future version; "
            "use `flavor='modern'` instead of `wrapper='op'` and "
            "`flavor='legacy'` instead of `wrapper='get'`.",
            DeprecationWarning,
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
            "`reload` is deprecated and ignored; DoRothEA versions are controlled "
            "with `version`.",
            DeprecationWarning,
            stacklevel=3,
        )
    if kwargs:
        warnings.warn(
            "additional DoRothEA wrapper keyword arguments are deprecated and "
            f"ignored: {', '.join(sorted(kwargs))}",
            DeprecationWarning,
            stacklevel=3,
        )

    if isinstance(version, str) and version.strip().lower() in ["", "latest"]:
        dorothea_db = _load_latest_dorothea(
            organism=organism,
            levels=levels,
            hcop_version=hcop_version,
            flavor=flavor,
        )
    else:
        dorothea_db = load_interactions_version(
            "dorothea",
            version=version,
            organism=organism,
            levels=levels,
            flavor=flavor,
            hcop_version=hcop_version,
        )

    if "confidence" in dorothea_db:
        dorothea_columns = list(dorothea_db.columns)
        dorothea_db = pd.DataFrame.from_records(
            (
                row
                for row in dorothea_db.to_dict("records")
                if row["confidence"] in levels
            ),
            columns=dorothea_columns,
        )
    if "weight" in dorothea_db:
        dorothea_db = dorothea_db.rename(columns={"weight": "sign"})
    dorothea_db["sign"] = dorothea_db["sign"].apply(lambda x: -1 if x < 0 else 1)

    grn = InfluenceGraph(
        nx.from_pandas_edgelist(
            df=dorothea_db,
            source="source",
            target="target",
            edge_attr=True,
            create_using=nx.MultiDiGraph,
        )
    )

    if genesyn is None:
        return grn

    genesyn(
        grn,
        input_identifier_type="name",
        output_identifier_type=gene_identifier_type,
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
) -> pd.DataFrame:
    if flavor == "legacy":
        return _load_latest_legacy_dorothea(organism=organism, levels=levels)

    return _load_latest_modern_dorothea(
        organism=organism,
        levels=levels,
        hcop_version=hcop_version,
    )


def _load_latest_modern_dorothea(
    organism: Union[str, int],
    levels: List[str],
    hcop_version: HcopVersion,
) -> pd.DataFrame:
    organism_name = _normalize_organism(organism)
    url = (
        OMNIPATH_INTERACTIONS_URL
        + f"datasets=dorothea&dorothea_levels={','.join(levels)}"
        + "&fields=dorothea_level&license=academic"
    )
    dorothea_db = pd.read_csv(url, sep="\t", low_memory=False)
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
                subset=["source_genesymbol", "dorothea_level", "target_genesymbol"]
            )
        ],
    ).copy()
    dorothea_level = cast(pd.Series, dorothea_db["dorothea_level"])
    dorothea_db["dorothea_level"] = dorothea_level.astype(str).str.split(";").str[0]

    stimulation = _truthy(dorothea_db["is_stimulation"])
    inhibition = _truthy(dorothea_db["is_inhibition"])
    consensus_stimulation = _truthy(dorothea_db["consensus_stimulation"])

    signs = []
    for is_stimulation, is_inhibition, is_consensus_stimulation in zip(
        stimulation,
        inhibition,
        consensus_stimulation,
    ):
        if is_stimulation and is_inhibition:
            signs.append(1 if is_consensus_stimulation else -1)
        elif is_stimulation:
            signs.append(1)
        elif is_inhibition:
            signs.append(-1)
        else:
            signs.append(1)

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

    return cast(
        pd.DataFrame,
        dorothea_db.drop_duplicates(subset=["source", "target"]).reset_index(drop=True),
    )


def _load_latest_legacy_dorothea(
    organism: Union[str, int],
    levels: List[str],
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
    dorothea_db = pd.read_csv(url, sep="\t", low_memory=False)
    dorothea_db = _filter_dataset_evidences(dorothea_db, dataset="dorothea")
    dorothea_db = _format_dorothea(
        dorothea_db,
        source_column="source_genesymbol",
        target_column="target_genesymbol",
        version_label="latest",
        levels=levels,
    )
    return _deduplicate_dorothea(dorothea_db)


def load_dorothea_grn(*args: Any, **kwargs: Any) -> InfluenceGraph:
    """
    Deprecated alias for `dorothea`.
    """

    warnings.warn(
        "`load_dorothea_grn` is deprecated and will be removed in a future "
        "version; use `dorothea` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return dorothea(*args, **kwargs)
