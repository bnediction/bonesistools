#!/usr/bin/env python

from __future__ import annotations

import json
import re
import sys
from datetime import date
from typing import Any, Iterable, List, Optional, Sequence, Tuple, Union
from urllib.request import urlopen

import pandas as pd

from ..._compat import Literal

OmnipathDateString = str
OmnipathVersion = Union[Literal["latest"], OmnipathDateString, date]
DorotheaFlavor = Literal["modern", "legacy"]

OMNIPATH_ARCHIVE_URL = "https://archive.omnipathdb.org"
OMNIPATH_INTERACTIONS_PREFIX = "omnipath_webservice_interactions__"
OMNIPATH_INTERACTIONS_PATTERN = re.compile(
    r"omnipath_webservice_interactions__([0-9]{8})-([0-9]{8})[.]tsv[.]xz"
)

ORGANISM_NAMES = {
    "human": "human",
    "homo sapiens": "human",
    "9606": "human",
    "mouse": "mouse",
    "mus musculus": "mouse",
    "10090": "mouse",
    "rat": "rat",
    "rattus norvegicus": "rat",
    "10116": "rat",
}


def resolve_interactions_archive(
    version: OmnipathVersion = "latest",
) -> Tuple[str, str]:
    """
    Resolve an OmniPath interactions archive URL from a version label or date.

    Parameters
    ----------
    version: str or date (default: "latest")
        Version to resolve. `"latest"` resolves the latest archive file. Dates
        are matched to the archive interval that served the web service content
        at that time.

    Returns
    -------
    tuple of str
        Archive URL and normalized version label.
    """

    version_label = _normalize_version(version)

    if version_label == "latest":
        return (
            f"{OMNIPATH_ARCHIVE_URL}/"
            f"{OMNIPATH_INTERACTIONS_PREFIX}latest.tsv.gz",
            "latest",
        )

    archives = _list_interactions_archives()
    for start, end, filename in archives:
        if start <= version_label <= end:
            return f"{OMNIPATH_ARCHIVE_URL}/{filename}", version_label

    raise ValueError(
        "no OmniPath interactions archive found for version "
        f"{version_label!r}; available date ranges: "
        f"{_format_archive_ranges(archives)}"
    )


def load_interactions_version(
    resource: str,
    version: OmnipathVersion,
    organism: Union[str, int],
    levels: Optional[Sequence[str]] = None,
    flavor: DorotheaFlavor = "modern",
) -> pd.DataFrame:
    """
    Load a signed regulatory resource from an OmniPath interactions version.

    When the archive contains interactions for the requested organism, those
    records are used directly. Otherwise, supported non-human organisms are
    translated from human records with decoupler's HCOP-based `op.translate`
    helper.
    """

    _validate_dorothea_flavor(flavor)
    target_organism = _normalize_organism(organism)

    url, version_label = resolve_interactions_archive(version)
    interactions = _read_interactions_archive(url)

    if resource not in interactions:
        raise ValueError(f"resource column not found in OmniPath archive: {resource}")

    interactions = interactions[_truthy(interactions[resource])]

    if resource == "dorothea" and flavor == "legacy":
        interactions, needs_translation = _filter_organism_interactions(
            interactions,
            target_organism,
        )
    else:
        interactions, needs_translation = _filter_human_interactions(
            interactions,
            target_organism,
        )

    source_column = (
        "source_genesymbol" if "source_genesymbol" in interactions else "source"
    )
    target_column = (
        "target_genesymbol" if "target_genesymbol" in interactions else "target"
    )

    if resource == "dorothea":
        interactions = _filter_dataset_evidences(interactions, dataset=resource)
        result = _format_dorothea(
            interactions,
            source_column=source_column,
            target_column=target_column,
            version_label=version_label,
            levels=levels,
        )
    else:
        result = _format_signed_interactions(
            interactions,
            source_column=source_column,
            target_column=target_column,
            version_label=version_label,
        )

    if needs_translation:
        result = _translate_hcop(result, target_organism)

    if resource == "dorothea":
        return _deduplicate_dorothea(result)

    return (
        result.drop_duplicates()
        .drop_duplicates(subset=["source", "target", "sign"])
        .reset_index(drop=True)
    )


def _format_signed_interactions(
    interactions: pd.DataFrame,
    source_column: str,
    target_column: str,
    version_label: str,
) -> pd.DataFrame:
    stimulation_column = (
        "consensus_stimulation"
        if "consensus_stimulation" in interactions
        else "is_stimulation"
    )
    inhibition_column = (
        "consensus_inhibition"
        if "consensus_inhibition" in interactions
        else "is_inhibition"
    )

    stimulation = _truthy(interactions[stimulation_column])
    inhibition = _truthy(interactions[inhibition_column])
    interactions = interactions[stimulation ^ inhibition].copy()
    stimulation = stimulation.loc[interactions.index]

    result = pd.DataFrame(
        {
            "source": interactions[source_column],
            "target": interactions[target_column],
            "sign": stimulation.map(lambda value: 1 if value else -1),
            "version": version_label,
        }
    )

    _copy_optional_columns(result, interactions, columns=["sources", "references"])
    return result


def _format_dorothea(
    interactions: pd.DataFrame,
    source_column: str,
    target_column: str,
    version_label: str,
    levels: Optional[Sequence[str]],
) -> pd.DataFrame:
    interactions = interactions.copy()
    interactions["dorothea_level"] = (
        interactions["dorothea_level"].astype(str).str.split(";").str[0]
    )

    if levels is not None:
        interactions = interactions[interactions["dorothea_level"].isin(levels)]

    stimulation = _truthy(interactions["is_stimulation"])
    inhibition = _truthy(interactions["is_inhibition"])
    consensus_stimulation = _truthy(interactions["consensus_stimulation"])

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

    result = pd.DataFrame(
        {
            "source": interactions[source_column],
            "target": interactions[target_column],
            "sign": signs,
            "confidence": interactions["dorothea_level"],
            "version": version_label,
        }
    )

    _copy_optional_columns(result, interactions, columns=["sources", "references"])
    return result


def _deduplicate_dorothea(dorothea: pd.DataFrame) -> pd.DataFrame:
    dorothea = dorothea[dorothea["source"] != dorothea["target"]]
    return (
        dorothea.drop_duplicates()
        .drop_duplicates(subset=["source", "target"])
        .reset_index(drop=True)
    )


def _copy_optional_columns(
    target: pd.DataFrame,
    source: pd.DataFrame,
    columns: Sequence[str],
) -> None:
    for column in columns:
        if column in source:
            target[column] = source[column]


def _filter_organism_interactions(
    interactions: pd.DataFrame,
    target_organism: str,
) -> Tuple[pd.DataFrame, bool]:
    tax_columns = ["ncbi_tax_id_source", "ncbi_tax_id_target"]
    if not all(column in interactions for column in tax_columns):
        return interactions, target_organism != "human"

    target_tax_id = _tax_id(target_organism)
    tax_source = interactions["ncbi_tax_id_source"].astype(str)
    tax_target = interactions["ncbi_tax_id_target"].astype(str)
    target_mask = (tax_source == target_tax_id) & (tax_target == target_tax_id)

    if target_mask.any():
        return interactions[target_mask], False

    human_mask = (tax_source == "9606") & (tax_target == "9606")
    return interactions[human_mask], target_organism != "human"


def _filter_human_interactions(
    interactions: pd.DataFrame,
    target_organism: str,
) -> Tuple[pd.DataFrame, bool]:
    tax_columns = ["ncbi_tax_id_source", "ncbi_tax_id_target"]
    if not all(column in interactions for column in tax_columns):
        return interactions, target_organism != "human"

    tax_source = interactions["ncbi_tax_id_source"].astype(str)
    tax_target = interactions["ncbi_tax_id_target"].astype(str)
    human_mask = (tax_source == "9606") & (tax_target == "9606")
    return interactions[human_mask], target_organism != "human"


def _filter_dataset_evidences(
    interactions: pd.DataFrame,
    dataset: str,
) -> pd.DataFrame:
    if "evidences" not in interactions:
        return interactions

    interactions = interactions.copy()
    filtered_evidences = [
        _filter_evidences(_parse_evidences(evidences), dataset=dataset)
        for evidences in interactions["evidences"]
    ]
    evidence_rows = [_flatten_evidences(evidences) for evidences in filtered_evidences]
    keep = [bool(evidences) for evidences in evidence_rows]

    interactions = interactions.loc[keep].copy()
    filtered_evidences = [
        evidences for evidences, keep_row in zip(filtered_evidences, keep) if keep_row
    ]
    evidence_rows = [
        evidences for evidences, keep_row in zip(evidence_rows, keep) if keep_row
    ]

    positive = [
        _flatten_evidences(evidences.get("positive", []))
        for evidences in filtered_evidences
    ]
    negative = [
        _flatten_evidences(evidences.get("negative", []))
        for evidences in filtered_evidences
    ]
    directed = [
        _flatten_evidences(evidences.get("directed", []))
        for evidences in filtered_evidences
    ]

    positive_effort = [_curation_effort(evidences) for evidences in positive]
    negative_effort = [_curation_effort(evidences) for evidences in negative]

    interactions["is_directed"] = [bool(evidences) for evidences in directed]
    interactions["is_stimulation"] = [bool(evidences) for evidences in positive]
    interactions["is_inhibition"] = [bool(evidences) for evidences in negative]
    interactions["consensus_stimulation"] = [
        pos >= neg for pos, neg in zip(positive_effort, negative_effort)
    ]
    interactions["consensus_inhibition"] = [
        pos <= neg for pos, neg in zip(positive_effort, negative_effort)
    ]
    interactions["sources"] = [
        _resources_from(evidences) for evidences in evidence_rows
    ]
    interactions["references"] = [
        _references_from(evidences) for evidences in evidence_rows
    ]

    return interactions


def _filter_evidences(evidences: Any, dataset: str) -> Any:
    if isinstance(evidences, dict):
        if "dataset" in evidences and "resource" in evidences:
            return evidences if evidences.get("dataset") == dataset else []
        return {
            key: _filter_evidences(value, dataset=dataset)
            for key, value in evidences.items()
        }

    if isinstance(evidences, list):
        return [
            evidence
            for evidence in evidences
            if evidence.get("dataset") == dataset
        ]

    return evidences


def _parse_evidences(evidences: Any) -> Any:
    return json.loads(evidences) if isinstance(evidences, str) else evidences


def _flatten_evidences(evidences: Any) -> List[dict]:
    if isinstance(evidences, dict):
        if "dataset" in evidences and "resource" in evidences:
            return [evidences]
        return [
            evidence
            for values in evidences.values()
            for evidence in _flatten_evidences(values)
        ]

    if isinstance(evidences, list):
        return [
            evidence
            for value in evidences
            for evidence in _flatten_evidences(value)
        ]

    return []


def _curation_effort(evidences: Iterable[dict]) -> int:
    return sum(len(evidence.get("references", [])) + 1 for evidence in evidences)


def _resources_from(evidences: Iterable[dict]) -> str:
    return ";".join(
        sorted(
            {
                f"{evidence['resource']}{'_' if evidence.get('via') else ''}"
                f"{evidence.get('via') or ''}"
                for evidence in evidences
            }
        )
    )


def _references_from(evidences: Iterable[dict]) -> str:
    return ";".join(
        sorted(
            {
                f"{evidence['resource']}:{reference}"
                for evidence in evidences
                for reference in evidence.get("references", [])
            }
        )
    )


def _normalize_version(version: OmnipathVersion) -> str:
    if isinstance(version, date):
        return version.strftime("%Y%m%d")

    if not isinstance(version, str):
        raise TypeError(
            f"unsupported argument type for 'version': "
            f"expected {str} or {date} but received {type(version)}"
        )

    version_label: str = version.strip().lower()

    if version_label in ["", "latest"]:
        return "latest"

    date_match = re.fullmatch(
        r"([0-9]{4})-?([0-9]{2})-?([0-9]{2})",
        version_label,
    )
    if date_match:
        return "".join(date_match.groups())

    raise ValueError(
        "invalid argument value for 'version': expected 'latest' or a date "
        "formatted as YYYY-MM-DD, but received "
        f"{version_label!r}"
    )


def _validate_dorothea_flavor(flavor: str) -> None:
    if flavor not in ["modern", "legacy"]:
        raise ValueError(
            "invalid argument value for 'flavor': expected 'modern' or "
            f"'legacy', but received {flavor!r}"
        )


def _format_archive_ranges(archives: Sequence[Tuple[str, str, str]]) -> str:
    if not archives:
        return "none"

    return ", ".join(
        f"{_format_archive_date(start)}..{_format_archive_date(end)}"
        for start, end, _filename in archives
    )


def _format_archive_date(version_label: str) -> str:
    return (
        f"{version_label[0:4]}-{version_label[4:6]}-{version_label[6:8]}"
        if re.fullmatch(r"[0-9]{8}", version_label)
        else version_label
    )


def _list_interactions_archives() -> List[Tuple[str, str, str]]:
    with urlopen(f"{OMNIPATH_ARCHIVE_URL}/") as response:
        index = response.read().decode("utf-8")

    archives = [
        (match.group(1), match.group(2), match.group(0))
        for match in OMNIPATH_INTERACTIONS_PATTERN.finditer(index)
    ]
    return sorted(set(archives))


def _read_interactions_archive(url: str) -> pd.DataFrame:
    columns = {
        "source",
        "target",
        "source_genesymbol",
        "target_genesymbol",
        "is_stimulation",
        "is_inhibition",
        "consensus_stimulation",
        "consensus_inhibition",
        "dorothea",
        "collectri",
        "dorothea_level",
        "evidences",
        "sources",
        "references",
        "ncbi_tax_id_source",
        "ncbi_tax_id_target",
    }

    return pd.read_csv(
        url,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in columns,
        low_memory=False,
    )


def _normalize_organism(organism: Union[str, int]) -> str:
    organism_key = str(organism).lower().replace("-", " ")
    organism_name = ORGANISM_NAMES.get(organism_key)

    if organism_name is None:
        raise ValueError(
            "invalid argument value for 'organism': archived OmniPath "
            "versions support human, mouse and rat, but received "
            f"{organism!r}"
        )

    return organism_name


def _tax_id(organism: str) -> str:
    return {
        "human": "9606",
        "mouse": "10090",
        "rat": "10116",
    }[organism]


def _translate_hcop(net: pd.DataFrame, target_organism: str) -> pd.DataFrame:
    import decoupler as _dc  # type: ignore

    dc = getattr(sys.modules.get("decoupler"), "op", None)
    dc = dc or getattr(_dc, "op", None)

    if dc is None or not hasattr(dc, "translate"):
        raise AttributeError(
            "OmniPath archive translation requires decoupler.op.translate "
            "with HCOP support"
        )

    translated = dc.translate(
        net,
        columns=["source", "target"],
        target_organism=target_organism,
    )

    return translated


def _truthy(values: Any) -> pd.Series:
    return pd.Series(values).astype(str).str.lower().isin(["true", "1", "yes"])
