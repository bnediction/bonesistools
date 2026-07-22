#!/usr/bin/env python

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple, Union, cast
from urllib.parse import urlsplit

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype

from ..._compat import Literal
from ..cache import _cached_download
from ..hcop import orthologs as hcop_orthologs

OmnipathDateString = str
OmnipathVersion = Union[Literal["latest"], OmnipathDateString, date]
DorotheaFlavor = Literal["modern", "legacy"]
HcopVersion = Union[Literal["bundled", "latest"], str, Path]
_EvidenceSummary = Tuple[bool, bool, bool, bool, bool, str, str]

OMNIPATH_ARCHIVE_URL = "https://archive.omnipathdb.org"
OMNIPATH_INTERACTIONS_PREFIX = "omnipath_webservice_interactions__"
OMNIPATH_INTERACTIONS_PATTERN = re.compile(
    r"omnipath_webservice_interactions__([0-9]{8})-([0-9]{8})[.]tsv[.]xz"
)
_LATEST_CACHE_MAX_AGE = 72 * 60 * 60

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
            f"{OMNIPATH_ARCHIVE_URL}/{OMNIPATH_INTERACTIONS_PREFIX}latest.tsv.gz",
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


def list_interactions_versions() -> List[str]:
    """
    List available OmniPath interaction versions without loading interactions.

    Returns
    -------
    list of str
        Accepted version labels. `"latest"` denotes the current OmniPath
        endpoint. Date ranges denote archived web-service snapshots; any date
        within one of these ranges can be passed as `version`.
    """

    return ["latest"] + [
        _format_archive_interval(start, end)
        for start, end, _filename in _list_interactions_archives()
    ]


def load_interactions_version(
    resource: str,
    version: OmnipathVersion,
    organism: Union[str, int],
    levels: Optional[Sequence[str]] = None,
    flavor: DorotheaFlavor = "modern",
    hcop_version: HcopVersion = "latest",
    compatibility: bool = False,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:
    """
    Load a signed regulatory resource from an OmniPath interactions version.

    When the archive contains interactions for the requested organism, those
    records are used directly. Otherwise, supported non-human organisms are
    translated from human records with bonesistools' HCOP translator using
    decoupler-compatible `one_to_many=5` expansion.
    """

    _validate_dorothea_flavor(flavor)
    target_organism = _normalize_organism(organism)

    url, version_label = resolve_interactions_archive(version)
    interactions = (
        _read_interactions_archive(url)
        if _downloads is None
        else _read_interactions_archive(url, _downloads=_downloads)
    )

    if resource not in interactions:
        raise ValueError(f"resource column not found in OmniPath archive: {resource}")

    interactions = cast(pd.DataFrame, interactions[_truthy(interactions[resource])])

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
        result = _translate_hcop(result, target_organism, hcop_version=hcop_version)

    if resource == "dorothea":
        return _deduplicate_dorothea(result, compatibility=compatibility)

    return result.drop_duplicates(
        subset=["source", "target", "sign"]
    ).reset_index(drop=True)


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
    interactions = cast(pd.DataFrame, interactions[stimulation ^ inhibition]).copy()
    stimulation = stimulation.loc[interactions.index]

    result = pd.DataFrame(
        {
            "source": interactions[source_column],
            "target": interactions[target_column],
            "sign": np.where(stimulation.to_numpy(), 1, -1),
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
        interactions = cast(
            pd.DataFrame,
            interactions[interactions["dorothea_level"].isin(levels)],
        )

    stimulation = _truthy(interactions["is_stimulation"])
    inhibition = _truthy(interactions["is_inhibition"])
    consensus_stimulation = _truthy(interactions["consensus_stimulation"])

    signs = np.where(
        inhibition & ~(stimulation & consensus_stimulation),
        -1,
        1,
    )

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


def _deduplicate_dorothea(
    dorothea: pd.DataFrame,
    compatibility: bool = False,
) -> pd.DataFrame:

    dorothea = cast(
        pd.DataFrame,
        dorothea[dorothea["source"] != dorothea["target"]],
    )
    duplicate_subset = (
        ["source", "target"]
        if compatibility
        else [
            "source",
            "target",
            "sign",
        ]
    )
    if compatibility:
        sort_columns = ["source", "target"]
        if "confidence" in dorothea:
            sort_columns.append("confidence")
        sort_columns.append("sign")
        dorothea = cast(
            pd.DataFrame,
            dorothea.sort_values(sort_columns, kind="mergesort"),
        )

    return dorothea.drop_duplicates(subset=duplicate_subset).reset_index(drop=True)


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
        return cast(pd.DataFrame, interactions[target_mask]), False

    human_mask = (tax_source == "9606") & (tax_target == "9606")
    return cast(pd.DataFrame, interactions[human_mask]), target_organism != "human"


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
    return cast(pd.DataFrame, interactions[human_mask]), target_organism != "human"


def _filter_dataset_evidences(
    interactions: pd.DataFrame,
    dataset: str,
) -> pd.DataFrame:

    if "evidences" not in interactions:
        return interactions

    keep: List[bool] = []
    is_directed: List[bool] = []
    is_stimulation: List[bool] = []
    is_inhibition: List[bool] = []
    consensus_stimulation: List[bool] = []
    consensus_inhibition: List[bool] = []
    sources: List[str] = []
    references: List[str] = []

    for evidences in interactions["evidences"]:
        summary = _summarize_dataset_evidences(evidences, dataset=dataset)
        keep.append(summary is not None)
        if summary is None:
            continue

        (
            directed_value,
            stimulation_value,
            inhibition_value,
            consensus_stimulation_value,
            consensus_inhibition_value,
            source_value,
            reference_value,
        ) = summary
        is_directed.append(directed_value)
        is_stimulation.append(stimulation_value)
        is_inhibition.append(inhibition_value)
        consensus_stimulation.append(consensus_stimulation_value)
        consensus_inhibition.append(consensus_inhibition_value)
        sources.append(source_value)
        references.append(reference_value)

    interactions = cast(pd.DataFrame, interactions.loc[keep]).copy()
    interactions["is_directed"] = is_directed
    interactions["is_stimulation"] = is_stimulation
    interactions["is_inhibition"] = is_inhibition
    interactions["consensus_stimulation"] = consensus_stimulation
    interactions["consensus_inhibition"] = consensus_inhibition
    interactions["sources"] = sources
    interactions["references"] = references

    return interactions


def _summarize_dataset_evidences(
    evidences: Any,
    dataset: str,
) -> Optional[_EvidenceSummary]:

    parsed_evidences = _parse_evidences(evidences)
    evidence_rows: List[dict] = []
    positive: List[dict] = []
    negative: List[dict] = []
    directed: List[dict] = []

    if isinstance(parsed_evidences, dict) and not _is_evidence(parsed_evidences):
        for group, values in parsed_evidences.items():
            group_evidences = _dataset_evidences(values, dataset=dataset)
            evidence_rows.extend(group_evidences)
            if group == "positive":
                positive.extend(group_evidences)
            elif group == "negative":
                negative.extend(group_evidences)
            elif group == "directed":
                directed.extend(group_evidences)
    else:
        evidence_rows = _dataset_evidences(parsed_evidences, dataset=dataset)

    if not evidence_rows:
        return None

    positive_effort = _curation_effort(positive)
    negative_effort = _curation_effort(negative)
    return (
        bool(directed),
        bool(positive),
        bool(negative),
        positive_effort >= negative_effort,
        positive_effort <= negative_effort,
        _resources_from(evidence_rows),
        _references_from(evidence_rows),
    )


def _parse_evidences(evidences: Any) -> Any:

    return json.loads(evidences) if isinstance(evidences, str) else evidences


def _dataset_evidences(evidences: Any, dataset: str) -> List[dict]:

    if isinstance(evidences, dict):
        if _is_evidence(evidences):
            return [evidences] if evidences.get("dataset") == dataset else []
        return [
            evidence
            for values in evidences.values()
            for evidence in _dataset_evidences(values, dataset=dataset)
        ]

    if isinstance(evidences, list):
        return [
            evidence
            for value in evidences
            for evidence in _dataset_evidences(value, dataset=dataset)
        ]

    return []


def _is_evidence(value: Any) -> bool:

    return isinstance(value, dict) and "dataset" in value and "resource" in value


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

    return ", ".join(_format_archive_interval(start, end) for start, end, _ in archives)


def _format_archive_interval(start: str, end: str) -> str:

    start_label = _format_archive_date(start)
    end_label = _format_archive_date(end)
    if start_label == end_label:
        return start_label
    return f"{start_label}..{end_label}"


def _format_archive_date(version_label: str) -> str:

    return (
        f"{version_label[0:4]}-{version_label[4:6]}-{version_label[6:8]}"
        if re.fullmatch(r"[0-9]{8}", version_label)
        else version_label
    )


def _list_interactions_archives() -> List[Tuple[str, str, str]]:

    index_file = _cached_download(
        f"{OMNIPATH_ARCHIVE_URL}/",
        resource="omnipath",
        category="queries",
        max_age=_LATEST_CACHE_MAX_AGE,
        suffix=".html",
    )
    index = index_file.read_text(encoding="utf-8")

    archives = [
        (match.group(1), match.group(2), match.group(0))
        for match in OMNIPATH_INTERACTIONS_PATTERN.finditer(index)
    ]
    return sorted(set(archives))


def _read_interactions_archive(
    url: str,
    *,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:

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

    filename = Path(urlsplit(url).path).name
    max_age = (
        _LATEST_CACHE_MAX_AGE
        if filename == f"{OMNIPATH_INTERACTIONS_PREFIX}latest.tsv.gz"
        else None
    )
    archive = _cached_download(
        url,
        resource="omnipath",
        category="archives",
        max_age=max_age,
    )
    if _downloads is not None:
        _downloads.append(archive)

    return pd.read_csv(
        archive,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in columns,
        low_memory=False,
    )


def _read_omnipath_query(
    url: str,
    *,
    _downloads: Optional[List[Path]] = None,
) -> pd.DataFrame:

    response = _cached_download(
        url,
        resource="omnipath",
        category="queries",
        max_age=_LATEST_CACHE_MAX_AGE,
        suffix=".tsv",
    )
    if _downloads is not None:
        _downloads.append(response)
    return pd.read_csv(response, sep="\t", low_memory=False)


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


def _translate_hcop(
    net: pd.DataFrame,
    target_organism: str,
    hcop_version: HcopVersion = "latest",
) -> pd.DataFrame:

    translated = hcop_orthologs(
        target_organism=target_organism,
        version=hcop_version,
    ).translate_dataframe(
        net,
        columns=["source", "target"],
        on_unmapped="drop",
        one_to_many=5,
    )
    return cast(pd.DataFrame, translated)


def _truthy(values: Any) -> pd.Series:

    values = values if isinstance(values, pd.Series) else pd.Series(values)
    if is_bool_dtype(values.dtype):
        return values.fillna(False).astype(bool)

    return values.astype(str).str.lower().isin(["true", "1", "yes"])


def _normalize_signed_weights(values: Any) -> np.ndarray:
    """Normalize signed numeric weights while preserving missing values."""

    weights = np.asarray(values)
    missing = pd.isna(weights)
    normalized = np.full(weights.shape, np.nan)
    normalized[~missing] = np.where(weights[~missing] < 0, -1, 1)
    return normalized
