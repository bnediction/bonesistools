#!/usr/bin/env python

from __future__ import annotations

import warnings
from typing import Any, Optional, Union

import networkx as nx
import pandas as pd

from ...boolpy.influence_graph import InfluenceGraph
from ..ncbi import GeneSynonyms, OutputIdentifierType
from ._archive import OmnipathVersion, load_interactions_version


def collectri(
    organism: Union[str, int] = "mouse",
    split_complexes: bool = False,
    remove_pmid: bool = False,
    genesyn: Optional[GeneSynonyms] = None,
    gene_identifier_type: OutputIdentifierType = "official_name",
    version: OmnipathVersion = "latest",
) -> InfluenceGraph:
    """
    Load a CollecTRI signed regulatory network from OmniPath archives.

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Archived OmniPath versions support human, mouse
        and rat. If a dated archive lacks records for the requested organism,
        non-human networks are translated from human interactions with
        decoupler's HCOP-based `op.translate` helper.
    split_complexes: bool (default: False)
        Whether to split regulatory complexes into subunits.
    remove_pmid: bool (default: False)
        Whether to remove publication-reference edge attributes.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert graph node identifiers.
    gene_identifier_type: OutputIdentifierType (default: "official_name")
        Output gene identifier type used when `genesyn` is provided.
    version: str or date (default: "latest")
        OmniPath resource version to load. Accepted values are `"latest"` or a
        date such as `"2024-01-01"`.

    Returns
    -------
    InfluenceGraph
        Signed regulatory network.

    Raises
    ------
    TypeError
        If `organism`, `split_complexes`, `remove_pmid` or `genesyn` has an
        unsupported type.

    References
    ----------
    [1] Müller-Dott et al. (2023). Expanding the coverage of regulons
    from high-confidence prior knowledge for accurate estimation of
    transcription factor activities. Nucleic Acids Research, 51(20),
    10934-10949 (https://doi.org/10.1093/nar/gkad841)
    """

    if not isinstance(organism, (str, int)):
        raise TypeError(
            f"unsupported argument type for 'organism': "
            f"expected {str} or {int} but received {type(organism)}"
        )
    if not isinstance(split_complexes, bool):
        raise TypeError(
            f"unsupported argument type for 'split_complexes': "
            f"expected {bool} but received {type(split_complexes)}"
        )
    if not isinstance(remove_pmid, bool):
        raise TypeError(
            f"unsupported argument type for 'remove_pmid': "
            f"expected {bool} but received {type(remove_pmid)}"
        )

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
        remaining_columns = [
            column for column in collectri_db.columns if column not in reference_columns
        ]
        collectri_db = pd.DataFrame.from_records(
            (
                {
                    column: value
                    for column, value in row.items()
                    if column in remaining_columns
                }
                for row in collectri_db.to_dict("records")
            ),
            columns=remaining_columns,
        )
    collectri_db["sign"] = collectri_db["sign"].apply(lambda x: -1 if x < 0 else 1)

    grn = InfluenceGraph(
        nx.from_pandas_edgelist(
            df=collectri_db,
            source="source",
            target="target",
            edge_attr=True,
            create_using=nx.MultiDiGraph,
        )
    )

    if genesyn is None:
        return grn

    if isinstance(genesyn, GeneSynonyms):
        genesyn(
            grn,
            input_identifier_type="name",
            output_identifier_type=gene_identifier_type,
            copy=False,
        )
        return grn

    raise TypeError(
        f"unsupported argument type for 'genesyn': "
        f"expected {GeneSynonyms} but received {type(genesyn)}"
    )


def load_collectri_grn(*args: Any, **kwargs: Any) -> InfluenceGraph:
    """
    Deprecated alias for `collectri`.
    """

    warnings.warn(
        "`load_collectri_grn` is deprecated and will be removed in a future "
        "version; use `collectri` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return collectri(*args, **kwargs)
