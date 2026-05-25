#!/usr/bin/env python

from __future__ import annotations

from typing import (
    Any,
    Optional,
    Union,
    cast,
)

import networkx as nx
import pandas as pd

from ..ncbi import GeneSynonyms, OutputIdentifierType


def load_collectri_grn(
    organism: Union[str, int] = "mouse",
    split_complexes: bool = False,
    remove_pmid: bool = False,
    genesyn: Optional[GeneSynonyms] = None,
    gene_identifier_type: OutputIdentifierType = "official_name",
    **kwargs: Any,
) -> nx.MultiDiGraph[Any]:
    """
    Load a signed regulatory network derived from the CollecTRI database [1].

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Accepted values depend on the selected decoupler
        wrapper.
    split_complexes: bool (default: False)
        Whether to split regulatory complexes into subunits.
    remove_pmid: bool (default: False)
        Whether to remove the PMID edge attribute returned by decoupler.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert graph node identifiers.
    gene_identifier_type: OutputIdentifierType (default: "official_name")
        Output gene identifier type used when `genesyn` is provided.
    **kwargs: Any
        Keyword arguments passed to the selected decoupler CollecTRI wrapper.

    Returns
    -------
    nx.MultiDiGraph
        Signed regulatory network. Edges contain the attributes returned by
        decoupler, with `weight` renamed to `sign` and converted to -1 or 1.

    Raises
    ------
    TypeError
        If `organism`, `split_complexes`, `remove_pmid` or `genesyn` has an
        unsupported type.

    References
    ----------
    [1] Müller-Dott et al. (2023). Expanding the coverage of regulons
    from high-confidence prior knowledge for accurate estimation of
    transcription factor activities.
    Nucleic Acids Research, 51(20), 10934-10949 (https://doi.org/10.1093/nar/gkad841)
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

    import decoupler as _dc  # type: ignore

    dc = cast(Any, _dc)

    try:
        collectri_db = dc.get_collectri(
            organism=organism,
            split_complexes=split_complexes,
            **kwargs,
        )
    except AttributeError:
        collectri_db = dc.op.collectri(
            organism=organism,
            remove_complexes=split_complexes,
            **kwargs,
        )
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

    grn = nx.from_pandas_edgelist(
        df=collectri_db,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph,
    )
    if genesyn is None:
        return grn
    elif isinstance(genesyn, GeneSynonyms):
        genesyn(
            grn,
            input_identifier_type="name",
            output_identifier_type=gene_identifier_type,
            copy=False,
        )
        return grn
    else:
        raise TypeError(
            f"unsupported argument type for 'genesyn': "
            f"expected {GeneSynonyms} but received {type(genesyn)}"
        )
