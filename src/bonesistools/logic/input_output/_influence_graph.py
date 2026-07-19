#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import (
    Any,
    Optional,
    Union,
)

import networkx as nx
import pandas as pd

from ..._warnings import _rename_deprecated_arguments
from ...resources.ncbi._typing import (
    GeneIdentifiersLike,
    InputIdentifierType,
    OutputIdentifierType,
)


@_rename_deprecated_arguments(
    genesyn="identifiers",
    gene_type="input_type",
    alias_gene="output_type",
    input_identifier_type="input_type",
    output_identifier_type="output_type",
)
def read_influence_graph(
    file: Union[str, Path],
    *,
    identifiers: Optional[GeneIdentifiersLike] = None,
    input_type: InputIdentifierType = "name",
    output_type: OutputIdentifierType = "symbol",
    sep: str = ",",
    **kwargs: Any,
) -> nx.MultiDiGraph[Any]:
    """
    Read an influence graph from a tabular file.

    The input file must contain at least three columns: 'source', 'target' and
    'sign'. Additional columns are kept as edge attributes.

    Parameters
    ----------
    file: str | Path
        Path to the tabular file containing the influence graph.
    identifiers: GeneIdentifiers, optional
        GeneIdentifiers object used to convert graph node identifiers.
    input_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
        (default: 'name')
        Input gene identifier type. Valid database-specific values are listed
        in `databases`.
    output_type: 'symbol' | 'ncbi_symbol' | 'gene_id' |
        'ensembl_id' | <database> (default: 'symbol')
        Output gene identifier type. Valid database-specific values are listed
        in `databases`.
    sep: str (default: ",")
        Field delimiter passed to pandas.read_csv.
    **kwargs: Any
        Additional keyword arguments passed to pandas.read_csv.

    Returns
    -------
    nx.MultiDiGraph
        Directed influence graph.
        Edges store all columns from the input file as attributes.

    Raises
    ------
    FileNotFoundError
        If `file` does not exist.
    ValueError
        If the input file is missing required columns or contains unsupported
        sign values.
    TypeError
        If `identifiers` is neither None nor a GeneIdentifiers object.

    Notes
    -----
    The 'sign' column must contain:
    - 1 for positive influence
    - -1 for negative influence
    - NaN for unsigned or unknown influence
    """

    file = Path(file)

    if not file.is_file():
        raise FileNotFoundError(f"file not found: {file}")

    influence_graph = pd.read_csv(file, sep=sep, **kwargs)

    required_columns = {"source", "target", "sign"}
    missing_columns = required_columns - set(influence_graph.columns)

    if missing_columns:
        raise ValueError(f"invalid file: missing columns {sorted(missing_columns)}")

    influence_graph["sign"] = influence_graph["sign"].astype(float)

    invalid_signs = set(influence_graph["sign"].dropna().unique()) - {-1.0, 1.0}

    if invalid_signs:
        raise ValueError(
            f"invalid file: unsupported sign values {sorted(invalid_signs)}"
        )

    grn = nx.from_pandas_edgelist(
        df=influence_graph,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph,
    )

    if identifiers is None:
        return grn

    if callable(identifiers):
        identifiers(
            grn,
            input_type=input_type,
            output_type=output_type,
            copy=False,
        )
        return grn

    raise TypeError(
        f"unsupported argument type for 'identifiers': "
        "expected a GeneIdentifiers-compatible callable but received "
        f"{type(identifiers)}"
    )
