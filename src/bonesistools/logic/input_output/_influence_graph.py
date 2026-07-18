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

from ...resources.ncbi._genesyn import support_legacy_gene_synonyms_args
from ...resources.ncbi._typing import (
    GeneIdentifiersLike,
    InputIdentifierType,
    OutputIdentifierType,
)


@support_legacy_gene_synonyms_args
def read_influence_graph(
    file: Union[str, Path],
    *,
    genesyn: Optional[GeneIdentifiersLike] = None,
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
    genesyn: GeneIdentifiers, optional
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
        If `genesyn` is neither None nor a GeneIdentifiers object.

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

    if genesyn is None:
        return grn

    if callable(genesyn):
        genesyn(
            grn,
            input_type=input_type,
            output_type=output_type,
            copy=False,
        )
        return grn

    raise TypeError(
        f"unsupported argument type for 'genesyn': "
        f"expected a GeneIdentifiers-compatible callable but received {type(genesyn)}"
    )
