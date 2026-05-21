#!/usr/bin/env python

from pathlib import Path
from typing import Any, Mapping, Optional, Union

import networkx as nx
import pandas as pd

from ...databases.ncbi._genesyn import GeneSynonyms
from ...databases.ncbi._typing import (
    InputIdentifierType,
    OutputIdentifierType,
)


def read_interaction_graph(
    infile: Union[str, Path],
    genesyn: Optional[GeneSynonyms] = None,
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    sep: str = ",",
    **kwargs: Mapping[str, Any],
) -> nx.MultiDiGraph:
    """
    Read an interaction graph from a tabular file.

    The input file must contain at least three columns: 'source', 'target' and
    'sign'. Additional columns are kept as edge attributes.

    Parameters
    ----------
    infile: str | Path
        Path to the tabular file containing the interaction graph.
    genesyn: GeneSynonyms (optional, default: None)
        If specified, convert node identifiers using a GeneSynonyms object.
    input_identifier_type: str (default: 'name')
        Gene identifier input format.
    output_identifier_type: str (default: 'official_name')
        Gene identifier output format.
    sep: str (optional, default: None)
        Field delimiter passed to pandas.read_csv.
        If None, pandas infers the delimiter.
    **kwargs: Mapping[str, Any]
        Additional keyword arguments passed to pandas.read_csv.

    Returns
    -------
    nx.MultiDiGraph
        Directed interaction graph.
        Edges store all columns from the input file as attributes.

    Raises
    ------
    FileNotFoundError
        If `infile` does not exist.
    ValueError
        If the input file is missing required columns or contains unsupported
        sign values.
    TypeError
        If `genesyn` is neither None nor a GeneSynonyms object.

    Notes
    -----
    The 'sign' column must contain:
    - 1 for positive influence
    - -1 for negative influence
    - NaN for unsigned or unknown influence
    """

    infile = Path(infile)

    if not infile.is_file():
        raise FileNotFoundError(f"file not found: {infile}")

    interaction_graph = pd.read_csv(infile, sep=sep, **kwargs)

    required_columns = {"source", "target", "sign"}
    missing_columns = required_columns - set(interaction_graph.columns)

    if missing_columns:
        raise ValueError(f"invalid file: missing columns {sorted(missing_columns)}")

    interaction_graph["sign"] = interaction_graph["sign"].astype(float)

    invalid_signs = set(interaction_graph["sign"].dropna().unique()) - {-1.0, 1.0}

    if invalid_signs:
        raise ValueError(
            f"invalid file: unsupported sign values {sorted(invalid_signs)}"
        )

    grn = nx.from_pandas_edgelist(
        df=interaction_graph,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph,
    )

    if genesyn is None:
        return grn

    if isinstance(genesyn, GeneSynonyms):
        genesyn(
            grn,
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            copy=False,
        )
        return grn

    raise TypeError(
        f"unsupported argument type for 'genesyn': "
        f"expected {GeneSynonyms} but received {type(genesyn)}"
    )
