#!/usr/bin/env python

"""
Deprecated parser module kept for internal compatibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import networkx as nx

from ...resources.ncbi._typing import InputIdentifierType, OutputIdentifierType
from ..input_output import _influence_graph

GeneSynonyms = _influence_graph.GeneSynonyms


def read_influence_graph(
    infile: Union[str, Path],
    genesyn: Optional[Any] = None,
    input_identifier_type: InputIdentifierType = "name",
    output_identifier_type: OutputIdentifierType = "official_name",
    sep: str = ",",
    **kwargs: Any,
) -> nx.MultiDiGraph[Any]:
    original_gene_synonyms = _influence_graph.GeneSynonyms
    _influence_graph.GeneSynonyms = GeneSynonyms

    try:
        return _influence_graph.read_influence_graph(
            infile,
            genesyn=genesyn,
            input_identifier_type=input_identifier_type,
            output_identifier_type=output_identifier_type,
            sep=sep,
            **kwargs,
        )
    finally:
        _influence_graph.GeneSynonyms = original_gene_synonyms


__all__ = ["read_influence_graph"]
