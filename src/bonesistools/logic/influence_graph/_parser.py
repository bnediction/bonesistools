#!/usr/bin/env python

"""
Deprecated parser module kept for internal compatibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import networkx as nx

from ...resources.ncbi._genesyn import support_legacy_gene_synonyms_args
from ...resources.ncbi._typing import InputIdentifierType, OutputIdentifierType
from ..input_output import _influence_graph


@support_legacy_gene_synonyms_args
def read_influence_graph(
    file: Union[str, Path],
    *,
    genesyn: Optional[Any] = None,
    input_type: InputIdentifierType = "name",
    output_type: OutputIdentifierType = "symbol",
    sep: str = ",",
    **kwargs: Any,
) -> nx.MultiDiGraph[Any]:
    return _influence_graph.read_influence_graph(
        file,
        genesyn=genesyn,
        input_type=input_type,
        output_type=output_type,
        sep=sep,
        **kwargs,
    )


__all__ = ["read_influence_graph"]
