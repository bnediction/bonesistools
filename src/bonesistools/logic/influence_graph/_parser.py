#!/usr/bin/env python

"""
Deprecated parser module kept for internal compatibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import networkx as nx

from ..._warnings import _rename_deprecated_arguments
from ...resources.ncbi._typing import (
    GeneIdentifiersLike,
    InputIdentifierType,
    OutputIdentifierType,
)
from ..input_output import _influence_graph


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
    return _influence_graph.read_influence_graph(
        file,
        identifiers=identifiers,
        input_type=input_type,
        output_type=output_type,
        sep=sep,
        **kwargs,
    )


__all__ = ["read_influence_graph"]
