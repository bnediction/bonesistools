#!/usr/bin/env python

from typing import Any

from typing_extensions import Literal, Protocol

InputIdentifierType = Literal["name", "gene_id", "ensembl_id"]
OutputIdentifierType = Literal[
    "symbol",
    "ncbi_symbol",
    "gene_id",
    "ensembl_id",
    "official_name",
    "ncbi_name",
]


class GeneSynonymsLike(Protocol):
    """
    Structural interface for objects that convert gene identifiers in place.
    """

    def __call__(self, *args: Any, **kwargs: Any) -> Any: ...
