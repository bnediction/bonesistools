#!/usr/bin/env python

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

InputIdentifierType = Literal["name", "gene_id", "ensembl_id"]
OutputIdentifierType = Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"]
