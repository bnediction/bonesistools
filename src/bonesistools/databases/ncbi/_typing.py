#!/usr/bin/env python

from ..._compat import Literal

InputIdentifierType = Literal["name", "gene_id", "ensembl_id"]
OutputIdentifierType = Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"]
