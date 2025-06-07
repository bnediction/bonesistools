#!/usr/bin/env python

from ._simple import (
    set_index,
    filter_obs,
    filter_var,
    regress_out,
    merge,
    transfer_layer,
    transfer_obs_sti,
    transfer_obs_its
)

from ._genename import (
    gene_synonyms_conversion,
    set_ncbi_reference_name,
    var_names_merge_duplicates
)
