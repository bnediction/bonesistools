#!/usr/bin/env python

from typing import Optional

from anndata import AnnData

import numpy as np

from .._check_anndata import _adata_arg_checking

@_adata_arg_checking
def barycenters(
    adata: AnnData,
    obs: str,
    obsm: str="X_umap",
    n_components: Optional[int]=None
):
    """
    Compute the barycenter with respect to an embedding projection in adata.obsm.

    Parameters
    ----------
    adata
        Annotated data matrix
    obs
        The classification is retrieved by .obs[`obs`], which must be categorical/qualitative values
    obsm
        The data points are retrieved by the first `n_components` rows in .obsm[`obsm`]
    n_components
        Dimension of the embedding projection (default: all components)

    Returns
    -------
    dict[cluster] = barycenter
    """

    if n_components is None:
        n_components = adata.obsm[obsm].shape[1]
    
    if not hasattr(adata.obs[obs], "cat"):
        raise TypeError("adata.obs[`obs`] has not attribute `.cat`")
    else:
        clusters = adata.obs[obs].cat.categories

    return {cluster: np.nanmean(adata.obsm[obsm][:,:n_components][adata.obs[obs] == cluster], axis=0) for cluster in clusters}
