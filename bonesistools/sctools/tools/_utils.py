#!/usr/bin/env python

from typing import Optional
from anndata import AnnData
from numpy import ndarray

from .._typing import (
    ScData,
    anndata_or_mudata_checker
)

@anndata_or_mudata_checker
def choose_representation(
    scdata: ScData, # type: ignore
    use_rep: Optional[str]="X_pca",
    n_components: Optional[str]=None
) -> ndarray:
    """
    Get a matrix representation of data.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    use_rep
        Use the indicated representation in scdata.obsm
    n_components
        Number of principal components or dimensions in the embedding space
        taken into account for each observation

    Returns
    -------
    Return Array object.
    """

    if use_rep is None:
        use_rep = "X_pca"

    if use_rep not in scdata.obsm_keys():
        if use_rep == "X_pca":
            raise KeyError("key 'X_pca' not found in scdata.obsm: please run scanpy.tl.pca")
        else:
            raise KeyError(f"key '{use_rep}' not found in scdata.obsm")

    if n_components is None:
        return scdata.obsm[use_rep]
    else:
        return scdata.obsm[use_rep][:,:n_components]

@anndata_or_mudata_checker
def _get_distances(
    scdata: ScData, # type: ignore
    obsp: Optional[str]=None,
    neighbors_key: Optional[str]=None
) -> ndarray:
    """
    Get distance matrix already computed.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    obsp
        If specified, return scdata.obsp['obsp']
    neighbors_key
        If specified, retrieve distance matrix using information in scdata.uns['neighbors_key']

    Returns
    -------
    Return Array object.
    """

    if obsp is not None and neighbors_key is not None:
        raise ValueError("arguments 'obsp' and 'neighbors_key' cannot be both specified")
    elif obsp is not None:
        return scdata.obsp[obsp]
    elif neighbors_key is not None:
        distances_key = scdata.uns[neighbors_key]["distances_key"]
        return scdata.obsp[distances_key]
    else:
        if "neighbors" in scdata.uns:
            distances_key = scdata.uns["neighbors"]["distances_key"]
            return scdata.obsp[distances_key]
        else:
            raise KeyError("distances not found in 'scdata': please run scanpy.pp.neighbors or specify 'obsp' or 'neighbors_key'")

@anndata_or_mudata_checker
def _get_connectivities(
    scdata: ScData, # type: ignore
    obsp: Optional[str]=None,
    neighbors_key: Optional[str]=None
) -> ndarray:
    """
    Get connectivity matrix already computed.

    Parameters
    ----------
    scdata
        AnnData or MuData object
    obsp
        If specified, return scdata.obsp['obsp']
    neighbors_key
        If specified, retrieve connectivity matrix using information in scdata.uns['neighbors_key']

    Returns
    -------
    Return Array object.
    """

    if obsp is not None and neighbors_key is not None:
        raise ValueError("arguments 'obsp' and 'neighbors_key' cannot be both specified")
    elif obsp is not None:
        return scdata.obsp[obsp]
    elif neighbors_key is not None:
        connectivities_key = scdata.uns[neighbors_key]["connectivities_key"]
        return scdata.obsp[connectivities_key]
    else:
        if "neighbors" in scdata.uns:
            connectivities_key = scdata.uns["neighbors"]["connectivities_key"]
            return scdata.obsp[connectivities_key]
        else:
            raise KeyError("connectivities not found in 'scdata': please run scanpy.pp.neighbors or specify 'obsp' or 'neighbors_key'")
