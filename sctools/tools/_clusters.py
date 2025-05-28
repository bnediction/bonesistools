#!/usr/bin/env python

import collections.abc
from typing import Optional, Union, Sequence
from .._typing import anndata_checker

from anndata import AnnData

import pandas as pd
import numpy as np

from . import barycenters

@anndata_checker
def subclusters_at_center(
    adata: AnnData,
    obs: str,
    obsm: str="X_umap",
    n_components: Optional[int]=None,
    key: str="subclusters",
    n_neighbors: int=30,
    include: Optional[Sequence[str]]=None,
    copy: bool=False
) -> Union[AnnData, None]:
    """
    Compute subclusters at center with respect to an embedding projection in adata.obsm.

    Parameters
    ----------
    adata
        Annotated data matrix
    obs
        The classification is retrieved by .obs[`obs`], which must be categorical/qualitative values
    obsm
        The data points are retrieved by the first `n_components` rows in .obsm[`obsm`]
    n_components
        Number of principal components or dimensions in the embedding space
        taken into account for each observation
    key
        Name of the column created in adata.obs
    n_neighbors
        Size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation
    include
        Cluster for which subclusters are computed (default: all clusters)
    copy
        return a copy instead of updating `adata` object

    Returns
    -------
    Depending on `copy`, updates or returns `adata`.
    """
    
    adata = adata.copy() if copy else adata

    if include is None:
        include = adata.obs[obs].cat.categories
    elif isinstance(include, collections.abc.MutableSequence) or isinstance(include, set):
        pass
    else:
        raise TypeError(f"`include` must be a sequence (value: {include}).")
    
    df = pd.DataFrame(index=adata.obs.index, columns=[key])

    _barycenters = barycenters(
        adata,
        obs=obs,
        obsm=obsm,
        n_components=n_components
    )

    for cluster in adata.obs[obs].cat.categories:
        if cluster not in include:
            continue
        cluster_data_points = adata.obsm[obsm][:,:n_components][adata.obs[obs] == cluster]
        if n_neighbors > cluster_data_points.shape[0]:
            barcodes = adata.obs[obs][adata.obs[obs] == cluster].index
        else:
            distances_from_barycenter = np.linalg.norm(cluster_data_points - _barycenters[cluster], axis=1)
            idx = distances_from_barycenter.argsort()[:n_neighbors]
            barcodes = adata.obs[obs][adata.obs[obs] == cluster].index[idx]
        for barcode in barcodes:
            df.at[barcode,key] = cluster
    
    adata.obs[key] = df.astype("category")

    return adata if copy else None

@anndata_checker
def subclusters_at_extremity(
    adata: AnnData,
    obs: str,
    obsm: str="X_umap",
    n_components: Optional[int]=None,
    key: str="subclusters",
    n_neighbors: int=30,
    include: Optional[Sequence[str]]=None,
    exclude_for_computation: Optional[Sequence[str]]=None,
    copy: bool=False
) -> Union[AnnData, None]:
    """
    Compute subclusters at end with respect to an embedding projection in adata.obsm.

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
    key
        Name of the column created in adata.obs
    n_neighbors
        Size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation
    include
        Cluster for which subclusters are computed (default: all clusters)
    exclude_for_computation
        Clusters not taken into account for computing other subclusters (default: None)
    copy
        return a copy instead of updating `adata` object

    Returns
    -------
    Depending on `copy`, updates or returns `adata`.
    """
    
    adata = adata.copy() if copy else adata

    if include is None:
        include = []
    elif isinstance(include, collections.abc.MutableSequence) or isinstance(include, set):
        pass
    else:
        raise TypeError(f"`include` must be a sequence (value: {include}).")
    
    df = pd.DataFrame(index=adata.obs.index, columns=[key])

    _barycenters = barycenters(
        adata,
        obs=obs,
        obsm=obsm,
        n_components=n_components
    )

    if exclude_for_computation is None:
        pass
    elif isinstance(exclude_for_computation, collections.abc.MutableSequence) or isinstance(exclude_for_computation, set):
        exclude_for_computation = set(exclude_for_computation) if isinstance(exclude_for_computation, collections.abc.MutableSequence) else exclude_for_computation
        for cluster in exclude_for_computation:
            del _barycenters[cluster]
    else:
        raise TypeError(f"`exclude_for_computation` must be a sequence (value: {exclude_for_computation}).")

    for cluster in adata.obs[obs].cat.categories:
        if cluster not in include:
            continue
        barycenters_wo_cluster = _barycenters.copy(); del barycenters_wo_cluster[cluster]
        distances = []
        cluster_data_points = adata.obsm[obsm][:,:n_components][adata.obs[obs] == cluster]
        if n_neighbors > cluster_data_points.shape[0]:
            barcodes = adata.obs[obs][adata.obs[obs] == cluster].index
        else:
            for data_point in cluster_data_points:
                distances_from_barycenters = []
                for value in barycenters_wo_cluster.values():
                    distances_from_barycenters.append(np.linalg.norm(data_point-value))
                distances.append(min(distances_from_barycenters))
            selected_cell_idx = max(range(len(distances)), key=distances.__getitem__)
            selected_cell_data_point = cluster_data_points[selected_cell_idx,:]
            distances_from_selected_cell = np.linalg.norm(cluster_data_points - selected_cell_data_point, axis=1)
            idx = distances_from_selected_cell.argsort()[:n_neighbors]
            barcodes = adata.obs[obs][adata.obs[obs] == cluster].index[idx]
        for barcode in barcodes:
            df.at[barcode,key] = cluster
    
    adata.obs[key] = df.astype("category")

    return adata if copy else None

@anndata_checker
def subclusters(
    adata: AnnData,
    obs: str,
    obsm: str="X_umap",
    n_components: Optional[int]=None,
    key: str="subclusters",
    n_neighbors: int=30,
    include_center: Optional[Sequence[str]]=None,
    include_extremity: Optional[Sequence[str]]=None,
    exclude_for_computation: Optional[Sequence[str]]=None,
    copy: bool=False
) -> Union[AnnData, None]:
    """
    Compute subclusters at center with respect to an embedding projection in adata.obsm.

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
    key
        Name of the column created in adata.obs
    n_neighbors
        Size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation
    include_center
        Cluster for which subclusters are computed with subclusters_at_center() (default: no cluster)
    include_extremity
        Cluster for which subclusters are computed with subclusters_at_extremity() (default: no cluster)
    exclude_for_computation
        Clusters not taken into account for computing other subclusters (default: None)
    copy
        return a copy instead of updating `adata` object

    Returns
    -------
    Depending on `copy`, updates or returns `adata`.
    """
    
    adata = adata.copy() if copy else adata

    if include_center is None and include_extremity is None:
        raise ValueError(f"at least one of the arguments `include_center` or `include_extremity` must be specified")
    elif include_center is not None and include_extremity is not None:
        include_center = set(include_center)
        include_extremity = set(include_extremity)
        if bool(include_center.intersection(include_extremity)):
            raise ValueError(f"`include_center` and `include_extremity` must not contains common values")
        else:
            key_center = key + "#0"
            key_extremity = key + "#1"
            subclusters_at_center(
                adata=adata,
                obs=obs,
                obsm=obsm,
                n_components=n_components,
                key=key_center,
                n_neighbors=n_neighbors,
                include=include_center,
                copy=False
            )
            subclusters_at_extremity(
                adata=adata,
                obs=obs,
                obsm=obsm,
                n_components=n_components,
                key=key_extremity,
                n_neighbors=n_neighbors,
                include=include_extremity,
                exclude_for_computation=exclude_for_computation,
                copy=False
            )
            subclusters_df = pd.concat([adata.obs[key_center].astype(str).replace("nan","",regex=True),adata.obs[key_extremity].astype(str).replace("nan","",regex=True)],axis=1)
            subclusters_df = subclusters_df[key_center] + subclusters_df[key_extremity]
            subclusters_df = subclusters_df.replace(r'^\s*$', np.nan, regex=True).astype("category")
            subclusters_df.name = key
            adata.obs = pd.concat([adata.obs, subclusters_df], axis=1)
            del adata.obs[key_center]; del adata.obs[key_extremity]
    elif include_center is not None and include_extremity is None:
        include_center = set(include_center)
        subclusters_at_center(
            adata=adata,
            obs=obs,
            obsm=obsm,
            n_components=n_components,
            key=key,
            n_neighbors=n_neighbors,
            include=include_center,
            copy=False
        )
    elif include_center is None and include_extremity is not None:
        include_extremity = set(include_extremity)
        subclusters_at_extremity(
            adata=adata,
            obs=obs,
            obsm=obsm,
            n_components=n_components,
            key=key,
            n_neighbors=n_neighbors,
            include=include_extremity,
            copy=False
        )
    
    return adata if copy else None