#!/usr/bin/env python

from typing import Optional

from anndata import AnnData

import pandas as pd
import numpy as np

from . import barycenters
from .._check_anndata import _adata_arg_checking

@_adata_arg_checking
def subclusters_at_extremity(
    adata: AnnData,
    obs: str,
    obsm: str="X_umap",
    n_components: Optional[int]=None,
    key: str="subclusters",
    n_neighbors: int=30,
    copy: bool=False
):
    
    adata = adata.copy() if copy else adata
    
    df = pd.DataFrame(index=adata.obs.index, columns=[key])

    _barycenters = barycenters(
        adata,
        obs=obs,
        obsm=obsm,
        n_components=n_components
    )

    for cluster in adata.obs[obs].cat.categories:
        barycenters_wo_cluster = _barycenters.copy(); del barycenters_wo_cluster[cluster]
        distances = []
        cluster_data_points = adata.obsm[obsm][:,:n_components][adata.obs[obs] == cluster]
        for data_point in cluster_data_points:
            distances_from_barycenters = []
            for value in barycenters_wo_cluster.values():
                distances_from_barycenters.append(np.linalg.norm(value-data_point))
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
