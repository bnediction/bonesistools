#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Any
)
from collections.abc import Mapping
from anndata import AnnData
from .._typing import (
    ScData,
    Metric,
    Metric_Function,
    anndata_checker,
    anndata_or_mudata_checker
)

from ._utils import (
    choose_representation
)

import numpy as np
from sklearn import metrics

@anndata_checker
def pairwise_distances(
    adata,
    n_components: Optional[int] = None,
    use_rep: Optional[str] = None,
    metric: Union[Metric, Metric_Function] = "euclidean",
    key_added: Optional[str] = None,
    n_jobs: int = 1,
    **metric_kwds: Mapping[str, Any]
) -> Union[np.ndarray, None]:
    """
    Calculate the distance matrix between observations with respect to a choosen representation.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    n_components: int (optional, default: None)
        Number of principal components or dimensions in the embedding space taken into account for each observation.
        If not specified, consider all principal components.
    use_rep: str (optional, default: None)
        Use the indicated representation in adata.obsm.
    metric: Metric | Metric_Function (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    key_added: str (optional, default: None)
        If not specified, return distance matrix.
        If specified, distance matrix is added to .obsp['key_added'].
    n_jobs: int (default: 1)
        Number of allocated processors.
    **metric_kwds
        Any further parameters passed to the distance function.

    Returns
    -------
    Depending on 'key_added', update 'adata' or return ndarray object.
    """
    
    X = choose_representation(adata, use_rep=use_rep, n_components=n_components)

    distances = metrics.pairwise_distances(
        X,
        metric=metric,
        n_jobs=n_jobs,
        **metric_kwds
    )

    if key_added is not None:
        adata.obsp[key_added] = distances
        adata.uns[key_added] = {
            "n_components": n_components,
            "use_rep": use_rep,
            "metric": metric,
            "metric_kwds": metric_kwds
        }
        return None
    else:
        return distances

@anndata_or_mudata_checker
def barycenters(
    scdata: ScData, # type: ignore
    obs: str,
    use_rep: Optional[str] = None,
    n_components: Optional[int] = None,
) -> Mapping[str, np.ndarray]:
    """
    Calculate the barycenter with respect to an embedding space.

    Parameters
    ----------
    scdata: ad.AnnData | md.MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        The classification is retrieved by .obs['obs'], which must be categorical/qualitative values.
    use_rep: str (optional, default: None)
        Use the indicated representation in scdata.obsm.
    n_components: int (optional, default: None)
        Number of principal components or dimensions in the embedding space
        taken into account for each observation.

    Returns
    -------
    Return dict-like mapping where keys are clusters and values are barycenter values.
    """

    X = choose_representation(scdata, use_rep=use_rep, n_components=n_components)
    
    if not hasattr(scdata.obs[obs], "cat"):
        raise AttributeError(f"scdata.obs[{obs}] object has not attribute '.cat'")
    else:
        clusters = scdata.obs[obs].cat.categories

    return {cluster: np.nanmean(X[scdata.obs[obs] == cluster], axis=0) for cluster in clusters}
