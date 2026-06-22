#!/usr/bin/env python

from __future__ import annotations

from typing import (
    Any,
    Mapping,
    Optional,
    Union,
    cast,
    overload,
)

import numpy as np

from .._typing import (
    AnnData,
    Metric,
    Metric_Function,
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
)
from ._utils import (
    _UNSET,
    _resolve_representation_argument,
    get_representation,
)


@overload
def pairwise_distances(
    adata: AnnData,
    n_components: Optional[int] = None,
    representation: Optional[str] = None,
    metric: Union[Metric, Metric_Function] = "euclidean",
    *,
    key_added: None = None,
    n_jobs: int = 1,
    use_rep: Any = _UNSET,
    **metric_kwds: Any,
) -> np.ndarray: ...


@overload
def pairwise_distances(
    adata: AnnData,
    n_components: Optional[int] = None,
    representation: Optional[str] = None,
    metric: Union[Metric, Metric_Function] = "euclidean",
    *,
    key_added: str,
    n_jobs: int = 1,
    use_rep: Any = _UNSET,
    **metric_kwds: Any,
) -> None: ...


@anndata_checker
def pairwise_distances(
    adata: AnnData,  # type: ignore
    n_components: Optional[int] = None,
    representation: Any = _UNSET,
    metric: Union[Metric, Metric_Function] = "euclidean",
    key_added: Optional[str] = None,
    n_jobs: int = 1,
    use_rep: Any = _UNSET,
    **metric_kwds: Any,
) -> Union[np.ndarray, None]:
    """
    Calculate pairwise distances between observations in a selected representation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.
    representation: str, optional
        Representation key in `adata.obsm`.
    metric: Metric | Metric_Function (default: 'euclidean')
        Metric used when calculating pairwise distances between observations.
    key_added: str, optional
        If specified, store the distance matrix in `adata.obsp[key_added]`.
        Otherwise, return the distance matrix.
    n_jobs: int (default: 1)
        Number of allocated processors.
    **metric_kwds: Any
        Additional keyword arguments passed to the distance function.
    use_rep: str, optional
        Deprecated alias for `representation`.

    Returns
    -------
    ndarray or None
        If `key_added` is None, returns the pairwise distance matrix.
        Otherwise, updates `adata` in place and returns None.

        Pairwise distance results are stored in:

        - `adata.obsp[key_added]`: pairwise distance matrix;
        - `adata.uns[key_added]`: distance metadata.
    """

    from sklearn.metrics import pairwise_distances

    representation = _resolve_representation_argument(
        representation,
        use_rep,
        default=None,
        stacklevel=2,
    )
    representation_mtx = get_representation(
        adata,
        obsm=representation,
        n_components=n_components,
    )

    distances = pairwise_distances(
        representation_mtx,
        metric=cast(Any, metric),
        n_jobs=n_jobs,
        **metric_kwds,
    )

    if key_added is not None:
        adata.obsp[key_added] = distances
        adata.uns[key_added] = {
            "n_components": n_components,
            "representation": representation,
            "metric": metric,
            "metric_kwds": metric_kwds,
        }
        return None
    else:
        return distances


@anndata_or_mudata_checker
def barycenters(
    scdata: ScData,  # type: ignore
    obs: str,
    representation: Any = _UNSET,
    n_components: Optional[int] = None,
    *,
    use_rep: Any = _UNSET,
) -> Mapping[str, np.ndarray]:
    """
    Calculate cluster barycenters in an embedding space.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Categorical observation key in `scdata.obs`.
    representation: str, optional
        Representation key in `scdata.obsm`.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.
    use_rep: str, optional
        Deprecated alias for `representation`.

    Returns
    -------
    Mapping[str, ndarray]
        Mapping from cluster names to barycenter coordinates.

    Raises
    ------
    AttributeError
        If `scdata.obs[obs]` does not expose categorical `.cat` access.
    """

    representation = _resolve_representation_argument(
        representation,
        use_rep,
        default=None,
        stacklevel=2,
    )
    representation_mtx = get_representation(
        scdata,
        obsm=representation,
        n_components=n_components,
    )

    if not hasattr(scdata.obs[obs], "cat"):
        raise AttributeError(f"scdata.obs[{obs!r}] object has no attribute 'cat'")
    else:
        clusters = scdata.obs[obs].cat.categories

    return {
        cluster: np.nanmean(
            cast(Any, representation_mtx)[scdata.obs[obs] == cluster],
            axis=0,
        )
        for cluster in clusters
    }
