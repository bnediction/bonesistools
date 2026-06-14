#!/usr/bin/env python

from __future__ import annotations

import importlib
import inspect
from typing import Any, Dict, Optional, Tuple, Union, cast

import numpy as np
from anndata import AnnData
from scipy import sparse

from ..._compat import Literal, get_args
from ..._typing import RandomStateSeed
from ..._validation import (
    _as_boolean,
    _as_literal,
    _as_positive_integer,
    _as_positive_number,
    _as_seed,
    _as_string,
)
from .._dependencies import require_sklearn
from .._metadata import _format_random_state
from .._typing import Metric, UMAPMetric, anndata_checker
from ._utils import choose_matrix_representation, choose_representation

EmbeddingMethod = Literal["spectral", "umap", "tsne"]
EigenSolver = Literal["arpack", "lobpcg", "amg"]
PCASolver = Literal["auto", "full", "arpack", "randomized"]
TruncatedSVDSolver = Literal["arpack", "randomized"]

EMBEDDING_METHODS: Tuple[EmbeddingMethod, ...] = ("spectral", "umap", "tsne")
EIGEN_SOLVERS: Tuple[EigenSolver, ...] = ("arpack", "lobpcg", "amg")
PCA_SOLVERS: Tuple[PCASolver, ...] = ("auto", "full", "arpack", "randomized")


@require_sklearn
@anndata_checker
def pca(
    adata: AnnData,
    n_components: int = 50,
    *,
    layer: Optional[str] = None,
    use_raw: bool = False,
    use_highly_variable: Optional[bool] = None,
    zero_center: bool = True,
    svd_solver: PCASolver = "auto",
    key_added: str = "X_pca",
    seed: RandomStateSeed = 0,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute principal components and store them in an AnnData object.

    PCA (Principal Component Analysis) is a linear dimensionality reduction
    method based on the Singular Value Decomposition (SVD).

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    n_components: int (default: 50)
        Number of principal components to compute.
    layer: str, optional
        Layer to use instead of `adata.X`.
    use_raw: bool (default: False)
        Use `adata.raw.X` instead of `adata.X`.
    use_highly_variable: bool, optional
        If True, use only genes marked by `adata.var["highly_variable"]`.
    zero_center: bool (default: True)
        If True, use centered PCA. If False, use truncated SVD.
    svd_solver: {'auto', 'full', 'arpack', 'randomized'} (default: 'auto')
        SVD solver passed to the sklearn estimator.
    key_added: str (default: 'X_pca')
        Key used to store principal component scores in `adata.obsm`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by the sklearn estimator.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Examples
    --------
    >>> bt.sct.tl.pca(adata, n_components=30)

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with PCA results added.
        Otherwise, updates `adata` in place and returns None.

        PCA results are stored in:

        - `adata.obsm[key_added]`: projected cell coordinates;
        - `adata.varm["PCs"]`: principal component loadings;
        - `adata.uns["pca"]`: PCA metadata.

    Notes
    -----
    A fixed `seed` makes the randomized parts deterministic, but it does not
    fully control multi-threaded numerical kernels. When PCA runs through a
    stack using several threads, for example with `n_jobs > 1` in a surrounding
    workflow or with multi-threaded BLAS/OpenMP backends, the order of
    floating-point operations can vary. The resulting differences are usually
    very small, but they may still propagate to neighbor graphs, embeddings,
    or clusters. Use a single-threaded numerical backend when strict
    reproducibility is required.
    """

    from sklearn.decomposition import PCA, TruncatedSVD

    n_components = _as_positive_integer(n_components, "n_components")

    zero_center = _as_boolean(zero_center, "zero_center")

    svd_solver = _as_literal(
        svd_solver,
        choices=PCA_SOLVERS,
        name="svd_solver",
    )

    if not zero_center and svd_solver == "full":
        raise ValueError(
            "invalid argument value for 'svd_solver': "
            "'full' is not supported when zero_center=False"
        )

    use_highly_variable = _as_boolean(
        use_highly_variable,
        "use_highly_variable",
        allow_none=True,
    )

    key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    matrix = cast(
        Any,
        choose_matrix_representation(
            adata,
            use_raw=use_raw,
            layer=layer,
            copy=True,
        ),
    )

    if sparse.issparse(matrix):
        X = cast(np.ndarray, matrix.toarray())
    else:
        X = np.asarray(matrix)
    mask = None
    if use_highly_variable is True:
        if "highly_variable" not in adata.var:
            raise KeyError(
                "column 'highly_variable' not found in adata.var: "
                "please run highly variable gene selection first"
            )

        mask = np.asarray(adata.var["highly_variable"], dtype=bool)
        if not bool(mask.any()):
            raise ValueError("no highly variable genes selected")

    if mask is not None:
        X = X[:, mask]

    max_components = min(X.shape)
    if n_components > max_components:
        raise ValueError(
            f"invalid argument value for 'n_components': "
            f"expected value <= {max_components} but received {n_components!r}"
        )

    resolved_random_state = _as_seed(seed)
    if zero_center:
        estimator = PCA(
            n_components=n_components,
            svd_solver=svd_solver,
            random_state=resolved_random_state,
        )
    else:
        truncated_solver: TruncatedSVDSolver
        if svd_solver == "arpack":
            truncated_solver = "arpack"
        else:
            truncated_solver = "randomized"

        estimator = TruncatedSVD(
            n_components=n_components,
            algorithm=truncated_solver,
            random_state=resolved_random_state,
        )

    scores = cast(np.ndarray, estimator.fit_transform(X))
    components = cast(np.ndarray, estimator.components_).T
    loadings = np.zeros((adata.n_vars, n_components), dtype=components.dtype)
    if mask is None:
        loadings[:, :] = components
    else:
        loadings[mask, :] = components

    adata.obsm[key_added] = scores
    adata.varm["PCs"] = loadings
    adata.uns["pca"] = {
        "params": {
            "zero_center": zero_center,
            "use_highly_variable": use_highly_variable,
            "layer": layer,
            "use_raw": use_raw,
            "svd_solver": svd_solver,
            "key_added": key_added,
            "seed": _format_random_state(seed),
        },
        "variance": cast(np.ndarray, estimator.explained_variance_),
        "variance_ratio": cast(np.ndarray, estimator.explained_variance_ratio_),
    }

    return adata if copy else None


@anndata_checker
def embedding(
    adata: AnnData,
    method: EmbeddingMethod = "umap",
    use_rep: Optional[str] = "X_pca",
    n_rep_components: Optional[int] = None,
    n_components: int = 2,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = None,
    n_jobs: int = 1,
    copy: bool = False,
    **kwargs: Any,
) -> Union[AnnData, None]:
    """
    Compute a low-dimensional embedding from an existing representation.

    `method="umap"` uses `umap-learn` on the selected representation.
    `method="spectral"` computes a spectral embedding with
    `sklearn.manifold.SpectralEmbedding`. `method="tsne"` uses
    `sklearn.manifold.TSNE`.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    method: {'umap', 'spectral', 'tsne'} (default: 'umap')
        Embedding method.
    use_rep: str, optional (default: 'X_pca')
        Representation key in `adata.obsm` used as input.
    n_rep_components: int, optional
        Number of input representation dimensions to use. If None, use all
        dimensions in `use_rep`.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_se`,
        `X_umap` or `X_tsne` depending on `method`.
    seed: int, np.random.RandomState, np.random or None, optional
        Random seed or random state used by the embedding estimator.
    n_jobs: int (default: 1)
        Number of allocated processors.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.
    **kwargs: Any
        Additional keyword arguments passed to the selected embedding function:
        `spectral`, `umap` or `tsne`.

    Examples
    --------
    >>> bt.sct.tl.embedding(adata, use_rep="X_pca")

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with embedding results
        added. Otherwise, updates `adata` in place and returns None.

        Embedding results are stored in:

        - `adata.obsm[key_added]`: cell coordinates;
        - `adata.uns[key_added]`: embedding metadata.
    """

    method = _as_literal(
        method,
        choices=EMBEDDING_METHODS,
        name="method",
    )

    if method == "spectral":
        return spectral(
            adata,
            use_rep=use_rep,
            n_rep_components=n_rep_components,
            n_components=n_components,
            key_added=key_added,
            seed=seed,
            n_jobs=n_jobs,
            copy=copy,
            **kwargs,
        )
    elif method == "umap":
        return umap(
            adata,
            use_rep=use_rep,
            n_rep_components=n_rep_components,
            n_components=n_components,
            key_added=key_added,
            seed=seed,
            n_jobs=n_jobs,
            copy=copy,
            **kwargs,
        )
    else:
        return tsne(
            adata,
            use_rep=use_rep,
            n_rep_components=n_rep_components,
            n_components=n_components,
            key_added=key_added,
            seed=seed,
            n_jobs=n_jobs,
            copy=copy,
            **kwargs,
        )


@require_sklearn
@anndata_checker
def spectral(
    adata: AnnData,
    use_rep: Optional[str] = "X_pca",
    n_rep_components: Optional[int] = None,
    n_components: int = 2,
    n_neighbors: int = 15,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = None,
    eigen_solver: Optional[EigenSolver] = None,
    n_jobs: int = 1,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a spectral embedding from an existing representation.

    Spectral embedding is a nonlinear dimensionality reduction method based on
    the eigenvectors of a graph Laplacian constructed from local neighborhood
    relationships.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    use_rep: str, optional (default: 'X_pca')
        Representation key in `adata.obsm` used as input.
    n_rep_components: int, optional
        Number of input representation dimensions to use. If None, use all
        dimensions in `use_rep`.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    n_neighbors: int (default: 15)
        Number of nearest neighbors used by spectral embedding.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_se`.
    seed: int, np.random.RandomState, np.random or None, optional
        Random seed or random state used by the embedding estimator.
    eigen_solver: {'arpack', 'lobpcg', 'amg'}, optional
        Eigen solver passed to `sklearn.manifold.SpectralEmbedding`.
    n_jobs: int (default: 1)
        Number of allocated processors.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with spectral embedding
        results added. Otherwise, updates `adata` in place and returns None.

        Spectral embedding results are stored in:

        - `adata.obsm[key_added]`: cell coordinates;
        - `adata.uns[key_added]`: embedding metadata.

    References
    ----------
    Belkin and Niyogi (2003). Laplacian Eigenmaps for Dimensionality
    Reduction and Data Representation. Neural Computation, 15(6), 1373-1396.
    """

    from sklearn.manifold import SpectralEmbedding

    n_components = _as_positive_integer(n_components, "n_components")
    n_neighbors = _as_positive_integer(n_neighbors, "n_neighbors")
    if not isinstance(n_jobs, int):
        raise TypeError(
            f"unsupported argument type for 'n_jobs': "
            f"expected {int} but received {type(n_jobs)}"
        )

    resolved_random_state = _as_seed(seed)

    eigen_solver = _as_literal(
        eigen_solver,
        choices=EIGEN_SOLVERS,
        name="eigen_solver",
        allow_none=True,
    )

    if n_neighbors >= adata.n_obs:
        raise ValueError(
            f"invalid argument value for 'n_neighbors': "
            f"expected value smaller than number of observations "
            f"({adata.n_obs}) but received {n_neighbors!r}"
        )

    if key_added is None:
        key_added = "X_se"
    else:
        key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    X = choose_representation(
        adata,
        use_rep=use_rep,
        n_components=n_rep_components,
    )

    embedding = cast(
        np.ndarray,
        SpectralEmbedding(
            n_components=n_components,
            n_neighbors=n_neighbors,
            random_state=resolved_random_state,
            eigen_solver=eigen_solver,
            n_jobs=n_jobs,
        ).fit_transform(X),
    )

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "spectral",
        "use_rep": use_rep,
        "n_rep_components": n_rep_components,
        "n_components": n_components,
        "n_neighbors": n_neighbors,
        "seed": _format_random_state(seed),
        "eigen_solver": eigen_solver,
        "n_jobs": n_jobs,
    }

    return adata if copy else None


@anndata_checker
def umap(
    adata: AnnData,
    use_rep: Optional[str] = "X_pca",
    n_rep_components: Optional[int] = None,
    n_components: int = 2,
    n_neighbors: int = 15,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    min_dist: float = 0.5,
    spread: float = 1.0,
    max_iter: Optional[int] = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: Union[str, np.ndarray] = "spectral",
    a: Optional[float] = None,
    b: Optional[float] = None,
    metric: UMAPMetric = "euclidean",
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a UMAP embedding from an existing representation.

    UMAP is a nonlinear dimensionality reduction method based on neighborhood
    graphs and manifold learning.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    use_rep: str, optional (default: 'X_pca')
        Representation key in `adata.obsm` used as input.
    n_rep_components: int, optional
        Number of input representation dimensions to use. If None, use all
        dimensions in `use_rep`.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    n_neighbors: int (default: 15)
        Number of nearest neighbors used by UMAP.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_umap`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by UMAP.
    n_jobs: int (default: 1)
        Number of allocated processors.
    min_dist: float (default: 0.5)
        Effective minimum distance between embedded points.
    spread: float (default: 1.0)
        Effective scale of embedded points.
    max_iter: int, optional
        Number of optimization epochs.
    alpha: float (default: 1.0)
        Initial learning rate.
    gamma: float (default: 1.0)
        Negative sample weighting.
    negative_sample_rate: int (default: 5)
        Number of negative samples per positive sample.
    init_pos: str or ndarray (default: 'spectral')
        Initialization used by UMAP.
    a: float, optional
        UMAP curve parameter. If None, UMAP determines it from `min_dist` and
        `spread`.
    b: float, optional
        UMAP curve parameter. If None, UMAP determines it from `min_dist` and
        `spread`.
    metric: UMAPMetric (default: 'euclidean')
        Distance metric.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with UMAP embedding results
        added. Otherwise, updates `adata` in place and returns None.

        UMAP embedding results are stored in:

        - `adata.obsm[key_added]`: cell coordinates;
        - `adata.uns[key_added]`: embedding metadata.

    References
    ----------
    McInnes et al. (2018). UMAP: Uniform Manifold Approximation and Projection
    for Dimension Reduction. arXiv:1802.03426.
    """

    n_components = _as_positive_integer(n_components, "n_components")
    n_neighbors = _as_positive_integer(n_neighbors, "n_neighbors")
    if not isinstance(n_jobs, int):
        raise TypeError(
            f"unsupported argument type for 'n_jobs': "
            f"expected {int} but received {type(n_jobs)}"
        )

    resolved_random_state = _as_seed(seed)

    if n_neighbors >= adata.n_obs:
        raise ValueError(
            f"invalid argument value for 'n_neighbors': "
            f"expected value smaller than number of observations "
            f"({adata.n_obs}) but received {n_neighbors!r}"
        )

    if max_iter is not None:
        max_iter = _as_positive_integer(max_iter, "max_iter")

    negative_sample_rate = _as_positive_integer(
        negative_sample_rate,
        "negative_sample_rate",
    )

    metric = _as_literal(
        metric,
        choices=get_args(UMAPMetric),
        name="metric",
    )

    if key_added is None:
        key_added = "X_umap"
    else:
        key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    X = choose_representation(
        adata,
        use_rep=use_rep,
        n_components=n_rep_components,
    )

    try:
        umap_module = importlib.import_module("umap")
    except ImportError as error:
        raise ImportError(
            "umap-learn is required for `bt.sct.tl.umap`. "
            "Install bonesistools with the sctools extra or install umap-learn."
        ) from error

    UMAP = cast(Any, getattr(umap_module, "UMAP"))
    estimator = cast(
        Any,
        UMAP(
            n_neighbors=n_neighbors,
            n_components=n_components,
            metric=metric,
            min_dist=min_dist,
            spread=spread,
            n_epochs=max_iter,
            learning_rate=alpha,
            repulsion_strength=gamma,
            negative_sample_rate=negative_sample_rate,
            init=init_pos,
            a=a,
            b=b,
            random_state=resolved_random_state,
            n_jobs=n_jobs,
        ),
    )
    embedding = cast(np.ndarray, estimator.fit_transform(X))

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "umap",
        "use_rep": use_rep,
        "n_rep_components": n_rep_components,
        "n_components": n_components,
        "seed": _format_random_state(seed),
        "n_neighbors": n_neighbors,
        "metric": metric,
        "min_dist": min_dist,
        "spread": spread,
        "max_iter": max_iter,
        "alpha": alpha,
        "gamma": gamma,
        "negative_sample_rate": negative_sample_rate,
        "init_pos": init_pos,
        "a": a,
        "b": b,
        "n_jobs": n_jobs,
    }

    return adata if copy else None


@require_sklearn
@anndata_checker
def tsne(
    adata: AnnData,
    use_rep: Optional[str] = "X_pca",
    n_rep_components: Optional[int] = None,
    n_components: int = 2,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    max_iter: Optional[int] = None,
    perplexity: float = 30.0,
    learning_rate: float = 1000.0,
    early_exaggeration: float = 12.0,
    metric: Metric = "euclidean",
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a t-SNE embedding from an existing representation.

    t-SNE (t-distributed Stochastic Neighbor Embedding) is a nonlinear
    dimensionality reduction method that preserves local neighborhood
    structure in a low-dimensional embedding.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    use_rep: str, optional (default: 'X_pca')
        Representation key in `adata.obsm` used as input.
    n_rep_components: int, optional
        Number of input representation dimensions to use. If None, use all
        dimensions in `use_rep`.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_tsne`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by t-SNE.
    n_jobs: int (default: 1)
        Number of allocated processors.
    max_iter: int, optional
        Maximum number of optimization iterations. Defaults to 1000.
    perplexity: float (default: 30.0)
        t-SNE perplexity.
    learning_rate: float (default: 1000.0)
        t-SNE learning rate.
    early_exaggeration: float (default: 12.0)
        t-SNE early exaggeration.
    metric: Metric (default: 'euclidean')
        Distance metric.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with t-SNE embedding results
        added. Otherwise, updates `adata` in place and returns None.

        t-SNE embedding results are stored in:

        - `adata.obsm[key_added]`: cell coordinates;
        - `adata.uns[key_added]`: embedding metadata.

    References
    ----------
    van der Maaten and Hinton (2008). Visualizing Data Using t-SNE. Journal
    of Machine Learning Research, 9(11).
    """

    from sklearn.manifold import TSNE

    n_components = _as_positive_integer(n_components, "n_components")
    if not isinstance(n_jobs, int):
        raise TypeError(
            f"unsupported argument type for 'n_jobs': "
            f"expected {int} but received {type(n_jobs)}"
        )

    resolved_random_state = _as_seed(seed)
    tsne_max_iter = 1000 if max_iter is None else max_iter
    tsne_max_iter = _as_positive_integer(tsne_max_iter, "max_iter")

    perplexity = _as_positive_number(perplexity, "perplexity")

    if perplexity >= adata.n_obs:
        raise ValueError(
            f"invalid argument value for 'perplexity': "
            f"expected value smaller than number of observations "
            f"({adata.n_obs}) but received {perplexity!r}"
        )

    metric = _as_literal(
        metric,
        choices=get_args(Metric),
        name="metric",
    )

    if key_added is None:
        key_added = "X_tsne"
    else:
        key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    X = choose_representation(
        adata,
        use_rep=use_rep,
        n_components=n_rep_components,
    )

    tsne_kwargs: Dict[str, Any] = {
        "n_components": n_components,
        "perplexity": perplexity,
        "early_exaggeration": early_exaggeration,
        "learning_rate": learning_rate,
        "metric": metric,
        "random_state": resolved_random_state,
        "n_jobs": n_jobs,
        "method": "barnes_hut" if n_components <= 3 else "exact",
    }
    if "max_iter" in inspect.signature(TSNE).parameters:
        tsne_kwargs["max_iter"] = tsne_max_iter
    else:
        tsne_kwargs["n_iter"] = tsne_max_iter

    embedding = cast(np.ndarray, TSNE(**tsne_kwargs).fit_transform(X))

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "tsne",
        "use_rep": use_rep,
        "n_rep_components": n_rep_components,
        "n_components": n_components,
        "seed": _format_random_state(seed),
        "perplexity": perplexity,
        "metric": metric,
        "early_exaggeration": early_exaggeration,
        "learning_rate": learning_rate,
        "max_iter": tsne_max_iter,
        "n_jobs": n_jobs,
    }

    return adata if copy else None
