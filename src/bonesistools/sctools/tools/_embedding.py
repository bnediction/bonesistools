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
from .._metadata import _format_random_state
from .._typing import Metric, VarSubset, anndata_checker
from .._validation import _as_var_subset
from ._utils import get_expression, get_representation

EigenSolver = Literal["arpack", "lobpcg", "amg"]
PCASolver = Literal["auto", "full", "arpack", "randomized"]
TruncatedSVDSolver = Literal["arpack", "randomized"]

EIGEN_SOLVERS: Tuple[EigenSolver, ...] = ("arpack", "lobpcg", "amg")
CENTERED_PCA_SOLVERS: Tuple[PCASolver, ...] = (
    "auto",
    "full",
    "arpack",
    "randomized",
)
UNCENTERED_PCA_SOLVERS: Tuple[TruncatedSVDSolver, ...] = ("arpack", "randomized")


def _format_var_subset(var_subset: VarSubset) -> Union[str, Tuple[str, ...], None]:

    if var_subset is None or isinstance(var_subset, str):
        return var_subset
    return tuple(sorted(var_subset))


@anndata_checker
def pca(
    adata: AnnData,
    n_components: int = 50,
    *,
    layer: Optional[str] = None,
    use_raw: bool = False,
    var_subset: VarSubset = None,
    zero_center: bool = True,
    svd_solver: PCASolver = "auto",
    key_added: str = "X_pca",
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a PCA representation of observations and store the results in an
    AnnData object.

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
    var_subset: str or collection of str, optional
        Variables used for PCA. If a string is provided, it is interpreted as a
        Boolean column in `adata.var`. If a collection is provided,
        it is interpreted as a collection of variable names.
    zero_center: bool (default: True)
        If True, use centered PCA. If False, use truncated SVD.
    svd_solver: {'auto', 'full', 'arpack', 'randomized'} (default: 'auto')
        SVD solver passed to the sklearn estimator. If `'auto'`, use
        `'arpack'` when `zero_center=True` and `'randomized'` when
        `zero_center=False`.
    key_added: str (default: 'X_pca')
        Key used to store PCA coordinates in `adata.obsm`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by the sklearn estimator.
    n_jobs: int (default: 1)
        Number of threads allocated to the PCA computation.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Examples
    --------
    Compute PCA using the variables marked as highly variable:

    >>> bt.sct.tl.pca(
    ...     adata,
    ...     n_components=30,
    ...     var_subset="highly_variable",
    ... )

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with PCA results added.
        Otherwise, updates `adata` in place and returns None.

        PCA results are stored in:

        - `adata.obsm[key_added]`: projected cell coordinates;
        - `adata.varm["PCs"]`: principal component loadings;
        - `adata.uns["pca"]`: PCA metadata.

    References
    ----------
    Jolliffe (1986). Principal Components in Regression Analysis. Principal
    Component Analysis. Springer, 129-155.
    """

    from sklearn.decomposition import PCA, TruncatedSVD
    from threadpoolctl import threadpool_limits

    n_components = _as_positive_integer(n_components, "n_components")
    n_jobs = _as_positive_integer(n_jobs, "n_jobs")

    zero_center = _as_boolean(zero_center, "zero_center")

    resolved_svd_solver = _as_literal(
        svd_solver,
        choices=CENTERED_PCA_SOLVERS,
        name="svd_solver",
    )

    if not zero_center and resolved_svd_solver == "full":
        raise ValueError(
            "invalid argument value for 'svd_solver': "
            "'full' is not supported when zero_center=False"
        )

    key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    mask = _as_var_subset(adata, var_subset)
    selected_expression_mtx = cast(
        Any,
        get_expression(
            adata,
            use_raw=use_raw,
            layer=layer,
            var_subset=var_subset,
            copy=False,
        ),
    )

    max_components = min(selected_expression_mtx.shape)
    if n_components > max_components:
        raise ValueError(
            f"invalid argument value for 'n_components': "
            f"expected value <= {max_components} but received {n_components!r}"
        )

    resolved_random_state = _as_seed(seed)
    if zero_center:
        if resolved_svd_solver == "auto":
            resolved_svd_solver = "arpack"

        if sparse.issparse(selected_expression_mtx):
            if resolved_svd_solver == "arpack":
                decomposition_mtx = selected_expression_mtx
            else:
                decomposition_mtx = cast(Any, selected_expression_mtx).toarray()
        else:
            decomposition_mtx = np.asarray(selected_expression_mtx)

        estimator = PCA(
            n_components=n_components,
            svd_solver=resolved_svd_solver,
            random_state=resolved_random_state,
        )
    else:
        if resolved_svd_solver == "auto":
            truncated_solver: TruncatedSVDSolver = "randomized"
        else:
            truncated_solver = _as_literal(
                resolved_svd_solver,
                choices=UNCENTERED_PCA_SOLVERS,
                name="svd_solver",
            )

        decomposition_mtx = selected_expression_mtx

        estimator = TruncatedSVD(
            n_components=n_components,
            algorithm=truncated_solver,
            random_state=resolved_random_state,
        )
        resolved_svd_solver = truncated_solver

    with threadpool_limits(limits=n_jobs):
        scores = cast(np.ndarray, estimator.fit_transform(decomposition_mtx))
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
            "var_subset": _format_var_subset(var_subset),
            "layer": layer,
            "use_raw": use_raw,
            "svd_solver": resolved_svd_solver,
            "key_added": key_added,
            "seed": _format_random_state(seed),
            "n_jobs": n_jobs,
        },
        "variance": cast(np.ndarray, estimator.explained_variance_),
        "variance_ratio": cast(np.ndarray, estimator.explained_variance_ratio_),
    }

    return adata if copy else None


def _neighbors_graph(
    adata: AnnData,
    neighbors_key: str,
) -> Tuple[str, Dict[str, Any], Any]:

    if neighbors_key not in adata.uns:
        raise KeyError(
            f"key {neighbors_key!r} not found in adata.uns: "
            f"please run `bt.sct.tl.neighbors(..., key_added={neighbors_key!r})`"
        )

    neighbors = cast(Dict[str, Any], adata.uns[neighbors_key])
    if "connectivities_key" not in neighbors:
        raise KeyError(
            f"key 'connectivities_key' not found in " f"adata.uns[{neighbors_key!r}]"
        )

    connectivities_key = _as_string(
        neighbors["connectivities_key"],
        "connectivities_key",
    )
    if connectivities_key not in adata.obsp:
        raise KeyError(
            f"key {connectivities_key!r} not found in adata.obsp: "
            f"expected connectivities from adata.uns[{neighbors_key!r}]"
        )

    connectivities = cast(Any, adata.obsp[connectivities_key])
    if sparse.issparse(connectivities):
        graph = connectivities.tocoo(copy=True)
    else:
        graph = sparse.coo_matrix(connectivities)

    neighbor_params = cast(Dict[str, Any], neighbors.get("params", {}))
    return connectivities_key, neighbor_params, graph


@anndata_checker
def spectral(
    adata: AnnData,
    n_components: int = 2,
    neighbors_key: str = "neighbors",
    eigen_solver: Optional[EigenSolver] = None,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = None,
    n_jobs: int = 1,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a spectral embedding from a precomputed neighborhood graph.

    Spectral embedding is a nonlinear dimensionality reduction method based on
    the eigenvectors of a graph Laplacian constructed from local neighborhood
    relationships.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    neighbors_key: str (default: 'neighbors')
        Key in `adata.uns` describing the precomputed neighborhood graph.
        Spectral embedding reads the graph from
        `adata.obsp[adata.uns[neighbors_key]["connectivities_key"]]`.
    eigen_solver: {'arpack', 'lobpcg', 'amg'}, optional
        Eigen solver passed to `sklearn.manifold.SpectralEmbedding`.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_se`.
    seed: int, np.random.RandomState, np.random or None, optional
        Random seed or random state used by the embedding estimator.
    n_jobs: int (default: 1)
        Number of threads allocated to the spectral embedding computation.
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
    from threadpoolctl import threadpool_limits

    n_components = _as_positive_integer(n_components, "n_components")
    n_jobs = _as_positive_integer(n_jobs, "n_jobs")

    resolved_random_state = _as_seed(seed)

    eigen_solver = _as_literal(
        eigen_solver,
        choices=EIGEN_SOLVERS,
        name="eigen_solver",
        allow_none=True,
    )

    neighbors_key = _as_string(neighbors_key, "neighbors_key")

    if key_added is None:
        key_added = "X_se"
    else:
        key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata
    connectivities_key, neighbor_params, graph = _neighbors_graph(
        adata,
        neighbors_key,
    )

    with threadpool_limits(limits=n_jobs):
        embedding = cast(
            np.ndarray,
            SpectralEmbedding(
                n_components=n_components,
                affinity="precomputed",
                random_state=resolved_random_state,
                eigen_solver=eigen_solver,
                n_jobs=n_jobs,
            ).fit_transform(graph),
        )

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "spectral",
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivities_key,
        "n_components": n_components,
        "n_neighbors": neighbor_params.get("n_neighbors"),
        "metric": neighbor_params.get("metric"),
        "seed": _format_random_state(seed),
        "eigen_solver": eigen_solver,
        "n_jobs": n_jobs,
    }

    return adata if copy else None


@anndata_checker
def umap(
    adata: AnnData,
    n_components: int = 2,
    neighbors_key: str = "neighbors",
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_iter: int = 500,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: Union[str, np.ndarray] = "spectral",
    a: Optional[float] = None,
    b: Optional[float] = None,
    key_added: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
    copy: bool = False,
) -> Union[AnnData, None]:
    """
    Compute a UMAP embedding from a precomputed neighborhood graph.

    UMAP is a nonlinear dimensionality reduction method based on neighborhood
    graphs and manifold learning.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    neighbors_key: str (default: 'neighbors')
        Key in `adata.uns` describing the precomputed neighborhood graph.
        UMAP reads the graph from
        `adata.obsp[adata.uns[neighbors_key]["connectivities_key"]]`.
    min_dist: float (default: 0.5)
        Effective minimum distance between embedded points.
    spread: float (default: 1.0)
        Effective scale of embedded points.
    n_iter: int (default: 500)
        Number of optimization iterations. Low values may produce poorly
        converged embeddings.
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
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_umap`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by UMAP.
    n_jobs: int (default: 1)
        Number of allocated processors.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with UMAP results added.
        Otherwise, updates `adata` in place and returns None.

        UMAP results are stored in:

        - `adata.obsm[key_added]`: UMAP coordinates;
        - `adata.uns[key_added]`: UMAP metadata.

    References
    ----------
    McInnes et al. (2018). UMAP: Uniform Manifold Approximation and Projection
    for Dimension Reduction. arXiv:1802.03426.
    """

    n_components = _as_positive_integer(n_components, "n_components")
    if not isinstance(n_jobs, int):
        raise TypeError(
            f"unsupported argument type for 'n_jobs': "
            f"expected {int} but received {type(n_jobs)}"
        )

    resolved_random_state = _as_seed(seed)

    neighbors_key = _as_string(neighbors_key, "neighbors_key")

    n_iter = _as_positive_integer(n_iter, "n_iter")

    negative_sample_rate = _as_positive_integer(
        negative_sample_rate,
        "negative_sample_rate",
    )

    if key_added is None:
        key_added = "X_umap"
    else:
        key_added = _as_string(key_added, "key_added")

    adata = adata.copy() if copy else adata

    connectivities_key, neighbor_params, graph = _neighbors_graph(
        adata,
        neighbors_key,
    )
    graph_metric = cast(str, neighbor_params.get("metric", "euclidean"))
    graph_metric_kwargs = cast(Dict[str, Any], neighbor_params.get("metric_kwds", {}))
    representation = cast(
        Optional[str],
        neighbor_params.get("representation", neighbor_params.get("use_rep")),
    )
    n_pcs = cast(Optional[int], neighbor_params.get("n_pcs"))
    representation_mtx = get_representation(
        adata,
        use_rep=representation,
        n_components=n_pcs,
    )

    try:
        umap_umap_module = importlib.import_module("umap.umap_")
    except ImportError as error:
        raise ImportError(
            "umap-learn is required for `bt.sct.tl.umap`. "
            "Install bonesistools with the sctools extra or install umap-learn."
        ) from error

    find_ab_params = cast(Any, getattr(umap_umap_module, "find_ab_params"))
    simplicial_set_embedding = cast(
        Any,
        getattr(umap_umap_module, "simplicial_set_embedding"),
    )
    if a is None or b is None:
        a, b = cast(Tuple[float, float], find_ab_params(spread, min_dist))

    if isinstance(init_pos, str) and init_pos in adata.obsm:
        init = np.asarray(adata.obsm[init_pos], dtype=np.float32)
    else:
        init = init_pos

    embedding, _ = cast(
        Tuple[np.ndarray, Any],
        simplicial_set_embedding(
            data=representation_mtx,
            graph=graph,
            n_components=n_components,
            initial_alpha=alpha,
            a=a,
            b=b,
            gamma=gamma,
            negative_sample_rate=negative_sample_rate,
            n_epochs=n_iter,
            init=init,
            random_state=resolved_random_state,
            metric=graph_metric,
            metric_kwds=graph_metric_kwargs,
            densmap=False,
            densmap_kwds={},
            output_dens=False,
            verbose=False,
        ),
    )

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "umap",
        "neighbors_key": neighbors_key,
        "connectivities_key": connectivities_key,
        "n_components": n_components,
        "seed": _format_random_state(seed),
        "n_neighbors": neighbor_params.get("n_neighbors"),
        "metric": graph_metric,
        "min_dist": min_dist,
        "spread": spread,
        "n_iter": n_iter,
        "alpha": alpha,
        "gamma": gamma,
        "negative_sample_rate": negative_sample_rate,
        "init_pos": init_pos,
        "a": a,
        "b": b,
        "n_jobs": n_jobs,
    }

    return adata if copy else None


@anndata_checker
def tsne(
    adata: AnnData,
    representation: Optional[str] = "X_pca",
    n_pcs: Optional[int] = None,
    n_components: int = 2,
    n_iter: int = 1000,
    perplexity: float = 30.0,
    learning_rate: float = 1000.0,
    early_exaggeration: float = 12.0,
    metric: Metric = "euclidean",
    key_added: Optional[str] = None,
    seed: RandomStateSeed = 0,
    n_jobs: int = 1,
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
    representation: str, optional (default: 'X_pca')
        Representation key in `adata.obsm` used as input.
    n_pcs: int, optional
        Number of input representation dimensions to use. If None, use all
        dimensions in `representation`.
    n_components: int (default: 2)
        Number of embedding dimensions to compute.
    n_iter: int (default: 1000)
        Number of gradient-descent optimization iterations. Low values may
        produce poorly converged embeddings.
    perplexity: float (default: 30.0)
        t-SNE perplexity.
    learning_rate: float (default: 1000.0)
        t-SNE learning rate.
    early_exaggeration: float (default: 12.0)
        t-SNE early exaggeration.
    metric: Metric (default: 'euclidean')
        Distance metric used to compare observations.
    key_added: str, optional
        Key used to store the embedding in `adata.obsm`. Defaults to `X_tsne`.
    seed: int, np.random.RandomState, np.random or None (default: 0)
        Random seed or random state used by t-SNE.
    n_jobs: int (default: 1)
        Number of allocated processors.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with t-SNE results added.
        Otherwise, updates `adata` in place and returns None.

        t-SNE results are stored in:

        - `adata.obsm[key_added]`: t-SNE coordinates;
        - `adata.uns[key_added]`: t-SNE metadata.

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
    n_iter = _as_positive_integer(n_iter, "n_iter")

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
    representation_mtx = get_representation(
        adata,
        use_rep=representation,
        n_components=n_pcs,
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
        tsne_kwargs["max_iter"] = n_iter
    else:
        tsne_kwargs["n_iter"] = n_iter

    embedding = cast(np.ndarray, TSNE(**tsne_kwargs).fit_transform(representation_mtx))

    adata.obsm[key_added] = embedding
    adata.uns[key_added] = {
        "method": "tsne",
        "representation": representation,
        "n_pcs": n_pcs,
        "n_components": n_components,
        "seed": _format_random_state(seed),
        "perplexity": perplexity,
        "metric": metric,
        "early_exaggeration": early_exaggeration,
        "learning_rate": learning_rate,
        "n_iter": n_iter,
        "n_jobs": n_jobs,
    }

    return adata if copy else None
