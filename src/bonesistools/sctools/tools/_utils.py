#!/usr/bin/env python

from collections.abc import Collection as CollectionInstance
from typing import Any, Optional, Tuple, cast, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from pandas import Index
from scipy import sparse

from ..._validation import _as_string
from ..._warnings import _warn_deprecated_argument
from .._typing import (
    AnnDataAxis,
    Matrix,
    ScData,
    VarSubset,
    anndata_checker,
    anndata_or_mudata_checker,
)
from .._validation import _as_anndata_axis, _as_var_subset

_UNSET = object()


def _resolve_representation_argument(
    representation: Any,
    use_rep: Any,
    *,
    default: Optional[str],
    stacklevel: int,
    obsm: Any = _UNSET,
) -> Optional[str]:

    if obsm is not _UNSET:
        _warn_deprecated_argument("obsm", "representation", stacklevel=stacklevel)
        if representation is not _UNSET:
            raise TypeError(
                "received both 'representation' and deprecated 'obsm'; "
                "please use only 'representation'"
            )
        representation = obsm

    if use_rep is not _UNSET:
        _warn_deprecated_argument("use_rep", "representation", stacklevel=stacklevel)
        if representation is not _UNSET:
            raise TypeError(
                "received both 'representation' and deprecated 'use_rep'; "
                "please use only 'representation'"
            )
        representation = use_rep
    elif representation is _UNSET:
        representation = default

    return cast(Optional[str], representation)


@anndata_checker
def get_expression(
    adata: AnnData,
    use_raw: bool = False,
    layer: Optional[str] = None,
    var_subset: VarSubset = None,
    copy: bool = True,
) -> Matrix:
    """
    Get the expression matrix from an AnnData object.

    The matrix is selected from `adata.X`, `adata.raw.X`, or a named layer.
    `use_raw` and `layer` are mutually exclusive.

    Parameters
    ----------
    adata: AnnData
        Unimodal or multimodal annotated data matrix.
    use_raw: bool, optional
        Use `adata.raw.X` instead of `adata.X`.
    layer: str, optional
        Layer to use instead of `adata.X`.
    var_subset: str or collection of str, optional
        Variables to select. If a string is provided, it is interpreted as a
        boolean column in `adata.var`. If a collection is provided, it is
        interpreted as variable names.
    copy: bool (default: True)
        Return a copy of the selected matrix.

    Returns
    -------
    ndarray or sparse matrix
        Selected expression matrix.

    Raises
    ------
    ValueError
        If `use_raw` and `layer` are both specified.
    """

    if use_raw and layer is not None:
        raise ValueError(
            "invalid argument combination: "
            "'use_raw' and 'layer' cannot be both specified"
        )
    elif layer is not None:
        matrix: Any = adata.layers[layer]
        var_names = adata.var_names
        var = adata.var
    elif use_raw:
        matrix = adata.raw.X
        var_names = adata.raw.var_names
        var = adata.raw.var
    else:
        matrix = adata.X
        var_names = adata.var_names
        var = adata.var

    if use_raw and var_subset is not None:
        if isinstance(var_subset, str):
            if var_subset not in var:
                raise KeyError(f"column {var_subset!r} not found in adata.raw.var")
            subset = var[var_subset]
            if not hasattr(subset, "dtype") or subset.dtype != bool:
                raise TypeError(
                    f"unsupported column dtype for 'var_subset': "
                    f"expected boolean values in adata.raw.var[{var_subset!r}]"
                )
            mask = subset.to_numpy(dtype=bool)
        elif isinstance(var_subset, CollectionInstance):
            variables = list(var_subset)
            if not variables:
                raise ValueError(
                    "invalid argument value for 'var_subset': "
                    "expected at least one variable name"
                )
            invalid_types = [
                variable for variable in variables if not isinstance(variable, str)
            ]
            if invalid_types:
                variable = invalid_types[0]
                raise TypeError(
                    f"unsupported element type in 'var_subset': "
                    f"expected {str} but received {type(variable)}"
                )
            missing = [variable for variable in variables if variable not in var_names]
            if missing:
                formatted_missing = ", ".join(repr(variable) for variable in missing)
                raise KeyError(
                    "variable(s) not found in adata.raw.var_names: "
                    f"{formatted_missing}"
                )
            mask = var_names.isin(variables)
        else:
            raise TypeError(
                f"unsupported argument type for 'var_subset': "
                f"expected {str} or a collection of variable names "
                f"but received {type(var_subset)}"
            )
    else:
        mask = _as_var_subset(adata, var_subset)
    if mask is not None:
        matrix = cast(Any, matrix)[:, mask]

    if copy:
        return cast(Matrix, matrix.copy())
    else:
        return cast(Matrix, matrix)


def _get_expression_with_gene_names(
    adata: AnnData,
    expression: Optional[str],
    var_subset: VarSubset,
) -> Tuple[Matrix, Index]:
    """
    Get an expression matrix and the matching gene names.

    This helper mirrors the public `expression` API used by tools that need
    both the selected matrix and gene labels. The returned gene names are
    aligned with the matrix columns after applying `var_subset`.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    expression: str, optional
        Expression source. If None or `"X"`, use `adata.X`. If `"raw.X"`, use
        `adata.raw.X`. Otherwise, interpret as a layer key in `adata.layers`.
    var_subset: str or collection of str, optional
        Variables to select. If a string is provided, it is interpreted as a
        boolean column in the selected `.var` annotation. If a collection is
        provided, it is interpreted as variable names.

    Returns
    -------
    tuple
        Selected expression matrix and matching gene names.
    """

    if expression is None or expression in {"X", ".X"}:
        expression_mtx: Any = adata.X
        gene_names = adata.var_names
        var = adata.var
        var_label = "adata.var"
        gene_names_label = "adata.var_names"
    elif expression in {"raw", "raw.X", ".raw.X"}:
        if adata.raw is None:
            raise ValueError(
                "invalid argument value for 'expression': "
                "adata.raw is required when expression='raw.X'"
            )
        expression_mtx = adata.raw.X
        gene_names = adata.raw.var_names
        var = adata.raw.var
        var_label = "adata.raw.var"
        gene_names_label = "adata.raw.var_names"
    else:
        expression = _as_string(expression, "expression")
        expression_mtx = adata.layers[expression]
        gene_names = adata.var_names
        var = adata.var
        var_label = "adata.var"
        gene_names_label = "adata.var_names"

    if var_subset is None:
        return cast(Matrix, expression_mtx), gene_names

    if isinstance(var_subset, str):
        if var_subset not in var:
            raise KeyError(f"column {var_subset!r} not found in {var_label}")

        subset = var[var_subset]
        if not hasattr(subset, "dtype") or not pd.api.types.is_bool_dtype(subset):
            raise TypeError(
                f"unsupported column dtype for 'var_subset': "
                f"expected boolean values in {var_label}[{var_subset!r}]"
            )
        mask = np.asarray(subset, dtype=bool)
        if not bool(mask.any()):
            raise ValueError(
                f"invalid argument value for 'var_subset': "
                f"{var_label}[{var_subset!r}] selects no variables"
            )
        return cast(Matrix, expression_mtx[:, mask]), gene_names[mask]

    if not isinstance(var_subset, CollectionInstance):
        raise TypeError(
            f"unsupported argument type for 'var_subset': "
            f"expected collection of {str} but received {type(var_subset)}"
        )

    variables = list(var_subset)
    if not variables:
        raise ValueError(
            "invalid argument value for 'var_subset': "
            "expected at least one variable name"
        )
    invalid_types = [
        variable for variable in variables if not isinstance(variable, str)
    ]
    if invalid_types:
        variable = invalid_types[0]
        raise TypeError(
            f"unsupported element type in 'var_subset': "
            f"expected {str} but received {type(variable)}"
        )

    positions = gene_names.get_indexer(Index(variables))
    missing = [
        variable for variable, position in zip(variables, positions) if position < 0
    ]
    if missing:
        formatted_missing = ", ".join(repr(variable) for variable in missing)
        raise KeyError(
            f"variable(s) not found in {gene_names_label}: {formatted_missing}"
        )

    mask = gene_names.isin(variables)
    return cast(Matrix, expression_mtx[:, mask]), gene_names[mask]


def _as_dense_matrix_chunk(expression_mtx: Matrix, start: int, end: int) -> np.ndarray:
    """
    Get a dense column chunk from a dense or sparse matrix.

    Parameters
    ----------
    expression_mtx: ndarray or sparse matrix
        Matrix to slice.
    start: int
        First column index.
    end: int
        Stop column index.

    Returns
    -------
    ndarray
        Dense two-dimensional matrix chunk.
    """

    if getattr(expression_mtx, "ndim", 2) != 2:
        raise ValueError("invalid expression matrix: expected a two-dimensional matrix")

    matrix_chunk = cast(Any, expression_mtx)[:, start:end]
    dense_chunk = (
        matrix_chunk.toarray()
        if sparse.issparse(matrix_chunk)
        else np.asarray(matrix_chunk)
    )
    if dense_chunk.ndim != 2:
        raise ValueError("invalid expression matrix: expected a two-dimensional matrix")
    return dense_chunk


@overload
def get_representation(
    scdata: ScData,
    obsm: Optional[str] = "X_pca",
    n_components: Optional[int] = None,
) -> Matrix: ...


@overload
def get_representation(
    scdata: ScData,
    obsm: Optional[str] = "X_pca",
    n_components: Optional[int] = None,
    *,
    use_rep: Any = _UNSET,
) -> Matrix: ...


@anndata_or_mudata_checker
def get_representation(
    scdata: ScData,  # type: ignore
    obsm: Any = _UNSET,
    n_components: Optional[int] = None,
    *,
    use_rep: Any = _UNSET,
) -> Matrix:
    """
    Get and optionally truncate an observation representation.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obsm: str (default: "X_pca")
        Representation key in `scdata.obsm`.
    use_rep: str, optional
        Deprecated alias for `obsm`.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.

    Returns
    -------
    ndarray or sparse matrix
        Selected representation.

    Raises
    ------
    KeyError
        If `obsm` is not found in `scdata.obsm`.
    """

    obsm = _resolve_obsm_argument(
        obsm,
        use_rep,
        default="X_pca",
        stacklevel=2,
    )
    if obsm is None:
        obsm = "X_pca"

    if obsm not in scdata.obsm:
        if obsm == "X_pca":
            raise KeyError(
                "key 'X_pca' not found in scdata.obsm: "
                "please run bonesistools.sct.tl.pca"
            )
        else:
            raise KeyError(f"key '{obsm}' not found in scdata.obsm")

    if n_components is None:
        return cast(Matrix, scdata.obsm[obsm])
    else:
        return cast(Matrix, scdata.obsm[obsm][:, :n_components])


def _resolve_obsm_argument(
    obsm: Any,
    use_rep: Any,
    *,
    default: Optional[str],
    stacklevel: int,
) -> Optional[str]:

    if use_rep is not _UNSET:
        _warn_deprecated_argument("use_rep", "obsm", stacklevel=stacklevel)
        if obsm is not _UNSET:
            raise TypeError(
                "received both 'obsm' and deprecated 'use_rep'; "
                "please use only 'obsm'"
            )
        obsm = use_rep
    elif obsm is _UNSET:
        obsm = default

    return cast(Optional[str], obsm)


@anndata_checker
def get_pairwise(
    adata: AnnData,
    key: str,
    axis: AnnDataAxis = "obs",
) -> Matrix:
    """
    Get a pairwise matrix from an AnnData object.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    key: str
        Pairwise matrix key.
    axis: {"obs", "var"} (default: "obs")
        Axis whose pairwise matrix is retrieved. If `"obs"`, read from
        `adata.obsp[key]`. If `"var"`, read from `adata.varp[key]`.

    Returns
    -------
    ndarray or sparse matrix
        Selected pairwise matrix.
    """

    axis = _as_anndata_axis(axis, allow_both=False)
    if axis == "obs":
        return cast(Matrix, adata.obsp[key])
    else:
        return cast(Matrix, adata.varp[key])


@anndata_or_mudata_checker
def _get_distances(
    scdata: ScData,  # type: ignore
    obsp: Optional[str] = None,
    neighbors_key: Optional[str] = None,
) -> Matrix:
    """
    Retrieve a precomputed distance matrix.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obsp: str, optional
        Key in `scdata.obsp`.
    neighbors_key: str, optional
        Key in `scdata.uns` storing neighborhood metadata.

    Returns
    -------
    ndarray or sparse matrix
        Distance matrix.

    Raises
    ------
    ValueError
        If `obsp` and `neighbors_key` are both specified.
    KeyError
        If no distances matrix can be found from the provided arguments or
        from `scdata.uns["neighbors"]`.
    """

    if obsp is not None and neighbors_key is not None:
        raise ValueError(
            "invalid argument combination: "
            "'obsp' and 'neighbors_key' cannot be both specified"
        )
    elif obsp is not None:
        return cast(Matrix, scdata.obsp[obsp])
    elif neighbors_key is not None:
        distances_key = scdata.uns[neighbors_key]["distances_key"]
        return cast(Matrix, scdata.obsp[distances_key])
    else:
        if "neighbors" in scdata.uns:
            distances_key = scdata.uns["neighbors"]["distances_key"]
            return cast(Matrix, scdata.obsp[distances_key])
        else:
            raise KeyError(
                "distances not found in 'scdata': "
                "please run scanpy.pp.neighbors or specify 'obsp' or 'neighbors_key'"
            )


@anndata_or_mudata_checker
def _get_connectivities(
    scdata: ScData,  # type: ignore
    obsp: Optional[str] = None,
    neighbors_key: Optional[str] = None,
) -> Matrix:
    """
    Retrieve a precomputed connectivity matrix.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obsp: str, optional
        Key in `scdata.obsp`.
    neighbors_key: str, optional
        Key in `scdata.uns` storing neighborhood metadata.

    Returns
    -------
    ndarray or sparse matrix
        Connectivity matrix.

    Raises
    ------
    ValueError
        If `obsp` and `neighbors_key` are both specified.
    KeyError
        If no connectivities matrix can be found from the provided arguments or
        from `scdata.uns["neighbors"]`.
    """

    if obsp is not None and neighbors_key is not None:
        raise ValueError(
            "invalid argument combination: "
            "'obsp' and 'neighbors_key' cannot be both specified"
        )
    elif obsp is not None:
        return cast(Matrix, scdata.obsp[obsp])
    elif neighbors_key is not None:
        connectivities_key = scdata.uns[neighbors_key]["connectivities_key"]
        return cast(Matrix, scdata.obsp[connectivities_key])
    else:
        if "neighbors" in scdata.uns:
            connectivities_key = scdata.uns["neighbors"]["connectivities_key"]
            return cast(Matrix, scdata.obsp[connectivities_key])
        else:
            raise KeyError(
                "connectivities not found in 'scdata': "
                "please run scanpy.pp.neighbors or specify 'obsp' or 'neighbors_key'"
            )
