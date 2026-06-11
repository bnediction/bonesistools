#!/usr/bin/env python

import warnings
from typing import Any, Optional, cast

from anndata import AnnData
from numpy import ndarray

from .._typing import (
    ScData,
    anndata_checker,
    anndata_or_mudata_checker,
)


@anndata_checker
def choose_matrix_representation(
    adata: AnnData,
    use_raw: bool = False,
    layer: Optional[str] = None,
    copy: bool = True,
) -> ndarray:
    """
    Select the expression matrix from an AnnData object.

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
    copy: bool (default: True)
        Return a copy of the selected matrix.

    Returns
    -------
    ndarray
        Selected matrix.

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
    elif use_raw:
        matrix = adata.raw.X
    else:
        matrix = adata.X

    if copy:
        return cast(ndarray, matrix.copy())
    else:
        return cast(ndarray, matrix)


def choose_mtx_representation(*args: Any, **kwargs: Any) -> ndarray:
    """
    Deprecated alias for `choose_matrix_representation`.
    """

    warnings.warn(
        "`bt.sct.tl.choose_mtx_representation` is deprecated; use "
        "`bt.sct.tl.choose_matrix_representation` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return choose_matrix_representation(*args, **kwargs)


@anndata_or_mudata_checker
def choose_representation(
    scdata: ScData,  # type: ignore
    use_rep: Optional[str] = "X_pca",
    n_components: Optional[int] = None,
) -> ndarray:
    """
    Select and optionally truncate an embedding representation.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    use_rep: str (default: "X_pca")
        Representation key in `scdata.obsm`.
    n_components: int, optional
        Number of dimensions to use. If None, use all dimensions.

    Returns
    -------
    ndarray
        Selected representation.

    Raises
    ------
    KeyError
        If `use_rep` is not found in `scdata.obsm`.
    """

    if use_rep is None:
        use_rep = "X_pca"

    if use_rep not in scdata.obsm:
        if use_rep == "X_pca":
            raise KeyError(
                "key 'X_pca' not found in scdata.obsm: please run scanpy.tl.pca"
            )
        else:
            raise KeyError(f"key '{use_rep}' not found in scdata.obsm")

    if n_components is None:
        return cast(ndarray, scdata.obsm[use_rep])
    else:
        return cast(ndarray, scdata.obsm[use_rep][:, :n_components])


@anndata_or_mudata_checker
def _get_distances(
    scdata: ScData,  # type: ignore
    obsp: Optional[str] = None,
    neighbors_key: Optional[str] = None,
) -> ndarray:
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
    ndarray
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
        return cast(ndarray, scdata.obsp[obsp])
    elif neighbors_key is not None:
        distances_key = scdata.uns[neighbors_key]["distances_key"]
        return cast(ndarray, scdata.obsp[distances_key])
    else:
        if "neighbors" in scdata.uns:
            distances_key = scdata.uns["neighbors"]["distances_key"]
            return cast(ndarray, scdata.obsp[distances_key])
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
) -> ndarray:
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
    ndarray
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
        return cast(ndarray, scdata.obsp[obsp])
    elif neighbors_key is not None:
        connectivities_key = scdata.uns[neighbors_key]["connectivities_key"]
        return cast(ndarray, scdata.obsp[connectivities_key])
    else:
        if "neighbors" in scdata.uns:
            connectivities_key = scdata.uns["neighbors"]["connectivities_key"]
            return cast(ndarray, scdata.obsp[connectivities_key])
        else:
            raise KeyError(
                "connectivities not found in 'scdata': "
                "please run scanpy.pp.neighbors or specify 'obsp' or 'neighbors_key'"
            )
