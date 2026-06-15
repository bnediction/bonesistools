#!/usr/bin/env python

from __future__ import annotations

import warnings
from collections import OrderedDict
from copy import deepcopy
from typing import Any, Dict, List, Optional, Union, cast, overload

import numpy as np
from anndata import AnnData
from pandas import DataFrame, Index
from scipy.sparse import csc_matrix, issparse

from ..._compat import Literal
from ..._warnings import _warn_deprecated
from .._typing import anndata_checker


@overload
def merge_duplicate_vars(
    adata: AnnData,
    keep: Literal["first", "consensus"] = "first",
    varm: Literal["nan", "first", "mean"] = "nan",
    varp: Literal["error", "drop"] = "error",
    copy: Literal[True] = True,
) -> AnnData: ...


@overload
def merge_duplicate_vars(
    adata: AnnData,
    keep: Literal["first", "consensus"] = "first",
    varm: Literal["nan", "first", "mean"] = "nan",
    varp: Literal["error", "drop"] = "error",
    *,
    copy: Literal[False],
) -> None: ...


@overload
def merge_duplicate_vars(
    adata: AnnData,
    keep: Literal["first", "consensus"] = "first",
    varm: Literal["nan", "first", "mean"] = "nan",
    varp: Literal["error", "drop"] = "error",
    copy: bool = True,
) -> Union[AnnData, None]: ...


@anndata_checker
def merge_duplicate_vars(
    adata: AnnData,
    keep: Literal["first", "consensus"] = "first",
    varm: Literal["nan", "first", "mean"] = "nan",
    varp: Literal["error", "drop"] = "error",
    copy: bool = True,
) -> Union[AnnData, None]:
    """
    Merge variables with duplicated `adata.var_names`.

    The resulting AnnData object contains one variable per variable name.
    Values in `.X` and in every `.layers` entry are summed across duplicate
    columns. The `.var` table is merged according to `keep`: with
    `keep="first"`, the first candidate row is kept; with `keep="consensus"`,
    each metadata column keeps its value only when all candidate rows agree,
    otherwise the merged value is set to NaN.

    Variable-level matrices stored in `.varm` do not have a universally correct
    aggregation rule, so their behavior is controlled by `varm`. Pairwise
    variable matrices stored in `.varp` are not merged silently; by default, the
    function raises an error when `.varp` is non-empty.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix containing duplicated variable names.
    keep: {'first', 'consensus'} (default: 'first')
        Strategy used to merge `.var` rows. With `'first'`, keep the first
        candidate row. With `'consensus'`, keep a value only if all candidate
        rows have the same value in that column; otherwise set it to NaN.
    varm: {'nan', 'first', 'mean'} (default: 'nan')
        Strategy used to reduce `.varm` arrays. With `'nan'`, rows for merged
        variables are filled with NaN and rows for non-duplicated variables are
        preserved. With `'first'`, the row corresponding to the kept `.var`
        metadata row is preserved when `keep='first'`; when `keep='consensus'`,
        the first candidate row is preserved. With `'mean'`, candidate rows are
        averaged.
    varp: {'error', 'drop'} (default: 'error')
        Strategy used when `.varp` is non-empty. With `'error'`, raise
        `NotImplementedError`. With `'drop'`, discard `.varp`.
    copy: bool (default: True)
        Return a copy instead of modifying `adata`.

    Examples
    --------
    >>> adata.var
         symbol         biotype
    Aars   Aars  protein_coding
    Aars  Aars1  protein_coding
    >>> merged = bt.sct.pp.merge_duplicate_vars(adata, keep="consensus")
    >>> merged.var
         symbol         biotype
    Aars    NaN  protein_coding

    Returns
    -------
    AnnData or None
        If `copy=True`, returns a copy of `adata` with duplicated variables
        merged. Otherwise, updates `adata` in place and returns None.

        Merged variable data are stored in:

        - `adata.var_names`: deduplicated variable names;
        - `adata.var`: merged variable annotations;
        - `adata.X`: summed expression matrix;
        - `adata.layers`: summed layer matrices;
        - `adata.varm`: merged variable mappings according to `varm`;
        - `adata.varp`: pairwise variable matrices handled according to `varp`.
    """

    if keep not in ["first", "consensus"]:
        raise ValueError(
            "invalid argument value for 'keep': expected 'first' or "
            f"'consensus' but received {keep!r}"
        )
    if varm not in ["nan", "first", "mean"]:
        raise ValueError(
            "invalid argument value for 'varm': expected 'nan', 'first' or "
            f"'mean' but received {varm!r}"
        )
    if varp not in ["error", "drop"]:
        raise ValueError(
            "invalid argument value for 'varp': expected 'error' or 'drop' "
            f"but received {varp!r}"
        )

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        source = adata.copy() if copy else adata

    groups = _group_var_positions(source)
    if all(len(indices) == 1 for indices in groups.values()):
        return source if copy else None

    if len(source.varp) > 0 and varp == "error":
        raise NotImplementedError(
            "`merge_duplicate_vars` does not merge `.varp`; pass "
            "`varp='drop'` to discard pairwise variable matrices."
        )

    selected_positions = [positions[0] for positions in groups.values()]
    group_matrix = _var_group_matrix(source.n_vars, groups)

    merged = AnnData(
        X=_sum_var_groups(cast(Any, source.X), group_matrix),
        obs=cast(DataFrame, source.obs).copy(),
        var=_merge_var(source, groups, selected_positions, keep=keep),
        uns=deepcopy(source.uns),
        obsm=_copy_axis_mapping(source.obsm),
        obsp=_copy_axis_mapping(source.obsp),
        layers={
            key: _sum_var_groups(cast(Any, value), group_matrix)
            for key, value in source.layers.items()
        },
        varm={
            key: _merge_varm_value(
                key,
                value,
                groups,
                selected_positions,
                strategy=varm,
            )
            for key, value in source.varm.items()
        },
    )

    if copy:
        return merged

    adata._init_as_actual(merged)
    return None


@overload
def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None,
    copy: Literal[True] = True,
) -> AnnData: ...


@overload
def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None,
    *,
    copy: Literal[False],
) -> None: ...


@overload
def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None,
    copy: bool = True,
) -> Union[AnnData, None]: ...


def var_names_merge_duplicates(
    adata: AnnData,
    var_names_column: Optional[str] = None,
    copy: bool = True,
) -> Union[AnnData, None]:
    """
    Deprecated. Merge duplicated variable names.

    Use `merge_duplicate_vars` instead.
    """

    _warn_deprecated(
        "`bt.sct.pp.var_names_merge_duplicates`",
        replacement="`bt.sct.pp.merge_duplicate_vars`",
        stacklevel=2,
    )
    if var_names_column is not None:
        _warn_deprecated(
            "`var_names_column`",
            replacement="`keep` to control `.var` merging",
            stacklevel=2,
        )
    return merge_duplicate_vars(
        adata=adata,
        copy=copy,
    )


def _group_var_positions(adata: AnnData) -> "OrderedDict[str, List[int]]":
    groups = OrderedDict()
    for position, var_name in enumerate(adata.var_names):
        groups.setdefault(str(var_name), []).append(position)
    return groups


def _var_group_matrix(
    n_vars: int,
    groups: "OrderedDict[str, List[int]]",
) -> csc_matrix:
    row_positions = list()
    column_positions = list()

    for column_position, positions in enumerate(groups.values()):
        row_positions.extend(positions)
        column_positions.extend([column_position] * len(positions))

    return csc_matrix(
        (
            np.ones(len(row_positions), dtype=np.int8),
            (row_positions, column_positions),
        ),
        shape=(n_vars, len(groups)),
    )


def _sum_var_groups(matrix: Any, group_matrix: csc_matrix) -> Any:
    if issparse(matrix):
        return (matrix @ group_matrix).tocsr()
    return np.asarray(matrix) @ group_matrix.toarray()


def _merge_var(
    adata: AnnData,
    groups: "OrderedDict[str, List[int]]",
    selected_positions: List[int],
    keep: Literal["first", "consensus"],
) -> DataFrame:
    source_var = cast(DataFrame, adata.var)
    merged_index = Index(list(groups.keys()))
    merged_index.name = (
        source_var.index.name if isinstance(source_var.index.name, str) else None
    )

    if keep == "first":
        merged_var = source_var.iloc[selected_positions].copy()
        merged_var.index = merged_index
        return merged_var

    merged_var = DataFrame(
        [
            {
                column: _consensus_value(source_var.iloc[positions][column])
                for column in source_var.columns
            }
            for positions in groups.values()
        ],
        index=merged_index,
        columns=source_var.columns,
    )
    merged_var.index = merged_index
    return merged_var


def _consensus_value(values: Any) -> Any:
    first = values.iloc[0]
    if all(_metadata_values_equal(first, value) for value in values.iloc[1:]):
        return first
    return np.nan


def _metadata_values_equal(left: Any, right: Any) -> bool:
    try:
        left_is_nan = bool(np.isscalar(left) and np.isnan(left))
        right_is_nan = bool(np.isscalar(right) and np.isnan(right))
    except TypeError:
        left_is_nan = False
        right_is_nan = False

    if left_is_nan or right_is_nan:
        return left_is_nan and right_is_nan

    return bool(left == right)


def _copy_axis_mapping(mapping: Any) -> Dict[str, Any]:
    return {
        key: value.copy() if hasattr(value, "copy") else deepcopy(value)
        for key, value in mapping.items()
    }


def _merge_varm_value(
    key: str,
    value: Any,
    groups: "OrderedDict[str, List[int]]",
    selected_positions: List[int],
    strategy: Literal["nan", "first", "mean"],
) -> Any:
    if isinstance(value, DataFrame):
        return _merge_varm_dataframe(
            key,
            value,
            groups,
            selected_positions,
            strategy=strategy,
        )
    return _merge_varm_array(
        key,
        np.asarray(value),
        groups,
        selected_positions,
        strategy=strategy,
    )


def _merge_varm_dataframe(
    key: str,
    value: DataFrame,
    groups: "OrderedDict[str, List[int]]",
    selected_positions: List[int],
    strategy: Literal["nan", "first", "mean"],
) -> DataFrame:
    if strategy == "mean":
        try:
            merged_values = [
                value.iloc[positions].astype(float).mean(axis=0).to_numpy()
                for positions in groups.values()
            ]
        except ValueError as error:
            raise TypeError(
                f"varm='mean' requires numeric values for `.varm[{key!r}]`"
            ) from error
        return DataFrame(
            merged_values,
            index=Index(list(groups.keys())),
            columns=value.columns,
        )

    merged = value.iloc[selected_positions].copy()
    merged.index = Index(list(groups.keys()))
    if strategy == "nan":
        for position, source_positions in enumerate(groups.values()):
            if len(source_positions) > 1:
                merged.iloc[position, :] = np.nan
    return merged


def _merge_varm_array(
    key: str,
    value: np.ndarray,
    groups: "OrderedDict[str, List[int]]",
    selected_positions: List[int],
    strategy: Literal["nan", "first", "mean"],
) -> np.ndarray:
    if strategy == "mean":
        if not np.issubdtype(value.dtype, np.number):
            raise TypeError(f"varm='mean' requires numeric values for `.varm[{key!r}]`")
        return np.asarray(
            [value[positions].mean(axis=0) for positions in groups.values()]
        )

    merged = value[selected_positions].copy()
    if strategy == "nan":
        if not (
            np.issubdtype(merged.dtype, np.floating)
            or np.issubdtype(merged.dtype, np.complexfloating)
        ):
            if np.issubdtype(merged.dtype, np.number) or np.issubdtype(
                merged.dtype,
                np.bool_,
            ):
                merged = merged.astype(float)
            else:
                merged = merged.astype(object)
        for position, source_positions in enumerate(groups.values()):
            if len(source_positions) > 1:
                merged[position] = np.nan
    return merged
