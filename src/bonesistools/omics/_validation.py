#!/usr/bin/env python

from collections import Counter
from typing import Optional, Union, overload

import numpy as np
import pandas as pd
from anndata import AnnData

from bonesistools._compat import Literal
from bonesistools._warnings import _warn_deprecated

from ._typing import AnnDataAxis, AnnDataAxisWithBoth, VarSubset


@overload
def _as_anndata_axis(
    axis: Union[int, str],
    *,
    allow_both: Literal[False] = False,
    allow_integer: bool = False,
) -> AnnDataAxis: ...


@overload
def _as_anndata_axis(
    axis: Union[int, str],
    *,
    allow_both: Literal[True],
    allow_integer: bool = False,
) -> AnnDataAxisWithBoth: ...


def _as_anndata_axis(
    axis: Union[int, str],
    *,
    allow_both: bool = False,
    allow_integer: bool = False,
) -> AnnDataAxisWithBoth:

    if isinstance(axis, str) and (axis == "obs" or axis == "var"):
        return axis

    if allow_both and isinstance(axis, str) and axis == "both":
        return axis

    if (
        allow_integer
        and isinstance(axis, int)
        and not isinstance(axis, bool)
        and axis in [0, 1]
    ):
        new_axis: AnnDataAxis = "obs" if axis == 0 else "var"
        _warn_deprecated(
            f"`axis={axis}`",
            replacement=f"`axis={new_axis!r}`",
            stacklevel=3,
        )
        return new_axis

    expected = "'obs', 'var' or 'both'" if allow_both else "'obs' or 'var'"
    raise ValueError(
        f"invalid argument value for 'axis': expected {expected} "
        f"but received {axis!r}"
    )


def _as_var_subset(adata: AnnData, var_subset: VarSubset) -> Optional[np.ndarray]:

    if var_subset is None:
        return None

    if isinstance(var_subset, str):
        if var_subset not in adata.var:
            raise KeyError(f"column {var_subset!r} not found in adata.var")

        subset = adata.var[var_subset]
        if not pd.api.types.is_bool_dtype(subset):
            raise TypeError(
                f"unsupported column dtype for 'var_subset': "
                f"expected boolean values in adata.var[{var_subset!r}]"
            )

        mask = np.asarray(subset, dtype=bool)
        if not bool(mask.any()):
            raise ValueError(
                f"invalid argument value for 'var_subset': "
                f"adata.var[{var_subset!r}] selects no variables"
            )
        return mask

    try:
        variables = list(var_subset)
    except TypeError as error:
        raise TypeError(
            f"unsupported argument type for 'var_subset': "
            f"expected {str} or a collection of variable names "
            f"but received {type(var_subset)}"
        ) from error

    invalid_types = [
        variable for variable in variables if not isinstance(variable, str)
    ]
    if invalid_types:
        variable = invalid_types[0]
        raise TypeError(
            f"unsupported element type in 'var_subset': "
            f"expected {str} but received {type(variable)}"
        )

    if not variables:
        raise ValueError(
            "invalid argument value for 'var_subset': "
            "expected at least one variable name"
        )

    requested_counts = Counter(variables)
    missing = [
        variable for variable in requested_counts if variable not in adata.var_names
    ]
    if missing:
        formatted_missing = ", ".join(repr(variable) for variable in missing)
        raise KeyError(f"variable(s) not found in adata.var_names: {formatted_missing}")

    mask = np.asarray(adata.var_names.isin(requested_counts.keys()), dtype=bool)
    if not bool(mask.any()):
        raise ValueError(
            "invalid argument value for 'var_subset': no variables selected"
        )
    return mask
