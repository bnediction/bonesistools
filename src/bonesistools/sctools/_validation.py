#!/usr/bin/env python

from typing import Union, overload

from bonesistools._compat import Literal
from bonesistools._warnings import _warn_deprecated

from ._typing import AnnDataAxis, AnnDataAxisWithBoth


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
