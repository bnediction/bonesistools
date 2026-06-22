#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from typing import (
    Any,
    Callable,
    Dict,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import numpy as np
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure

from ..._warnings import _warn_deprecated_argument
from .._typing import ScData

LegendArgument = Union[bool, Dict[str, Any]]


def set_window_title(fig: Figure, title: str) -> None:
    manager = fig.canvas.manager

    if manager is not None:
        manager.set_window_title(title)


def figure_from_axes(ax: Axes) -> Figure:
    return cast(Figure, ax.figure)


def deprecated_bool_kwarg(
    kwargs: Dict[str, Any],
    old_name: str,
    new_name: str,
    new_value: Any,
    default_value: Any,
    *,
    stacklevel: int = 3,
) -> Any:
    if old_name not in kwargs:
        return new_value

    old_value = kwargs.pop(old_name)
    _warn_deprecated_argument(old_name, new_name, stacklevel=stacklevel)

    if new_value != default_value:
        raise TypeError(
            f"invalid argument combination: use either '{old_name}' "
            f"or '{new_name}', not both"
        )

    return old_value


def _resolve_bool_kwarg(
    kwargs: Dict[str, Any],
    name: str,
    deprecated_name: str,
    default: Optional[bool],
    *,
    stacklevel: int = 3,
) -> Optional[bool]:
    if deprecated_name in kwargs:
        _warn_deprecated_argument(deprecated_name, name, stacklevel=stacklevel)
        if name in kwargs:
            raise TypeError(
                f"invalid argument combination: use either '{deprecated_name}' "
                f"or '{name}', not both"
            )
        kwargs[name] = kwargs.pop(deprecated_name)

    if name not in kwargs:
        return default

    return cast(Optional[bool], kwargs.pop(name))


def rename_deprecated_kwarg(
    kwargs: Dict[str, Any],
    old_name: str,
    new_name: str,
    *,
    stacklevel: int = 3,
) -> None:
    if old_name not in kwargs:
        return None

    if new_name in kwargs:
        raise TypeError(
            f"invalid argument combination: use either '{old_name}' "
            f"or '{new_name}', not both"
        )

    _warn_deprecated_argument(old_name, new_name, stacklevel=stacklevel)
    kwargs[new_name] = kwargs.pop(old_name)
    return None


def _resolve_legend_argument(
    legend: LegendArgument,
    kwargs: Dict[str, Any],
    *,
    stacklevel: int = 3,
) -> Tuple[bool, Dict[str, Any]]:
    for deprecated_name in ("showlegend", "show_legend"):
        if deprecated_name not in kwargs:
            continue

        _warn_deprecated_argument(deprecated_name, "legend", stacklevel=stacklevel)
        if legend is not True:
            raise TypeError(
                f"invalid argument combination: use either '{deprecated_name}' "
                "or 'legend', not both"
            )
        legend = bool(kwargs.pop(deprecated_name))

    for deprecated_name in ("lgd_params", "legend_params"):
        if deprecated_name not in kwargs:
            continue

        _warn_deprecated_argument(deprecated_name, "legend", stacklevel=stacklevel)
        if legend is not True:
            raise TypeError(
                f"invalid argument combination: use either '{deprecated_name}' "
                "or 'legend', not both"
            )
        legend = cast(Dict[str, Any], kwargs.pop(deprecated_name))

    if isinstance(legend, bool):
        return legend, {}
    if isinstance(legend, dict):
        return True, dict(legend)

    raise TypeError(
        f"unsupported argument type for 'legend': "
        f"expected {bool} or {dict} but received {type(legend)}"
    )


def _resolve_legend_options(
    show_legend: bool,
    legend: Optional[Dict[str, Any]],
    kwargs: Dict[str, Any],
    *,
    deprecated_show_name: str = "showlegend",
    stacklevel: int = 3,
) -> Tuple[bool, Optional[Dict[str, Any]]]:
    show_legend = cast(
        bool,
        _resolve_bool_kwarg(
            kwargs,
            "show_legend",
            deprecated_show_name,
            show_legend,
            stacklevel=stacklevel,
        ),
    )

    for deprecated_name in ("lgd_params", "legend_params"):
        if deprecated_name not in kwargs:
            continue

        _warn_deprecated_argument(deprecated_name, "legend", stacklevel=stacklevel)
        if legend is not None:
            raise TypeError(
                f"invalid argument combination: use either '{deprecated_name}' "
                "or 'legend', not both"
            )
        legend = cast(Dict[str, Any], kwargs.pop(deprecated_name))

    if "legend" in kwargs:
        if legend is not None:
            raise TypeError(
                "invalid argument combination: use either deprecated "
                "'lgd_params'/'legend_params' or 'legend', not both"
            )
        legend = cast(Dict[str, Any], kwargs.pop("legend"))

    return show_legend, legend


def colormap_colors(colors: object) -> Sequence[object]:
    if isinstance(colors, ListedColormap):
        return cast(Sequence[object], colors.colors)

    return cast(Sequence[object], colors)


def normalize_color(color: object) -> object:
    if isinstance(color, str):
        return color

    normalized_color = (
        color.tolist() if isinstance(color, np.ndarray) else list(cast(Any, color))
    )
    if len(normalized_color) == 3:
        normalized_color.append(1)

    return normalized_color


def qualitative_color_values(
    count: int,
    palette: Sequence[object],
    generate_colors: Callable[..., Union[Sequence[object], Colormap]],
) -> Union[Sequence[object], Colormap]:
    if len(palette) >= count:
        return palette[0:count]

    return generate_colors(color_number=count)


def colors_from_uns(
    scdata: ScData,
    key: str,
    values: Sequence[object],
) -> Optional[Mapping[object, object]]:
    for uns_key in (f"{key}_color", f"{key}_colors"):
        if uns_key not in scdata.uns:
            continue

        colors = scdata.uns[uns_key]

        if isinstance(colors, MappingABC):
            return cast(Mapping[object, object], colors)

        return {
            value: cast(Sequence[Any], colors)[index]
            for index, value in enumerate(values)
            if index < len(colors)
        }

    return None
