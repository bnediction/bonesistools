#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from typing import Any, Callable, Mapping, Optional, Sequence, Union, cast

import numpy as np
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure

from .._typing import ScData


def set_window_title(fig: Figure, title: str) -> None:
    manager = fig.canvas.manager

    if manager is not None:
        manager.set_window_title(title)


def figure_from_axes(ax: Axes) -> Figure:
    return cast(Figure, ax.figure)


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
