#!/usr/bin/env python

"""
Plotting utilities for single-cell annotated data.

The `pl` namespace provides embedding, density, composition, distribution and
graph plotting helpers, along with colour palettes and Matplotlib styling
utilities.
"""

from typing import List as _List

from ._barplot import composition
from ._colors import (
    CLASSIC_COLORS,
    EARTH_COLORS,
    HEX_COLORS,
    LIGHT_COLORS,
    PALETTES,
    QUALITATIVE_COLORS,
    RGB_COLORS,
    Palette,
    classic_cm,
    color_cycle,
    earth_cm,
    generate_colormap,
    get_color,
    get_colormap,
    get_palette,
    light_cm,
    rgb,
    rgb2hex,
    rgba,
)
from ._density import cdf, density, ecdf_plot, kde_plot
from ._distribution import boxplot, distribution
from ._figure import set_default_axis, set_default_params
from ._graph import add_graph, draw_paga, graph_overlay, paga, trajectory
from ._scatterplot import embedding, embedding_plot

__all__ = [
    "embedding",
    "embedding_plot",
    "distribution",
    "boxplot",
    "composition",
    "density",
    "cdf",
    "trajectory",
    "graph_overlay",
    "add_graph",
    "paga",
    "draw_paga",
    "kde_plot",
    "ecdf_plot",
    "rgb",
    "rgba",
    "rgb2hex",
    "get_color",
    "get_palette",
    "get_colormap",
    "generate_colormap",
    "CLASSIC_COLORS",
    "LIGHT_COLORS",
    "EARTH_COLORS",
    "QUALITATIVE_COLORS",
    "HEX_COLORS",
    "RGB_COLORS",
    "PALETTES",
    "Palette",
    "color_cycle",
    "classic_cm",
    "light_cm",
    "earth_cm",
    "set_default_params",
    "set_default_axis",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
