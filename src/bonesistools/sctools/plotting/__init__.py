#!/usr/bin/env python

"""
Plotting utilities for single-cell annotated data.

The `pl` namespace provides embedding, density, composition, boxplot and graph
plotting helpers, along with colour palettes and Matplotlib styling utilities.
"""

from typing import List as _List

from ._barplot import composition
from ._boxplot import boxplot
from ._colors import (
    COLORS,
    LIGHT_COLORS,
    QUALITATIVE_COLORS,
    bonesis_cm,
    color_cycle,
    generate_colormap,
    get_color,
    rgb,
    rgb2hex,
)
from ._density import cdf, density, ecdf_plot, kde_plot
from ._figure import set_default_axis, set_default_params
from ._graphplot import draw_paga
from ._scatterplot import embedding, embedding_plot

__all__ = [
    "embedding",
    "embedding_plot",
    "boxplot",
    "composition",
    "density",
    "cdf",
    "draw_paga",
    "kde_plot",
    "ecdf_plot",
    "rgb",
    "rgb2hex",
    "get_color",
    "generate_colormap",
    "COLORS",
    "LIGHT_COLORS",
    "QUALITATIVE_COLORS",
    "color_cycle",
    "bonesis_cm",
    "set_default_params",
    "set_default_axis",
]


def __dir__() -> _List[str]:
    return sorted(set(globals()) | set(__all__))
