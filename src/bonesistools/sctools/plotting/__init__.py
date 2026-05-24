#!/usr/bin/env python

"""
Plotting utilities for single-cell annotated data.

The `pl` namespace provides embedding, density, boxplot and graph plotting
helpers, along with colour palettes and Matplotlib styling utilities.
"""

from typing import List

from ._scatterplot import embedding_plot
from ._boxplot import boxplot
from ._graphplot import draw_paga
from ._density import kde_plot, ecdf_plot

from ._colors import (
    rgb,
    rgb2hex,
    get_color,
    generate_colormap,
    COLORS,
    LIGHT_COLORS,
    QUALITATIVE_COLORS,
    color_cycle,
    bonesis_cm,
)

from ._figure import set_default_params, set_default_axis

__all__ = [
    "embedding_plot",
    "boxplot",
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


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
