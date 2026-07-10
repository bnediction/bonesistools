#!/usr/bin/env python

"""
Plotting utilities for single-cell annotated data.

The `pl` namespace provides embedding, density, composition, distribution and
graph plotting helpers, along with color accessors and Matplotlib styling
utilities.
"""

from typing import List as _List

from ._barplot import composition
from ._colors import (
    generate_colormap,
    get_color,
    get_colormap,
    get_palette,
    rgb,
    rgb2hex,
    rgba,
)
from ._density import cdf, density
from ._density import ecdf_plot as ecdf_plot
from ._density import kde_plot as kde_plot
from ._distribution import boxplot as boxplot
from ._distribution import distribution
from ._figure import set_default_axis, set_default_params
from ._graph import add_graph as add_graph
from ._graph import draw_paga as draw_paga
from ._graph import graph_overlay, paga, trajectory
from ._scatterplot import embedding
from ._scatterplot import embedding_plot as embedding_plot

__all__ = [
    "embedding",
    "distribution",
    "composition",
    "density",
    "cdf",
    "trajectory",
    "graph_overlay",
    "paga",
    "rgb",
    "rgba",
    "rgb2hex",
    "get_color",
    "get_palette",
    "get_colormap",
    "generate_colormap",
    "set_default_params",
    "set_default_axis",
]


def __dir__() -> _List[str]:
    hidden = {
        "add_graph",
        "boxplot",
        "draw_paga",
        "ecdf_plot",
        "embedding_plot",
        "kde_plot",
    }
    return sorted(
        name for name in (set(globals()) | set(__all__)) - hidden if name[0] != "_"
    )
