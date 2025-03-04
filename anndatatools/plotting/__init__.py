#!/usr/bin/env python

from . import _colors as color
from . import _figure as fig
from ._scatterplot import embedding_plot
from ._graphplot import draw_paga
from ._colors import (
    COLORS,
    LIGHT_COLORS,
    QUALITATIVE_COLORS,
    color_cycle,
    bonesis_cm,
    rgb2hex,
    generate_colormap
)
from ._figure import set_default
from ._kde import kde_plot