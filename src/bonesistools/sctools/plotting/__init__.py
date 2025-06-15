#!/usr/bin/env python

from ._scatterplot import embedding_plot
from ._boxplot import boxplot
from ._graphplot import draw_paga
from ._kde import kde_plot

from ._colors import (
    rgb,
    rgb2hex,
    get_color,
    generate_colormap,
    COLORS,
    LIGHT_COLORS,
    QUALITATIVE_COLORS,
    color_cycle,
    bonesis_cm
)

from ._figure import (
    set_default_params,
    set_default_axis
)
