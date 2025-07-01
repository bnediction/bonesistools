#!/usr/bin/env python

from typing import Optional
from ._typing import RGB

import math
import numpy as np

from itertools import cycle
from matplotlib.colors import Colormap, ListedColormap

def rgb(color):
    return list(map(lambda x: x/255, color))

def rgb2hex(rgb: RGB):
    """
    Convert color from rgb format to hexadecimal format.

    Parameters
    ----------
    color: RGB
        Color in RGB format.
    
    Returns
    -------
    Return the color in hexadecimal format.
    """

    r, g, b = rgb
    if isinstance(r, float) and isinstance(g, float) and isinstance(b, float):
        r = round(r*255) if 0<r<1 else int(r)
        g = round(g*255) if 0<g<1 else int(g)
        b = round(b*255) if 0<b<1 else int(b)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

black       = rgb([  0,   0,   0])
white       = rgb([255, 255, 255])
blue        = rgb([  0,  20, 255])
red         = rgb([255,  80,  50])
green       = rgb([ 20, 200,  80])
violet      = rgb([255,  51, 255])
lightgreen  = rgb([ 20, 250,  80])
coral       = rgb([255, 127,  80])
yellow      = rgb([255, 255,   0])
darkred     = rgb([139,   0,   0])
darkyellow  = rgb([204, 204,   0])
lightyellow = rgb([128, 128,   0])
darkorange  = rgb([255, 105,   0])
lightorange = rgb([255, 165,  90])
limegreen   = rgb([ 50, 255,  50])
pink        = rgb([255, 182, 193])
orchid      = rgb([218, 112, 214])
magenta     = rgb([255,   0, 255])
purple      = rgb([128,   0, 128])
indigo      = rgb([ 75,   0, 130])
slateblue   = rgb([ 71,  60, 139])
lightgray   = rgb([211, 211, 211])
gray        = rgb([112, 128, 144])
darkgreen   = rgb([  0, 100,   0])
gold        = rgb([238, 201,   0])
orange      = rgb([255, 128,   0])
salmon      = rgb([198, 113, 113])
maroon      = rgb([128,   0,   0])
beet        = rgb([142,  56, 142])
teal        = rgb([ 56, 142, 142])
olive       = rgb([142, 142,  56])
navy        = rgb([  0,   0, 128])
skyblue     = rgb([135, 206, 235])
beige       = rgb([255, 255, 204])
burgundy    = rgb([128,   0,  32])

COLORS = [
    blue,
    red,
    green,
    orange,
    purple,
    skyblue,
    teal,
    pink,
    magenta,
    darkgreen,
    darkorange,
    darkred,
    maroon,
    olive,
    orchid,
    beet,
    indigo,
    gold,
    navy,
    salmon,
    black,
    lightgreen,
    coral,
    yellow,
    limegreen,
    slateblue,
    darkyellow,
    darkorange,
    lightyellow,
    lightorange,
    burgundy
]

LIGHT_COLORS = [
    skyblue,
    red,
    green,
    orange,
    orchid,
    teal,
    magenta,
    violet,
    olive,
    beet,
    indigo,
    gold,
    navy,
    salmon
]

QUALITATIVE_COLORS = [
    blue,
    red,
    green,
    orange,
    purple,
    indigo,
    pink,
    darkred,
    darkgreen,
    gold,
    maroon
]

color_cycle = cycle(COLORS)

bonesis_cm = ListedColormap(
    colors=COLORS,
    name="bonesis"
)

def get_color(color:str):

    if color in [
        "black", "white", "blue", "red", "green",
        "violet", "lightgreen", "coral", "yellow", "darkyellow",
        "lightyellow", "darkorange", "darkred", "lightorange", "limegreen", "pink",
        "orchid", "magenta", "purple", "indigo", "slateblue",
        "lightgray", "gray", "darkgreen", "gold", "orange",
        "salmon", "maroon", "beet", "teal", "olive",
        "navy", "skyblue", "beige", "burgundy"
    ]:
        return eval(color)
    else:
        raise ValueError(f"color not found: {color}")

def generate_colormap(
    color_number: int = 80,
    shade_number: Optional[int] = None,
    cm: Colormap = bonesis_cm
) -> ListedColormap:
    """
    Create a colormap from another colormap by adding some new colors.

    Parameters
    ----------
    color_number: int (default: 80)
        Number of colors in the returned colormap.
    shade_number: int (optional, default: None)
        Number of shades in the returned colormap.
    cm: matplotlib.colors.Colormap (default: bonesis_cm)
        Initial colormap to use for creating new colormap.

    Returns
    -------
    Return a ListedColormap.
    """

    if color_number <= 0:
        raise ValueError(f"invalid argument value for 'color_number': expected non-null positive value but received '{color_number}'")
    elif color_number <= cm.N:
        return ListedColormap(cm.colors[0:color_number])
    
    if not isinstance(cm, Colormap):
        raise TypeError(f"unsupported argument type for 'cm': expected {Colormap} but received {type(cm)}")

    if shade_number is None:
        shade_number = cm.N
    elif shade_number <= 0:
        raise ValueError(f"invalid argument value for 'shade_number': expected non-null positive value but received '{shade_number}'")
    
    color_number_with_multiply_of_shades = int(math.ceil(color_number / shade_number) * shade_number)
    linearly_uniform_floats = np.arange(color_number_with_multiply_of_shades) / color_number_with_multiply_of_shades
    reorganised_array = linearly_uniform_floats.reshape(shade_number, color_number_with_multiply_of_shades // shade_number).transpose()
    partition_number = reorganised_array.shape[0]

    flatten_reorganised_array = reorganised_array.reshape(-1)

    initial_cm = cm(flatten_reorganised_array)

    lower_partitions_half = partition_number // 2
    upper_partitions_half = partition_number - lower_partitions_half

    lower_half = lower_partitions_half * shade_number
    for i in range(3):
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8/lower_half)

    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = np.ones(shade_number) - initial_cm[lower_half + j * shade_number: lower_half + (j + 1) * shade_number, i]
            modifier = j * modifier / upper_partitions_half
            initial_cm[lower_half + j * shade_number: lower_half + (j + 1) * shade_number, i] += modifier

    return ListedColormap(initial_cm)
