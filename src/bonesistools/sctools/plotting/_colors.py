#!/usr/bin/env python

from typing import Optional

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore
from ._typing import RGB

import math
import numpy as np

from itertools import cycle
from matplotlib.colors import Colormap, ListedColormap


def rgb(color):
    """
    Convert 0-255 RGB channels to 0-1 RGB channels.

    Parameters
    ----------
    color: sequence of numbers
        RGB channels in the 0-255 range.

    Returns
    -------
    list
        RGB channels scaled to the 0-1 range.
    """

    return list(map(lambda x: x / 255, color))


def rgb2hex(rgb):
    """
    Convert color from RGB format to hexadecimal format.

    Parameters
    ----------
    rgb: RGB
        Color in RGB format.

    Returns
    -------
    str
        Color in hexadecimal format.

    Raises
    ------
    TypeError
        If one of the RGB channels is not numeric.
    """

    r, g, b = rgb

    if all(isinstance(x, (int, float)) for x in (r, g, b)):
        if all(0 <= x <= 1 for x in (r, g, b)):
            r, g, b = [round(x * 255) for x in (r, g, b)]
        else:
            r, g, b = [int(x) for x in (r, g, b)]
    else:
        raise TypeError(
            f"unsupported argument type for 'rgb': "
            f"expected numeric values but received {[type(x) for x in (r, g, b)]}"
        )

    return "#{:02x}{:02x}{:02x}".format(r, g, b)


black = rgb([0, 0, 0])
white = rgb([255, 255, 255])
blue = rgb([0, 20, 255])
red = rgb([255, 80, 50])
green = rgb([20, 200, 80])
violet = rgb([255, 51, 255])
lightgreen = rgb([20, 250, 80])
coral = rgb([255, 127, 80])
yellow = rgb([255, 255, 0])
darkred = rgb([139, 0, 0])
darkyellow = rgb([204, 204, 0])
lightyellow = rgb([128, 128, 0])
darkorange = rgb([255, 105, 0])
lightorange = rgb([255, 165, 90])
limegreen = rgb([50, 255, 50])
pink = rgb([255, 182, 193])
orchid = rgb([218, 112, 214])
magenta = rgb([255, 0, 255])
purple = rgb([128, 0, 128])
indigo = rgb([75, 0, 130])
slateblue = rgb([71, 60, 139])
lightgray = rgb([211, 211, 211])
gray = rgb([112, 128, 144])
darkgreen = rgb([0, 100, 0])
gold = rgb([238, 201, 0])
orange = rgb([255, 128, 0])
salmon = rgb([198, 113, 113])
maroon = rgb([128, 0, 0])
beet = rgb([142, 56, 142])
teal = rgb([56, 142, 142])
olive = rgb([142, 142, 56])
navy = rgb([0, 0, 128])
skyblue = rgb([135, 206, 235])
beige = rgb([255, 255, 204])
burgundy = rgb([128, 0, 32])

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
    burgundy,
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
    salmon,
]

QUALITATIVE_COLORS = [
    blue,
    red,
    green,
    orange,
    purple,
    pink,
    darkred,
    darkgreen,
    gold,
    indigo,
    maroon,
]

color_cycle = cycle(COLORS)

bonesis_cm = ListedColormap(colors=COLORS, name="bonesis")


def get_color(color: str, color_type: Literal["rgb", "hex"] = "rgb"):
    """
    Return a named color in RGB or hexadecimal format.

    Parameters
    ----------
    color: str
        Name of a color defined in this module.
    color_type: {"rgb", "hex"} (default: "rgb")
        Output color format.

    Returns
    -------
    list or str
        Color value in the requested format.

    Raises
    ------
    ValueError
        If `color` is unknown or `color_type` is not `"rgb"` or `"hex"`.
    """

    if color in [
        "black",
        "white",
        "blue",
        "red",
        "green",
        "violet",
        "lightgreen",
        "coral",
        "yellow",
        "darkyellow",
        "lightyellow",
        "darkorange",
        "darkred",
        "lightorange",
        "limegreen",
        "pink",
        "orchid",
        "magenta",
        "purple",
        "indigo",
        "slateblue",
        "lightgray",
        "gray",
        "darkgreen",
        "gold",
        "orange",
        "salmon",
        "maroon",
        "beet",
        "teal",
        "olive",
        "navy",
        "skyblue",
        "beige",
        "burgundy",
    ]:
        if color_type == "rgb":
            return eval(color)
        elif color_type == "hex":
            return rgb2hex(eval(color))
        else:
            raise ValueError(
                f"invalid argument value for 'color_type': "
                f"expected 'rgb' or 'hex' but received {color_type!r}"
            )
    else:
        raise ValueError(f"invalid argument value for 'color': color not found: {color!r}")


def generate_colormap(
    color_number: int = 80,
    shade_number: Optional[int] = None,
    cm: Colormap = bonesis_cm,
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

    Raises
    ------
    TypeError
        If `cm` is not a matplotlib Colormap.
    ValueError
        If `color_number` or `shade_number` is not strictly positive.
    """

    if not isinstance(cm, Colormap):
        raise TypeError(
            f"unsupported argument type for 'cm': "
            f"expected {Colormap} but received {type(cm)}"
        )

    if color_number <= 0:
        raise ValueError(
            f"invalid argument value for 'color_number': "
            f"expected non-null positive value but received {color_number!r}"
        )
    elif color_number <= cm.N:
        return ListedColormap(cm.colors[0:color_number])

    if shade_number is None:
        shade_number = cm.N
    elif shade_number <= 0:
        raise ValueError(
            f"invalid argument value for 'shade_number': "
            f"expected non-null positive value but received {shade_number!r}"
        )

    color_number_with_multiply_of_shades = int(
        math.ceil(color_number / shade_number) * shade_number
    )
    linearly_uniform_floats = (
        np.arange(color_number_with_multiply_of_shades)
        / color_number_with_multiply_of_shades
    )
    reorganised_array = linearly_uniform_floats.reshape(
        shade_number, color_number_with_multiply_of_shades // shade_number
    ).transpose()
    partition_number = reorganised_array.shape[0]

    flatten_reorganised_array = reorganised_array.reshape(-1)

    initial_cm = cm(flatten_reorganised_array)

    lower_partitions_half = partition_number // 2
    upper_partitions_half = partition_number - lower_partitions_half

    lower_half = lower_partitions_half * shade_number
    for i in range(3):
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8 / lower_half)

    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = (
                np.ones(shade_number)
                - initial_cm[
                    lower_half + j * shade_number : lower_half + (j + 1) * shade_number,
                    i,
                ]
            )
            modifier = j * modifier / upper_partitions_half
            initial_cm[
                lower_half + j * shade_number : lower_half + (j + 1) * shade_number, i
            ] += modifier

    return ListedColormap(initial_cm)
