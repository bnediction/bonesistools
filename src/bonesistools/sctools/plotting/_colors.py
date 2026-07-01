#!/usr/bin/env python

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import cycle as _cycle
from types import MappingProxyType
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import numpy as np
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.colors import to_rgba as _to_rgba

from ..._compat import Literal
from ..._validation import _as_literal, _as_positive_integer, _as_probability

RGB = Tuple[float, float, float]
RGBLike = Sequence[float]


def rgb(color: RGBLike) -> RGB:
    """
    Convert 0-255 RGB channels to 0-1 RGB channels.

    Parameters
    ----------
    color: sequence of numbers
        RGB channels in the 0-255 range.

    Returns
    -------
    tuple
        RGB channels scaled to the 0-1 range.
    """

    r, g, b = _as_rgb_channels(color)
    return (r / 255, g / 255, b / 255)


def rgba(color: Union[str, RGBLike], alpha: Optional[float] = None) -> List[float]:
    """
    Return an RGBA color.

    Parameters
    ----------
    color: str or RGB
        Color specification. RGB sequences in the 0-1 range are kept as-is.
        RGB sequences outside that range are interpreted as 0-255 channels and
        converted with `rgb`. Strings are interpreted by Matplotlib.
    alpha: float, optional
        Alpha channel in the 0-1 range. If None, use the input color alpha or
        1.0 for RGB sequences.

    Returns
    -------
    list
        RGBA channels scaled to the 0-1 range.
    """

    resolved_alpha = None if alpha is None else _as_probability(alpha, "alpha")
    if isinstance(color, str):
        return list(_to_rgba(color, alpha=resolved_alpha))

    channels = _as_rgb_channels(color)
    if not all(0 <= x <= 1 for x in channels):
        channels = rgb(channels)

    return list(channels) + [1.0 if resolved_alpha is None else resolved_alpha]


def rgb2hex(rgb: RGBLike) -> str:
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


def _hex2rgb(color: str) -> RGB:
    if color.startswith("#"):
        color = color[1:]
    return rgb([int(color[i : i + 2], 16) for i in range(0, 6, 2)])


def _as_rgb_channels(color: RGBLike) -> RGB:
    try:
        r, g, b = color
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "invalid argument value for 'color': "
            f"expected 3 RGB channels but received {color!r}"
        ) from exc

    if not all(isinstance(x, (int, float)) for x in (r, g, b)):
        raise TypeError(
            f"unsupported argument type for 'color': "
            f"expected numeric values but received {[type(x) for x in (r, g, b)]}"
        )

    return float(r), float(g), float(b)


@dataclass(frozen=True)
class Palette:
    """
    Named color palette.

    Parameters
    ----------
    name: str
        Palette name.
    colors: sequence
        RGB colors with channels scaled between 0 and 1.
    """

    name: str
    colors: Sequence[RGB]

    @property
    def hex(self) -> List[str]:
        """
        Return palette colors in hexadecimal format.
        """

        return [rgb2hex(color) for color in self.colors]

    @property
    def cmap(self) -> ListedColormap:
        """
        Return palette as a Matplotlib listed colormap.
        """

        return ListedColormap(colors=self.colors, name=self.name)

    def cycle(self) -> Iterator[RGB]:
        """
        Return an infinite color cycle over the palette.
        """

        return _cycle(self.colors)


HEX_COLORS: Dict[str, str] = {
    "black": "#000000",
    "white": "#ffffff",
    "blue": "#0014ff",
    "red": "#ff5032",
    "green": "#14c850",
    "violet": "#ff33ff",
    "lightgreen": "#14fa50",
    "coral": "#ff7f50",
    "yellow": "#ffff00",
    "darkyellow": "#cccc00",
    "lightyellow": "#808000",
    "darkorange": "#ff6900",
    "darkred": "#8b0000",
    "lightorange": "#ffa55a",
    "limegreen": "#32ff32",
    "pink": "#ffb6c1",
    "orchid": "#da70d6",
    "magenta": "#ff00ff",
    "purple": "#800080",
    "indigo": "#4b0082",
    "slateblue": "#473c8b",
    "lightgray": "#d3d3d3",
    "gray": "#708090",
    "charcoal": "#3a3a3a",
    "brown": "#8a5a2b",
    "darkgreen": "#006400",
    "gold": "#eec900",
    "amber": "#e6a100",
    "orange": "#ff8000",
    "rust": "#b84a2a",
    "salmon": "#c67171",
    "maroon": "#800000",
    "beet": "#8e388e",
    "plum": "#692a70",
    "mauve": "#8a4d91",
    "teal": "#388e8e",
    "olive": "#8e8e38",
    "navy": "#000080",
    "skyblue": "#87ceeb",
    "beige": "#ffffcc",
    "burgundy": "#6f1d1b",
    "sage": "#8a9a32",
    "moss": "#4f8a5b",
}

RGB_COLORS: Dict[str, RGB] = {
    name: _hex2rgb(color) for name, color in HEX_COLORS.items()
}

black = RGB_COLORS["black"]
white = RGB_COLORS["white"]
blue = RGB_COLORS["blue"]
red = RGB_COLORS["red"]
green = RGB_COLORS["green"]
violet = RGB_COLORS["violet"]
lightgreen = RGB_COLORS["lightgreen"]
coral = RGB_COLORS["coral"]
yellow = RGB_COLORS["yellow"]
darkred = RGB_COLORS["darkred"]
darkyellow = RGB_COLORS["darkyellow"]
lightyellow = RGB_COLORS["lightyellow"]
darkorange = RGB_COLORS["darkorange"]
lightorange = RGB_COLORS["lightorange"]
limegreen = RGB_COLORS["limegreen"]
pink = RGB_COLORS["pink"]
orchid = RGB_COLORS["orchid"]
magenta = RGB_COLORS["magenta"]
purple = RGB_COLORS["purple"]
indigo = RGB_COLORS["indigo"]
slateblue = RGB_COLORS["slateblue"]
lightgray = RGB_COLORS["lightgray"]
gray = RGB_COLORS["gray"]
charcoal = RGB_COLORS["charcoal"]
brown = RGB_COLORS["brown"]
darkgreen = RGB_COLORS["darkgreen"]
gold = RGB_COLORS["gold"]
amber = RGB_COLORS["amber"]
orange = RGB_COLORS["orange"]
rust = RGB_COLORS["rust"]
salmon = RGB_COLORS["salmon"]
maroon = RGB_COLORS["maroon"]
beet = RGB_COLORS["beet"]
plum = RGB_COLORS["plum"]
mauve = RGB_COLORS["mauve"]
teal = RGB_COLORS["teal"]
olive = RGB_COLORS["olive"]
navy = RGB_COLORS["navy"]
skyblue = RGB_COLORS["skyblue"]
beige = RGB_COLORS["beige"]
burgundy = RGB_COLORS["burgundy"]
sage = RGB_COLORS["sage"]
moss = RGB_COLORS["moss"]


def _palette(*names: str) -> List[RGB]:
    return [RGB_COLORS[name] for name in names]


_QUALITATIVE_COLOR_NAMES = (
    "blue",
    "red",
    "green",
    "orange",
    "purple",
    "pink",
    "darkred",
    "darkgreen",
    "gold",
    "indigo",
    "maroon",
)

QUALITATIVE_COLORS = _palette(*_QUALITATIVE_COLOR_NAMES)

CLASSIC_COLORS = _palette(
    *_QUALITATIVE_COLOR_NAMES,
    "skyblue",
    "teal",
    "magenta",
    "darkorange",
    "olive",
    "orchid",
    "beet",
    "navy",
    "salmon",
    "black",
    "lightgreen",
    "coral",
    "yellow",
    "limegreen",
    "slateblue",
    "darkyellow",
    "lightyellow",
    "lightorange",
    "burgundy",
)

LIGHT_COLORS = _palette(
    "skyblue",
    "red",
    "green",
    "orange",
    "orchid",
    "teal",
    "magenta",
    "violet",
    "olive",
    "beet",
    "indigo",
    "gold",
    "navy",
    "salmon",
)

EARTH_COLORS = _palette(
    "black",
    "charcoal",
    "gray",
    "brown",
    "plum",
    "mauve",
    "burgundy",
    "darkred",
    "rust",
    "orange",
    "amber",
    "gold",
    "sage",
    "moss",
    "teal",
)

_PALETTES: Dict[str, Palette] = {
    "classic": Palette(name="classic", colors=CLASSIC_COLORS),
    "light": Palette(name="light", colors=LIGHT_COLORS),
    "earth": Palette(name="earth", colors=EARTH_COLORS),
}

PALETTES: Mapping[str, Palette] = MappingProxyType(_PALETTES)

color_cycle = PALETTES["classic"].cycle()

classic_cm = PALETTES["classic"].cmap
light_cm = PALETTES["light"].cmap
earth_cm = PALETTES["earth"].cmap


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

    color_type = _as_literal(
        color_type,
        choices=("rgb", "hex"),
        name="color_type",
    )

    if color not in RGB_COLORS:
        raise ValueError(
            f"invalid argument value for 'color': color not found: {color!r}"
        )

    if color_type == "rgb":
        return RGB_COLORS[color]
    if color_type == "hex":
        return HEX_COLORS[color]


def get_palette(name: str) -> Palette:
    """
    Return a named palette.

    Parameters
    ----------
    name: str
        Palette name.

    Returns
    -------
    Palette
        Requested color palette.

    Raises
    ------
    ValueError
        If `name` is unknown.
    """

    if name not in PALETTES:
        raise ValueError(
            f"invalid argument value for 'name': palette not found: {name!r}"
        )

    return PALETTES[name]


def get_colormap(name: str) -> ListedColormap:
    """
    Return a named palette as a Matplotlib listed colormap.

    Parameters
    ----------
    name: str
        Palette name.

    Returns
    -------
    matplotlib.colors.ListedColormap
        Requested colormap.
    """

    return get_palette(name).cmap


def generate_colormap(
    color_number: int = 80,
    shade_number: Optional[int] = None,
    cm: Colormap = classic_cm,
) -> ListedColormap:
    """
    Create a colormap from another colormap by adding some new colors.

    Parameters
    ----------
    color_number: int (default: 80)
        Number of colors in the returned colormap.
    shade_number: int (optional, default: None)
        Number of shades in the returned colormap.
    cm: matplotlib.colors.Colormap (default: classic_cm)
        Initial colormap to use for creating new colormap.

    Returns
    -------
    Return a ListedColormap.

    """

    if not isinstance(cm, Colormap):
        raise TypeError(
            f"unsupported argument type for 'cm': "
            f"expected {Colormap} but received {type(cm)}"
        )

    color_number = _as_positive_integer(color_number, "color_number")
    if color_number <= cm.N:
        if isinstance(cm, ListedColormap):
            return ListedColormap(cast(Sequence[Any], cm.colors)[0:color_number])

        return ListedColormap(cm(np.linspace(0, 1, color_number)))

    if shade_number is None:
        shade_number = cm.N
    else:
        shade_number = _as_positive_integer(shade_number, "shade_number")

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
