#!/usr/bin/env python

import math
from typing import Any, cast

import pytest
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

import bonesistools as bt
from bonesistools.sctools.plotting import _colors


def test_rgb_and_rgb2hex_convert_channels():
    assert bt.sct.pl.rgb([0, 128, 255]) == (0.0, 128 / 255, 1.0)
    assert bt.sct.pl.rgba([0, 128, 255], alpha=0.5) == [
        0.0,
        128 / 255,
        1.0,
        0.5,
    ]
    assert bt.sct.pl.rgba([0.0, 0.5, 1.0]) == [0.0, 0.5, 1.0, 1.0]
    assert bt.sct.pl.rgba("red", alpha=0.5) == [1.0, 0.0, 0.0, 0.5]
    assert bt.sct.pl.rgb2hex([1.0, 0.5, 0.0]) == "#ff8000"
    assert bt.sct.pl.rgb2hex([255, 128, 0]) == "#ff8000"


def test_rgb2hex_and_rgba_reject_invalid_channels():
    with pytest.raises(TypeError, match="unsupported argument type for 'rgb'"):
        bt.sct.pl.rgb2hex(cast(Any, ["red", 0, 0]))
    with pytest.raises(TypeError, match="unsupported argument type for 'color'"):
        bt.sct.pl.rgba(cast(Any, ["red", 0, 0]))
    with pytest.raises(ValueError, match="expected 3 RGB channels"):
        bt.sct.pl.rgba([0, 0, 0, 0])
    with pytest.raises(ValueError, match="expected value between 0 and 1"):
        bt.sct.pl.rgba([0, 0, 0], alpha=2.0)


def test_get_color_returns_rgb_or_hex_and_rejects_invalid_values():
    assert bt.sct.pl.get_color("black") == (0.0, 0.0, 0.0)
    assert bt.sct.pl.get_color("black", color_type="hex") == "#000000"
    assert bt.sct.pl.get_color("charcoal", color_type="hex") == "#3a3a3a"
    assert bt.sct.pl.get_color("burgundy", color_type="hex") == "#6f1d1b"
    assert _colors.HEX_COLORS["burgundy"] == "#6f1d1b"
    assert bt.sct.pl.rgb2hex(_colors.RGB_COLORS["burgundy"]) == "#6f1d1b"

    with pytest.raises(ValueError, match="invalid argument value for 'color_type'"):
        bt.sct.pl.get_color("black", color_type=cast(Any, "hsl"))

    with pytest.raises(ValueError, match="invalid argument value for 'color'"):
        bt.sct.pl.get_color("not_a_color")


def test_earth_palette_matches_expected_hex_values():
    expected_hex = [
        "#000000",
        "#3a3a3a",
        "#708090",
        "#8a5a2b",
        "#692a70",
        "#8a4d91",
        "#6f1d1b",
        "#8b0000",
        "#b84a2a",
        "#ff8000",
        "#e6a100",
        "#eec900",
        "#8a9a32",
        "#4f8a5b",
        "#388e8e",
    ]

    assert [
        bt.sct.pl.rgb2hex(color) for color in _colors.EARTH_COLORS
    ] == expected_hex
    assert bt.sct.pl.get_palette("earth").hex == expected_hex
    assert bt.sct.pl.get_colormap("earth").N == len(_colors.EARTH_COLORS)


def test_classic_palette_starts_with_qualitative_colors_without_duplicates():
    classic_hex = bt.sct.pl.get_palette("classic").hex
    qualitative_hex = [
        bt.sct.pl.rgb2hex(color) for color in _colors.QUALITATIVE_COLORS
    ]

    assert classic_hex[: len(qualitative_hex)] == qualitative_hex
    assert len(classic_hex) == len(set(classic_hex))


def test_palette_registry_returns_palette_and_colormap():
    palette = bt.sct.pl.get_palette("earth")
    colormap = bt.sct.pl.get_colormap("earth")

    assert isinstance(palette, _colors.Palette)
    assert palette.name == "earth"
    assert palette.colors == _colors.EARTH_COLORS
    assert next(palette.cycle()) == _colors.EARTH_COLORS[0]
    assert isinstance(colormap, ListedColormap)
    assert colormap.name == "earth"
    assert colormap.N == len(_colors.EARTH_COLORS)
    assert bt.sct.pl.get_palette("classic").colors == _colors.CLASSIC_COLORS
    assert bt.sct.pl.get_palette("light").colors == _colors.LIGHT_COLORS
    assert bt.sct.pl.get_colormap("light").N == len(_colors.LIGHT_COLORS)

    with pytest.raises(ValueError, match="invalid argument value for 'name'"):
        bt.sct.pl.get_palette("unknown")


def test_palette_internals_are_not_exposed_publicly():
    assert "PALETTES" not in dir(bt.sct.pl)
    assert not hasattr(bt.sct.pl, "PALETTES")
    for name in [
        "CLASSIC_COLORS",
        "LIGHT_COLORS",
        "EARTH_COLORS",
        "QUALITATIVE_COLORS",
        "HEX_COLORS",
        "RGB_COLORS",
        "Palette",
        "color_cycle",
        "classic_cm",
        "light_cm",
        "earth_cm",
    ]:
        assert name not in dir(bt.sct.pl)
        assert not hasattr(bt.sct.pl, name)


def test_generate_colormap_shorter_and_longer_than_base_colormap():
    shorter = bt.sct.pl.generate_colormap(color_number=3)
    continuous_cm = LinearSegmentedColormap.from_list(
        "continuous",
        ["black", "white"],
    )
    shorter_from_continuous = bt.sct.pl.generate_colormap(
        color_number=3,
        cm=continuous_cm,
    )
    longer = bt.sct.pl.generate_colormap(color_number=40, shade_number=5)
    longer_with_default_shades = bt.sct.pl.generate_colormap(color_number=40)

    assert isinstance(shorter, ListedColormap)
    assert shorter.N == 3
    assert isinstance(shorter_from_continuous, ListedColormap)
    assert shorter_from_continuous.N == 3
    assert isinstance(longer, ListedColormap)
    assert longer.N == 40
    assert isinstance(longer_with_default_shades, ListedColormap)
    assert longer_with_default_shades.N == (
        math.ceil(40 / _colors.classic_cm.N) * _colors.classic_cm.N
    )


def test_generate_colormap_validates_arguments():
    with pytest.raises(TypeError, match="unsupported argument type for 'cm'"):
        bt.sct.pl.generate_colormap(cm=cast(Any, "not a colormap"))

    with pytest.raises(ValueError, match="invalid argument value for 'color_number'"):
        bt.sct.pl.generate_colormap(color_number=0)

    with pytest.raises(ValueError, match="invalid argument value for 'shade_number'"):
        bt.sct.pl.generate_colormap(color_number=40, shade_number=0)
