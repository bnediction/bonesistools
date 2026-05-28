#!/usr/bin/env python

import math

import pytest
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

import bonesistools as bt


def test_rgb_and_rgb2hex_convert_channels():
    assert bt.sct.pl.rgb([0, 128, 255]) == [0.0, 128 / 255, 1.0]
    assert bt.sct.pl.rgb2hex([1.0, 0.5, 0.0]) == "#ff8000"
    assert bt.sct.pl.rgb2hex([255, 128, 0]) == "#ff8000"


def test_rgb2hex_rejects_non_numeric_channels():
    with pytest.raises(TypeError, match="unsupported argument type for 'rgb'"):
        bt.sct.pl.rgb2hex(["red", 0, 0])


def test_get_color_returns_rgb_or_hex_and_rejects_invalid_values():
    assert bt.sct.pl.get_color("black") == [0.0, 0.0, 0.0]
    assert bt.sct.pl.get_color("black", color_type="hex") == "#000000"

    with pytest.raises(ValueError, match="invalid argument value for 'color_type'"):
        bt.sct.pl.get_color("black", color_type="hsl")

    with pytest.raises(ValueError, match="invalid argument value for 'color'"):
        bt.sct.pl.get_color("not_a_color")


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
        math.ceil(40 / bt.sct.pl.bonesis_cm.N) * bt.sct.pl.bonesis_cm.N
    )


def test_generate_colormap_validates_arguments():
    with pytest.raises(TypeError, match="unsupported argument type for 'cm'"):
        bt.sct.pl.generate_colormap(cm="not a colormap")

    with pytest.raises(ValueError, match="invalid argument value for 'color_number'"):
        bt.sct.pl.generate_colormap(color_number=0)

    with pytest.raises(ValueError, match="invalid argument value for 'shade_number'"):
        bt.sct.pl.generate_colormap(color_number=40, shade_number=0)
