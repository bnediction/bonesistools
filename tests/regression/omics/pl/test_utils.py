#!/usr/bin/env python

from typing import Any, Dict, cast

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import ListedColormap

from bonesistools.omics.plotting import _utils


def test_set_axis_label_accepts_public_label_forms_and_validates_errors():
    fig, ax = plt.subplots()

    _utils.set_axis_label(ax, "xlabel", "x")
    assert ax.get_xlabel() == "x"

    _utils.set_axis_label(ax, "ylabel", {"label": "y", "fontsize": 13})
    assert ax.get_ylabel() == "y"
    assert ax.yaxis.label.get_fontsize() == 13

    _utils.set_axis_label(ax, "xlabel", {"xlabel": "x2"})
    assert ax.get_xlabel() == "x2"

    _utils.set_axis_label(ax, "xlabel", {"label": None})
    assert ax.get_xlabel() == ""

    with pytest.raises(ValueError, match="axis label name"):
        _utils.set_axis_label(ax, "title", "bad")

    with pytest.raises(TypeError, match="unsupported argument type"):
        _utils.set_axis_label(ax, "xlabel", cast(Any, 1))

    with pytest.raises(TypeError, match="use either 'label' or 'xlabel'"):
        _utils.set_axis_label(ax, "xlabel", {"label": "x", "xlabel": "x"})

    with pytest.raises(ValueError, match="'label' key"):
        _utils.set_axis_label(ax, "xlabel", {"fontsize": 12})

    with pytest.raises(TypeError, match="unsupported label type"):
        _utils.set_axis_label(ax, "xlabel", {"label": 1})

    plt.close(fig)


def test_apply_legend_draws_with_handles_and_can_remove_existing_legend():
    fig, ax = plt.subplots()
    handle = ax.scatter([0], [0], label="A")

    _utils.apply_legend(ax, {"title": "group"}, handles=[handle], labels=["A"])
    legend = ax.get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "group"

    _utils.apply_legend(ax, False, remove=True)
    assert ax.get_legend() is None

    plt.close(fig)


def test_deprecated_bool_kwarg_routes_value_and_rejects_conflicts():
    assert (
        _utils.deprecated_bool_kwarg(
            {},
            "old",
            "new",
            new_value=True,
            default_value=True,
        )
        is True
    )

    kwargs: Dict[str, Any] = {"old": False}

    with pytest.warns(FutureWarning, match="`old` is deprecated"):
        value = _utils.deprecated_bool_kwarg(
            kwargs,
            "old",
            "new",
            new_value=True,
            default_value=True,
        )

    assert value is False
    assert kwargs == {}

    with pytest.warns(FutureWarning, match="`old` is deprecated"):
        with pytest.raises(TypeError, match="use either 'old' or 'new'"):
            _utils.deprecated_bool_kwarg(
                {"old": False},
                "old",
                "new",
                new_value=False,
                default_value=True,
            )


def test_deprecated_kwarg_helpers_validate_conflicting_names():
    empty_kwargs: Dict[str, Any] = {}
    assert _utils.rename_deprecated_kwarg(empty_kwargs, "old", "new") is None
    assert empty_kwargs == {}

    kwargs: Dict[str, Any] = {"old": 1}

    with pytest.warns(FutureWarning, match="`old` is deprecated"):
        _utils.rename_deprecated_kwarg(kwargs, "old", "new")

    assert kwargs == {"new": 1}

    with pytest.raises(TypeError, match="use either 'old' or 'new'"):
        _utils.rename_deprecated_kwarg({"old": 1, "new": 2}, "old", "new")

    bool_kwargs: Dict[str, Any] = {"old_bool": True}
    with pytest.warns(FutureWarning, match="`old_bool` is deprecated"):
        assert (
            _utils._resolve_bool_kwarg(
                bool_kwargs,
                "new_bool",
                "old_bool",
                default=False,
            )
            is True
        )

    with pytest.warns(FutureWarning, match="`old_bool` is deprecated"):
        with pytest.raises(TypeError, match="use either 'old_bool' or 'new_bool'"):
            _utils._resolve_bool_kwarg(
                {"old_bool": True, "new_bool": False},
                "new_bool",
                "old_bool",
                default=False,
            )


def test_legend_and_toggle_mapping_resolution_validate_public_contract():
    with pytest.warns(DeprecationWarning, match="'show_legend' is deprecated"):
        assert (
            _utils._resolve_legend_argument(
                True,
                {"show_legend": False},
            )
            is True
        )

    with pytest.warns(FutureWarning, match="`legend_params` is deprecated"):
        assert _utils._resolve_legend_argument(
            True,
            {"legend_params": {"loc": "upper left"}},
        ) == {"loc": "upper left"}

    with pytest.warns(FutureWarning, match="`legend_params` is deprecated"):
        with pytest.raises(TypeError, match="legend"):
            _utils._resolve_legend_argument(
                {"loc": "upper right"},
                {"legend_params": {"loc": "upper left"}},
            )

    with pytest.raises(TypeError, match="unsupported argument type for 'legend'"):
        _utils._resolve_legend_argument(cast(Any, 1), {})

    with pytest.warns(FutureWarning, match="`show_points` is deprecated"):
        enabled, params = _utils._resolve_toggle_mapping_argument(
            False,
            {"show_points": {"s": 3}},
            name="points",
            deprecated_names=("show_points",),
            default=False,
        )

    assert enabled is True
    assert params == {"s": 3}

    with pytest.warns(FutureWarning, match="`show_points` is deprecated"):
        with pytest.raises(TypeError, match="points"):
            _utils._resolve_toggle_mapping_argument(
                True,
                {"show_points": False},
                name="points",
                deprecated_names=("show_points",),
                default=False,
            )

    with pytest.raises(TypeError, match="unsupported argument type for 'points'"):
        _utils._resolve_toggle_mapping_argument(cast(Any, 1), {}, name="points")


def test_color_utility_helpers_handle_palettes_and_uns_metadata(mini_adata):
    colormap = ListedColormap(["red", "blue"])
    assert list(_utils.colormap_colors(colormap)) == ["red", "blue"]
    assert list(_utils.colormap_colors(["black"])) == ["black"]

    assert _utils.normalize_color((1.0, 0.0, 0.0)) == [1.0, 0.0, 0.0, 1]
    assert _utils.normalize_color("red") == "red"

    generated = _utils.qualitative_color_values(
        3,
        ["red"],
        lambda color_number: ["c0", "c1", "c2"][:color_number],
    )
    assert generated == ["c0", "c1", "c2"]
    assert _utils.qualitative_color_values(1, ["red"], cast(Any, None)) == ["red"]

    mini_adata.uns["cluster_color"] = {"A": "red", "B": "blue"}
    assert _utils.colors_from_uns(mini_adata, "cluster", ["A", "B"]) == {
        "A": "red",
        "B": "blue",
    }

    mini_adata.uns.pop("cluster_color")
    mini_adata.uns["cluster_colors"] = ["red"]
    assert _utils.colors_from_uns(mini_adata, "cluster", ["A", "B"]) == {"A": "red"}
    assert _utils.colors_from_uns(mini_adata, "missing", ["A"]) is None
