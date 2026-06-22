#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
    overload,
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes._axes import Axes
from matplotlib.colors import Colormap, ListedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from ..._compat import Literal
from ..._warnings import _warn_deprecated
from .._typing import ScData, anndata_or_mudata_checker
from ..tools import barycenters
from ..tools._utils import (
    _UNSET,
    _resolve_representation_argument,
    get_representation,
)
from ._colors import (
    QUALITATIVE_COLORS,
    generate_colormap,
    gray,
    lightgray,
)
from ._utils import (
    _resolve_legend_options,
    colormap_colors,
    colors_from_uns,
    deprecated_bool_kwarg,
    figure_from_axes,
    normalize_color,
    qualitative_color_values,
    set_window_title,
)

Colors = Union[
    str,
    Sequence[object],
    Iterator[object],
    Colormap,
    Mapping[object, object],
]
ContinuousColors = Union[str, Colormap]
MarkerParameter = Union[float, Tuple[float, float]]


def __split_marker_parameter(
    value: MarkerParameter,
    name: str,
) -> Tuple[float, float]:

    if isinstance(value, tuple):
        if len(value) != 2:
            raise ValueError(
                f"invalid argument value for {name!r}: "
                f"expected a float or a pair of floats but received {value!r}"
            )
        return value

    return value, value


def __scatter_alpha(color: object, kwargs: Mapping[str, Any]) -> float:

    if "alpha" in kwargs:
        return cast(float, kwargs["alpha"])
    if color == lightgray:
        return 0.3
    return 1.0


def __deprecated_graph_kwarg(
    kwargs: Dict[str, Any],
    old_name: str,
    target: str,
) -> bool:

    if old_name not in kwargs:
        return False

    old_value = kwargs.pop(old_name)
    _warn_deprecated(f"`{old_name}`", replacement=f"`{target}`", stacklevel=2)

    return cast(bool, old_value)


def __continuous_colormap(colors: Optional[Colors]) -> Colormap:

    if colors is None:
        return plt.get_cmap("autumn_r")
    if isinstance(colors, str):
        return plt.get_cmap(colors)
    if isinstance(colors, Colormap):
        return colors

    raise TypeError(
        "unsupported colors for continuous observations: "
        "expected a matplotlib colormap or colormap name"
    )


def __default_plot(plot: Callable[..., Tuple[Figure, Axes]]):

    def wrapper(
        scdata: ScData,  # type: ignore
        obs: str,
        representation: str,
        colors: Optional[Colors] = None,
        n_components: Optional[int] = 2,
        ax: Optional[Axes] = None,
        **kwargs: Any,
    ) -> Tuple[Figure, Axes]:

        n_components = 2 if n_components is None else n_components
        fig, ax = cast(
            Tuple[Figure, Axes],
            plot(scdata, obs, representation, colors, n_components, ax=ax, **kwargs),
        )

        if "xlabel" in kwargs:
            ax.set_xlabel("" if kwargs["xlabel"] is None else kwargs["xlabel"])
        if "ylabel" in kwargs:
            ax.set_ylabel("" if kwargs["ylabel"] is None else kwargs["ylabel"])
        if "zlabel" in kwargs and n_components > 2:
            cast(Any, ax).set_zlabel(
                "" if kwargs["zlabel"] is None else kwargs["zlabel"]
            )

        if "tick_params" in kwargs:
            ax.tick_params(**kwargs["tick_params"])
        else:
            if "xtick_params" in kwargs:
                ax.tick_params(axis="x", **kwargs["xtick_params"])
            if "ytick_params" in kwargs:
                ax.tick_params(axis="y", **kwargs["ytick_params"])
            if n_components == 3 and "ztick_params" in kwargs:
                cast(Any, ax).tick_params(axis="z", **kwargs["ztick_params"])

        plt.sca(ax)
        (
            ax.xaxis.set_major_formatter(kwargs["formatter"])
            if "formatter" in kwargs
            else ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        )
        (
            ax.yaxis.set_major_formatter(kwargs["formatter"])
            if "formatter" in kwargs
            else ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        )
        if n_components == 3:
            ax3d = cast(Any, ax)
            (
                ax3d.zaxis.set_major_formatter(kwargs["formatter"])
                if "formatter" in kwargs
                else ax3d.zaxis.set_major_formatter(FormatStrFormatter("%g"))
            )

        if n_components == 3 and "background_visible" in kwargs:
            if kwargs["background_visible"] is False:
                ax3d = cast(Any, ax)
                ax3d.xaxis.pane.fill = False
                ax3d.yaxis.pane.fill = False
                ax3d.zaxis.pane.fill = False
                ax3d.xaxis.pane.set_edgecolor("w")
                ax3d.yaxis.pane.set_edgecolor("w")
                ax3d.zaxis.pane.set_edgecolor("w")

        return fig, ax

    return wrapper


@__default_plot
def __scatterplot_discrete(
    scdata: ScData,  # type: ignore
    obs: str,
    representation: Optional[str] = None,
    colors: Optional[Colors] = None,
    n_components: int = 2,
    legend: Optional[Dict[str, Any]] = None,
    show_legend: bool = False,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]:

    categories = set(scdata.obs[obs].unique()) & set(
        scdata.obs[obs].astype("category").cat.categories
    )

    category_values = list(scdata.obs[obs].astype("category").cat.categories)

    if colors is None:
        colors = colors_from_uns(scdata, obs, category_values)

    if colors is None:
        colors = qualitative_color_values(
            len(category_values),
            QUALITATIVE_COLORS,
            generate_colormap,
        )
    elif isinstance(colors, Mapping):
        colors = [colors[cluster] for cluster in category_values]
    if isinstance(colors, ListedColormap):
        colors = colormap_colors(colors)

    representation_mtx = cast(
        np.ndarray,
        get_representation(
            scdata,
            representation=representation,
            n_components=n_components,
        ),
    )

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(
            111, projection="rectilinear" if n_components == 2 else "3d"
        )
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
        fig.set_figwidth(
            kwargs["figwidth"]
            if "figwidth" in kwargs
            else 6 if n_components == 2 else 8
        )
    else:
        fig = figure_from_axes(ax)

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}

    if scdata.obs[obs].isna().any():
        idx = scdata.obs[obs].isna()
        if n_components == 2:
            ax.scatter(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 2.0,
                facecolors=(
                    kwargs["nan"]["facecolor"] if "facecolor" in kwargs["nan"] else gray
                ),
                edgecolors=(
                    kwargs["nan"]["edgecolor"]
                    if "edgecolor" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 1.0,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                representation_mtx[idx, 2],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 2.0,
                facecolors=(
                    kwargs["nan"]["facecolors"]
                    if "facecolors" in kwargs["nan"]
                    else gray
                ),
                edgecolors=(
                    kwargs["nan"]["facecolors"]
                    if "edgecolors" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 1.0,
            )

    color_values = cast(Sequence[object], colors)

    for _cluster, _color in zip(
        scdata.obs[obs].astype("category").cat.categories, color_values
    ):
        _color = normalize_color(_color)

        if _cluster not in categories:
            continue
        idx = np.where(scdata.obs[obs] == _cluster)[0]
        if n_components == 2:
            ax.scatter(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                s=kwargs["s"] if "s" in kwargs else 2.0,
                facecolors=_color,
                edgecolors="none",
                alpha=__scatter_alpha(_color, kwargs),
                label=_cluster,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                representation_mtx[idx, 2],
                s=kwargs["s"] if "s" in kwargs else 2.0,
                facecolors=_color,
                edgecolors="none",
                alpha=__scatter_alpha(_color, kwargs),
                label=_cluster,
            )

    if show_legend:
        box = ax.get_position()
        ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))
        legend_kwargs = {} if legend is None else dict(legend)
        if legend_kwargs:
            if "loc" not in legend_kwargs and "bbox_to_anchor" not in legend_kwargs:
                if n_components == 3:
                    fig.tight_layout()
                    fig.subplots_adjust(right=0.8)
                legend_kwargs["loc"] = "center left"
                legend_kwargs["bbox_to_anchor"] = (
                    (1.04, 0.5) if n_components == 2 else (1.09, 0.5)
                )
            else:
                pass
        else:
            legend_kwargs = {"loc": "center left", "bbox_to_anchor": (1.04, 0.5)}

        handles, labels = ax.get_legend_handles_labels()
        index = sorted(range(len(labels)), key=lambda idx: labels[idx])
        handles = [handles[i] for i in index]
        labels = [labels[i] for i in index]
        ax.legend(handles, labels, **legend_kwargs)

    return fig, ax


@__default_plot
def __scatterplot_continuous(
    scdata: ScData,  # type: ignore
    obs: str,
    representation: Optional[str] = None,
    colors: Optional[ContinuousColors] = None,
    n_components: int = 2,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]:

    cmap = __continuous_colormap(colors)

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}
    if "facecolor" in kwargs["nan"] and "color" not in kwargs["nan"]:
        kwargs["nan"]["color"] = kwargs["nan"]["facecolor"]
    if not np.all(cmap.get_bad() != 0) and "facecolor" not in kwargs["nan"]:
        kwargs["nan"]["color"] = cmap.get_bad()

    representation_mtx = cast(
        np.ndarray,
        get_representation(
            scdata,
            representation=representation,
            n_components=n_components,
        ),
    )

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(
            111, projection="rectilinear" if n_components == 2 else "3d"
        )
        fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
        fig.set_figwidth(
            kwargs["figwidth"]
            if "figwidth" in kwargs
            else 6 if n_components == 2 else 8
        )
    else:
        fig = figure_from_axes(ax)

    if scdata.obs[obs].isna().any():
        idx = scdata.obs[obs].isna()
        if n_components == 2:
            ax.scatter(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 2.0,
                facecolors=(
                    kwargs["nan"]["facecolor"]
                    if "facecolor" in kwargs["nan"]
                    else lightgray
                ),
                edgecolors=(
                    kwargs["nan"]["edgecolor"]
                    if "edgecolor" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 1.0,
            )
        elif n_components == 3:
            ax3d = cast(Any, ax)
            ax3d.scatter3D(
                representation_mtx[idx, 0],
                representation_mtx[idx, 1],
                representation_mtx[idx, 2],
                s=kwargs["nan"]["s"] if "s" in kwargs["nan"] else 2.0,
                facecolors=(
                    kwargs["nan"]["facecolors"]
                    if "facecolors" in kwargs["nan"]
                    else gray
                ),
                edgecolors=(
                    kwargs["nan"]["facecolors"]
                    if "edgecolors" in kwargs["nan"]
                    else "none"
                ),
                alpha=kwargs["nan"]["alpha"] if "alpha" in kwargs["nan"] else 1.0,
            )

    if n_components == 2:
        sc = ax.scatter(
            representation_mtx[:, 0],
            representation_mtx[:, 1],
            s=kwargs["s"] if "s" in kwargs else 2.0,
            c=scdata.obs[obs],
            cmap=cmap,
            vmin=kwargs["vmin"] if "vmin" in kwargs else None,
            vmax=kwargs["vmax"] if "vmax" in kwargs else None,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1.0,
        )
        cb = fig.colorbar(
            sc, shrink=kwargs["colorbar_scale"] if "colorbar_scale" in kwargs else 1
        )
        cb.update_ticks()
    elif n_components == 3:
        ax3d = cast(Any, ax)
        ax3d.scatter3D(
            representation_mtx[:, 0],
            representation_mtx[:, 1],
            representation_mtx[:, 2],
            s=kwargs["s"] if "s" in kwargs else 2.0,
            c=scdata.obs[obs],
            cmap=cmap,
            vmin=kwargs["vmin"] if "vmin" in kwargs else None,
            vmax=kwargs["vmax"] if "vmax" in kwargs else None,
            edgecolors="none",
            alpha=kwargs["alpha"] if "alpha" in kwargs else 1.0,
        )

    return fig, ax


def __draw_labels(
    scdata: ScData,  # type: ignore
    obs: str,
    representation: Optional[str] = None,
    ax: Optional[Axes] = None,
    dim: int = 2,
    **kwargs: Any,
) -> None:

    barycenter_values = cast(
        Dict[object, np.ndarray],
        barycenters(scdata=scdata, obs=obs, representation=representation),
    )
    keys_to_remove = []
    for k, v in barycenter_values.items():
        if np.isnan(v).any():
            keys_to_remove.append(k)
    for k in keys_to_remove:
        barycenter_values.pop(k, None)

    if ax is None:
        ax = plt.gca()

    if "verticalalignment" not in kwargs:
        kwargs["verticalalignment"] = "center"
    else:
        pass

    if dim == 2:
        for label, value in barycenter_values.items():
            ax.text(x=value[0], y=value[1], s=str(label), **kwargs)
    elif dim == 3:
        ax3d = cast(Any, ax)
        for label, value in barycenter_values.items():
            ax3d.text(x=value[0], y=value[1], z=value[2], s=label, **kwargs)


@overload
def embedding(
    scdata: ScData,
    obs: str,
    representation: Any = _UNSET,
    *,
    n_components: Literal[2, 3] = 2,
    s: MarkerParameter = (2.0, 2.0),
    alpha: MarkerParameter = (1.0, 1.0),
    colors: Optional[Colors] = None,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Optional[Dict[str, Any]] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    ax: Optional[Axes] = None,
    outfile: None = None,
    use_rep: Any = _UNSET,
    obsm: Any = _UNSET,
    **kwargs: Any,
) -> Tuple[Figure, Axes]: ...


@overload
def embedding(
    scdata: ScData,
    obs: str,
    representation: Any = _UNSET,
    *,
    n_components: Literal[2, 3] = 2,
    s: MarkerParameter = (2.0, 2.0),
    alpha: MarkerParameter = (1.0, 1.0),
    colors: Optional[Colors] = None,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Optional[Dict[str, Any]] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    ax: Optional[Axes] = None,
    outfile: Path,
    use_rep: Any = _UNSET,
    obsm: Any = _UNSET,
    **kwargs: Any,
) -> None: ...


@overload
def embedding(
    scdata: ScData,
    obs: str,
    representation: Any = _UNSET,
    *,
    n_components: Literal[2, 3] = 2,
    s: MarkerParameter = (2.0, 2.0),
    alpha: MarkerParameter = (1.0, 1.0),
    colors: Optional[Colors] = None,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Optional[Dict[str, Any]] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    use_rep: Any = _UNSET,
    obsm: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]: ...


@anndata_or_mudata_checker
def embedding(
    scdata: ScData,  # type: ignore
    obs: str,
    representation: Any = _UNSET,
    *,
    n_components: Literal[2, 3] = 2,
    s: MarkerParameter = (2.0, 2.0),
    alpha: MarkerParameter = (1.0, 1.0),
    colors: Optional[Colors] = None,
    title: Optional[Union[str, Dict[str, Any]]] = None,
    legend: Optional[Dict[str, Any]] = None,
    show_legend: bool = True,
    show_labels: bool = False,
    ax: Optional[Axes] = None,
    outfile: Optional[Path] = None,
    use_rep: Any = _UNSET,
    obsm: Any = _UNSET,
    **kwargs: Any,
) -> Optional[Tuple[Figure, Axes]]:
    """
    Draw an embedding stored in `scdata.obsm`.

    Parameters
    ----------
    scdata: AnnData or MuData
        Unimodal or multimodal annotated data matrix.
    obs: str
        Observation column in `scdata.obs` defining groups or values.
    representation: str
        Representation key in `scdata.obsm`.
    n_components: 2 or 3 (default: 2)
        Number of plotted dimensions.
    s: float or tuple[float, float] (default: (2, 2))
        Marker size. If a pair is provided, the first value applies to valid
        observations and the second to observations with missing values in
        `obs`.
    alpha: float or tuple[float, float] (default: (1, 1))
        Marker transparency. If a pair is provided, the first value applies to
        valid observations and the second to observations with missing values
        in `obs`.
    colors: matplotlib.colors.Colormap, optional
        Colormap or colors used to draw observations.
    title: str or dict, optional
        Figure title, or keyword arguments passed to `Axes.set_title`.
    legend: dict, optional
        Keyword arguments passed to `Axes.legend` when `show_legend=True`.
    show_legend: bool (default: True)
        Draw the legend when `scdata.obs[obs]` contains discrete values.
    show_labels: bool (default: False)
        Draw labels retrieved from `scdata.obs[obs]` on the embedding.
    ax: matplotlib.axes.Axes, optional
        Axes on which to draw the embedding. If not provided, create a new
        figure and axes.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    use_rep: str, optional
        Deprecated alias for `representation`.
    obsm: str, optional
        Deprecated alias for `representation`.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
        - xlabel[str]: set the label for the x-axis
        - ylabel[str]: set the label for the y-axis
        - zlabel[str]: set the label for the z-axis
        - formatter[matplotlib.ticker.FormatStrFormatter]: specify the major
          formatter on x-, y- and z-axis
        - tick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines following the syntax of matplotlib.axes.Axes.tick_params
        - xtick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on x-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - ytick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on y-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - ztick_params[dict]: change the appearance of ticks, tick labels, and
          gridlines on z-axis following the syntax of
          matplotlib.axes.Axes.tick_params
        - text[dict]: change the appearance of text in figure following the
          syntax of matplotlib.text
        - background_visible[bool]: specify if background color is visible or
          not in case of 3D plotting
        - automatic_resize[bool]: resize figure to accommodate large legends
        - vmin[float]: lower bound for continuous color scaling
        - vmax[float]: upper bound for continuous color scaling

    Returns
    -------
    tuple[Figure, Axes] or None
        Figure and axes if `outfile` is None; otherwise None after saving the
        figure.

    Raises
    ------
    ValueError
        If `n_components` is not 2 or 3.
    """

    representation = _resolve_representation_argument(
        representation,
        use_rep,
        default=None,
        stacklevel=2,
        obsm=obsm,
    )
    if representation is None:
        raise TypeError("missing required argument: 'representation'")

    if n_components not in [2, 3]:
        raise ValueError(
            f"invalid argument value for 'n_components': "
            f"expected 2 or 3 but received {n_components!r}"
        )
    automatic_resize = bool(kwargs.pop("automatic_resize", False))
    show_legend, legend = _resolve_legend_options(
        show_legend,
        legend,
        kwargs,
        deprecated_show_name="add_legend",
    )

    if "alpha" in kwargs:
        raise TypeError(
            "invalid argument combination: use either 'alpha' "
            "or '**kwargs[\"alpha\"]', not both"
        )
    if "s" in kwargs:
        raise TypeError(
            "invalid argument combination: use either 's' "
            "or '**kwargs[\"s\"]', not both"
        )
    valid_s, missing_s = __split_marker_parameter(s, "s")
    valid_alpha, missing_alpha = __split_marker_parameter(alpha, "alpha")
    kwargs["s"] = valid_s
    kwargs["alpha"] = valid_alpha
    kwargs.setdefault("nan", {})
    kwargs["nan"].setdefault("s", missing_s)
    kwargs["nan"].setdefault("alpha", missing_alpha)

    component_number = n_components
    show_labels = deprecated_bool_kwarg(
        kwargs,
        "add_labels",
        "show_labels",
        show_labels,
        False,
    )
    draw_deprecated_graph = __deprecated_graph_kwarg(
        kwargs,
        "add_graph",
        "bt.sct.pl.trajectory(..., graph_key=...)",
    )
    draw_deprecated_graph_labels = __deprecated_graph_kwarg(
        kwargs,
        "add_labels_to_graph",
        "bt.sct.pl.trajectory(..., show_labels=...)",
    )

    if pd.api.types.is_float_dtype(scdata.obs[obs]):
        fig, ax = __scatterplot_continuous(
            scdata,
            obs,
            representation,
            cast(Optional[ContinuousColors], colors),
            component_number,
            ax=ax,
            **kwargs,
        )
    elif (
        pd.api.types.is_integer_dtype(scdata.obs[obs])
        or pd.api.types.is_bool_dtype(scdata.obs[obs])
        or pd.api.types.is_string_dtype(scdata.obs[obs])
        or isinstance(scdata.obs[obs].dtype, pd.CategoricalDtype)
    ):
        fig, ax = __scatterplot_discrete(
            scdata,
            obs,
            representation,
            colors,
            component_number,
            legend=legend,
            ax=ax,
            show_legend=show_legend,
            **kwargs,
        )
    else:
        raise TypeError(
            f"unsupported dtype for observation {obs!r}: " f"{scdata.obs[obs].dtype!r}"
        )

    if show_labels:
        _kwargs = {} if "text" not in kwargs else kwargs["text"]
        __draw_labels(
            scdata,
            obs=obs,
            representation=representation,
            ax=ax,
            dim=component_number,
            **_kwargs,
        )

    if draw_deprecated_graph or draw_deprecated_graph_labels:
        from ._graph import graph_overlay

        _kwargs = {"linewidth": 2.5, "zorder": 10}
        _kwargs.update({} if "graph" not in kwargs else kwargs["graph"])
        graph_overlay(
            scdata,
            ax=ax,
            n_components=component_number,
            show_labels=draw_deprecated_graph_labels,
            z_offset=kwargs["graph_z_offset"] if "graph_z_offset" in kwargs else 0.0,
            label_kwargs=kwargs["text"] if "text" in kwargs else None,
            **_kwargs,
        )

    if automatic_resize:
        fig.tight_layout(pad=1.2, rect=(0, 0, 0.84, 1))

    if title:
        if isinstance(title, str):
            set_window_title(fig, title)
            ax.set_title(title)
        elif isinstance(title, dict):
            set_window_title(fig, title["label"])
            ax.set_title(**title)
        else:
            raise TypeError(
                f"unsupported argument type for 'title': "
                f"expected {str} or {dict} but received {type(title)}"
            )

    if outfile:
        plt.savefig(outfile, bbox_inches="tight")
        plt.close()
        return None
    else:
        return fig, ax


def embedding_plot(*args: Any, **kwargs: Any) -> Optional[Tuple[Figure, Axes]]:
    """
    Deprecated alias for `embedding()`.
    """

    _warn_deprecated(
        "`bt.sct.pl.embedding_plot`",
        replacement="`bt.sct.pl.embedding`",
        stacklevel=2,
    )

    return embedding(*args, **kwargs)
