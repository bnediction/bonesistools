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

from ..._compat import Literal
from ..._warnings import _warn_deprecated, _warn_deprecated_argument
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
    _resolve_legend_argument,
    _resolve_toggle_mapping_argument,
    apply_legend,
    colormap_colors,
    colors_from_uns,
    normalize_color,
    qualitative_color_values,
    save_figure,
    set_axis_formatters,
    set_axis_label,
    set_title,
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
    legend: Union[bool, Mapping[str, Any]] = True,
    labels: Union[bool, Mapping[str, Any]] = False,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    zlabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    legend: Union[bool, Mapping[str, Any]] = True,
    labels: Union[bool, Mapping[str, Any]] = False,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    zlabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    legend: Union[bool, Mapping[str, Any]] = True,
    labels: Union[bool, Mapping[str, Any]] = False,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    zlabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    legend: Union[bool, Mapping[str, Any]] = True,
    labels: Union[bool, Mapping[str, Any]] = False,
    xlabel: Optional[Union[str, Mapping[str, Any]]] = None,
    ylabel: Optional[Union[str, Mapping[str, Any]]] = None,
    zlabel: Optional[Union[str, Mapping[str, Any]]] = None,
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
    legend: bool or mapping (default: True)
        Legend configuration. `False` disables the legend. `True` draws the legend
        using default Matplotlib parameters. If a mapping is provided, it is
        forwarded as keyword arguments to `Axes.legend`.
    labels: bool or mapping (default: False)
        Label configuration. `False` disables labels. `True` draws labels
        retrieved from `scdata.obs[obs]` using default Matplotlib text
        parameters. If a mapping is provided, it is forwarded as keyword
        arguments to `Axes.text`.
    xlabel: str or mapping, optional
        X-axis label. If a mapping is provided, it is forwarded as keyword
        arguments to `Axes.set_xlabel` and must contain a `label` key.
    ylabel: str or mapping, optional
        Y-axis label. If a mapping is provided, it is forwarded as keyword
        arguments to `Axes.set_ylabel` and must contain a `label` key.
    zlabel: str or mapping, optional
        Z-axis label for three-dimensional embeddings. If a mapping is
        provided, it is forwarded as keyword arguments to `Axes.set_zlabel` and
        must contain a `label` key.
    ax: matplotlib.axes.Axes, optional
        Axes on which to draw the embedding. If not provided, create a new
        figure and axes.
    outfile: Path, optional
        If specified, save the figure instead of returning it.
    **kwargs: Any
        Supplemental features for figure plotting:
        - figheight[float]: specify the figure height
        - figwidth[float]: specify the figure width
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
    legend = _resolve_legend_argument(legend, kwargs, stacklevel=2)
    draw_labels, label_kwargs = _resolve_toggle_mapping_argument(
        labels,
        kwargs,
        name="labels",
        deprecated_names=("show_labels", "add_labels"),
        stacklevel=2,
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
    if xlabel is not None:
        kwargs["xlabel"] = xlabel
    if ylabel is not None:
        kwargs["ylabel"] = ylabel
    if zlabel is not None:
        kwargs["zlabel"] = zlabel

    component_number = n_components
    draw_deprecated_graph = __deprecated_graph_kwarg(
        kwargs,
        "add_graph",
        "bt.sct.pl.trajectory(..., graph_key=...)",
    )
    draw_deprecated_graph_labels = __deprecated_graph_kwarg(
        kwargs,
        "add_labels_to_graph",
        "bt.sct.pl.trajectory(..., labels=...)",
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
            **kwargs,
        )
    else:
        raise TypeError(
            f"unsupported dtype for observation {obs!r}: " f"{scdata.obs[obs].dtype!r}"
        )

    if draw_labels:
        if "text" in kwargs:
            _warn_deprecated_argument("text", "labels", stacklevel=2)
            label_kwargs = {
                **cast(Mapping[str, Any], kwargs.pop("text")),
                **label_kwargs,
            }
        __draw_labels(
            scdata,
            obs=obs,
            representation=representation,
            ax=ax,
            dim=component_number,
            **label_kwargs,
        )

    if draw_deprecated_graph or draw_deprecated_graph_labels:
        from ._graph import graph_overlay

        _kwargs = {"linewidth": 2.5, "zorder": 10}
        _kwargs.update({} if "graph" not in kwargs else kwargs["graph"])
        graph_overlay(
            scdata,
            ax=ax,
            n_components=component_number,
            labels=draw_deprecated_graph_labels,
            z_offset=kwargs["graph_z_offset"] if "graph_z_offset" in kwargs else 0.0,
            **_kwargs,
        )

    if automatic_resize:
        fig.tight_layout(pad=1.2, rect=(0, 0, 0.84, 1))

    set_title(fig, ax, title)

    if save_figure(fig, outfile):
        return None

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


def __discrete_category_values(series: pd.Series) -> Sequence[object]:
    categorical = series.astype("category")
    values = list(categorical.cat.categories)
    if isinstance(series.dtype, pd.CategoricalDtype) and series.cat.ordered:
        return values

    try:
        return sorted(values)
    except TypeError:
        return sorted(values, key=lambda value: str(value))


def __embedding_axes(
    ax: Optional[Axes],
    n_components: int,
    kwargs: Mapping[str, Any],
) -> Tuple[Figure, Axes]:

    if ax is not None:
        return cast(Figure, ax.figure), ax

    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        projection="rectilinear" if n_components == 2 else "3d",
    )
    fig.set_figheight(kwargs["figheight"] if "figheight" in kwargs else 5)
    fig.set_figwidth(
        kwargs["figwidth"] if "figwidth" in kwargs else 6 if n_components == 2 else 8
    )
    return fig, ax


def __embedding_matrix(
    scdata: ScData,
    representation: Optional[str],
    n_components: int,
) -> np.ndarray:

    return cast(
        np.ndarray,
        get_representation(
            scdata,
            obsm=representation,
            n_components=n_components,
        ),
    )


def __scatter_values(
    ax: Axes,
    representation_mtx: np.ndarray,
    rows: Any,
    n_components: int,
    **kwargs: Any,
) -> Any:

    if n_components == 2:
        return ax.scatter(
            representation_mtx[rows, 0],
            representation_mtx[rows, 1],
            **kwargs,
        )

    ax3d = cast(Any, ax)
    return ax3d.scatter3D(
        representation_mtx[rows, 0],
        representation_mtx[rows, 1],
        representation_mtx[rows, 2],
        **kwargs,
    )


def __scatter_missing_values(
    scdata: ScData,
    obs: str,
    representation_mtx: np.ndarray,
    ax: Axes,
    n_components: int,
    kwargs: Mapping[str, Any],
    default_color: object,
) -> None:

    if not scdata.obs[obs].isna().any():
        return None

    nan_kwargs = cast(Mapping[str, Any], kwargs["nan"])
    idx = scdata.obs[obs].isna()
    if n_components == 2:
        __scatter_values(
            ax,
            representation_mtx,
            idx,
            n_components,
            s=nan_kwargs["s"] if "s" in nan_kwargs else 2.0,
            facecolors=(
                nan_kwargs["facecolor"] if "facecolor" in nan_kwargs else default_color
            ),
            edgecolors=(
                nan_kwargs["edgecolor"] if "edgecolor" in nan_kwargs else "none"
            ),
            alpha=nan_kwargs["alpha"] if "alpha" in nan_kwargs else 1.0,
        )
        return None

    __scatter_values(
        ax,
        representation_mtx,
        idx,
        n_components,
        s=nan_kwargs["s"] if "s" in nan_kwargs else 2.0,
        facecolors=(
            nan_kwargs["facecolors"] if "facecolors" in nan_kwargs else default_color
        ),
        edgecolors=nan_kwargs["facecolors"] if "edgecolors" in nan_kwargs else "none",
        alpha=nan_kwargs["alpha"] if "alpha" in nan_kwargs else 1.0,
    )
    return None


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
            set_axis_label(ax, "xlabel", kwargs["xlabel"])
        if "ylabel" in kwargs:
            set_axis_label(ax, "ylabel", kwargs["ylabel"])
        if "zlabel" in kwargs and n_components > 2:
            set_axis_label(ax, "zlabel", kwargs["zlabel"])

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
        set_axis_formatters(
            ax,
            kwargs.get("formatter"),
            n_components=n_components,
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
    legend: Union[bool, Mapping[str, Any]] = False,
    ax: Optional[Axes] = None,
    **kwargs: Any,
) -> Tuple[Figure, Axes]:

    category_values = __discrete_category_values(scdata.obs[obs])
    categories = set(scdata.obs[obs].unique()) & set(category_values)

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

    representation_mtx = __embedding_matrix(scdata, representation, n_components)
    fig, ax = __embedding_axes(ax, n_components, kwargs)

    kwargs["nan"] = kwargs["nan"] if "nan" in kwargs else {}
    __scatter_missing_values(
        scdata,
        obs,
        representation_mtx,
        ax,
        n_components,
        kwargs,
        gray,
    )

    color_values = cast(Sequence[object], colors)

    for _cluster, _color in zip(category_values, color_values):
        _color = normalize_color(_color)

        if _cluster not in categories:
            continue
        idx = np.where(scdata.obs[obs] == _cluster)[0]
        __scatter_values(
            ax,
            representation_mtx,
            idx,
            n_components,
            s=kwargs["s"] if "s" in kwargs else 2.0,
            facecolors=_color,
            edgecolors="none",
            alpha=__scatter_alpha(_color, kwargs),
            label=_cluster,
        )

    if legend is not False:
        box = ax.get_position()
        ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))

        handles, labels = ax.get_legend_handles_labels()
        apply_legend(ax, legend, handles=handles, labels=labels)

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

    representation_mtx = __embedding_matrix(scdata, representation, n_components)
    fig, ax = __embedding_axes(ax, n_components, kwargs)
    __scatter_missing_values(
        scdata,
        obs,
        representation_mtx,
        ax,
        n_components,
        kwargs,
        lightgray if n_components == 2 else gray,
    )

    if n_components == 2:
        sc = __scatter_values(
            ax,
            representation_mtx,
            slice(None),
            n_components,
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
        __scatter_values(
            ax,
            representation_mtx,
            slice(None),
            n_components,
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
