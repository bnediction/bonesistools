#!/usr/bin/env python

from __future__ import annotations

import importlib
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Mapping,
    Optional,
)

from ..._compat import Literal

if TYPE_CHECKING:
    from graphviz import Digraph  # pyright: ignore[reportMissingImports]

GraphvizProgram = Literal[
    "dot",
    "neato",
    "fdp",
    "sfdp",
    "circo",
    "twopi",
    "osage",
    "patchwork",
]


def _graphviz_attributes(
    attributes: Optional[Mapping[str, Any]],
) -> Mapping[str, str]:

    if attributes is None:
        return {}

    return {
        str(key): str(value) for key, value in attributes.items() if value is not None
    }


def _graph_attributes(
    graph_attr: Optional[Mapping[str, Any]],
    kwargs: Mapping[str, Any],
) -> Mapping[str, str]:

    attributes: Dict[str, Any] = {}
    if graph_attr is not None:
        attributes.update(graph_attr)
    attributes.update(kwargs)

    return _graphviz_attributes(attributes)


def _new_graphviz_digraph(
    program: GraphvizProgram,
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> "Digraph":

    try:
        graphviz = importlib.import_module("graphviz")

    except ImportError:
        raise ImportError(
            "to_graphviz() requires the 'graphviz' package. "
            "Install it with `pip install graphviz` or `bonesistools[graphviz]`."
        )

    graph = graphviz.Digraph(engine=program)

    graph_attributes = _graph_attributes(graph_attr, kwargs)
    node_attributes = _graphviz_attributes(node_attr)
    edge_attributes = _graphviz_attributes(edge_attr)
    if _uses_orthogonal_dot_edges(program, graph_attributes):
        edge_attributes = _externalized_edge_attributes(edge_attributes)

    if graph_attributes:
        graph.attr(**graph_attributes)
    if node_attributes:
        graph.attr("node", **node_attributes)
    if edge_attributes:
        graph.attr("edge", **edge_attributes)

    return graph


def _set_pydot_defaults(
    dot: Any,
    program: GraphvizProgram = "dot",
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> None:

    graph_attributes = _graph_attributes(graph_attr, kwargs)
    for key, value in graph_attributes.items():
        dot.set(key, value)

    node_attributes = _graphviz_attributes(node_attr)
    edge_attributes = _graphviz_attributes(edge_attr)
    if _uses_orthogonal_dot_edges(program, graph_attributes):
        edge_attributes = _externalized_edge_attributes(edge_attributes)

    if node_attributes:
        dot.set_node_defaults(**node_attributes)
    if edge_attributes:
        dot.set_edge_defaults(**edge_attributes)


def _networkx_to_graphviz(
    networkx_graph: Any,
    program: GraphvizProgram = "dot",
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> "Digraph":

    _externalize_orthogonal_edge_labels(
        networkx_graph,
        program=program,
        graph_attr=graph_attr,
        **kwargs,
    )

    graphviz_graph = _new_graphviz_digraph(
        program=program,
        graph_attr=graph_attr,
        node_attr=node_attr,
        edge_attr=edge_attr,
        **kwargs,
    )

    for node, data in networkx_graph.nodes(data=True):
        graphviz_graph.node(str(node), **_graphviz_attributes(data))

    for source, target, data in networkx_graph.edges(data=True):
        graphviz_graph.edge(
            str(source),
            str(target),
            **_graphviz_attributes(data),
        )

    return graphviz_graph


def _externalize_orthogonal_edge_labels(
    networkx_graph: Any,
    *,
    program: GraphvizProgram,
    graph_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> None:
    """Use external labels when dot routes edges orthogonally."""

    graph_attributes = _graph_attributes(graph_attr, kwargs)
    if not _uses_orthogonal_dot_edges(program, graph_attributes):
        return

    for _source, _target, data in networkx_graph.edges(data=True):
        if "label" not in data:
            continue

        label = data.pop("label")
        if data.get("xlabel") is None:
            data["xlabel"] = label


def _uses_orthogonal_dot_edges(
    program: GraphvizProgram,
    graph_attributes: Mapping[str, Any],
) -> bool:

    return (
        str(program).strip().lower() == "dot"
        and str(graph_attributes.get("splines", "")).strip().lower() == "ortho"
    )


def _externalized_edge_attributes(
    attributes: Mapping[str, str],
) -> Mapping[str, str]:

    externalized = dict(attributes)
    label = externalized.pop("label", None)
    if label is not None:
        externalized.setdefault("xlabel", label)
    return externalized
