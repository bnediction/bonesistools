#!/usr/bin/env python

import importlib
from typing import (
    Any,
    Dict,
    Mapping,
    Optional,
)


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
    program: str,
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> Any:

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

    if graph_attributes:
        graph.attr(**graph_attributes)
    if node_attributes:
        graph.attr("node", **node_attributes)
    if edge_attributes:
        graph.attr("edge", **edge_attributes)

    return graph


def _set_pydot_defaults(
    dot: Any,
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> None:

    for key, value in _graph_attributes(graph_attr, kwargs).items():
        dot.set(key, value)

    node_attributes = _graphviz_attributes(node_attr)
    edge_attributes = _graphviz_attributes(edge_attr)

    if node_attributes:
        dot.set_node_defaults(**node_attributes)
    if edge_attributes:
        dot.set_edge_defaults(**edge_attributes)


def _networkx_to_graphviz(
    networkx_graph: Any,
    program: str = "dot",
    graph_attr: Optional[Mapping[str, Any]] = None,
    node_attr: Optional[Mapping[str, Any]] = None,
    edge_attr: Optional[Mapping[str, Any]] = None,
    **kwargs: Any,
) -> Any:

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
