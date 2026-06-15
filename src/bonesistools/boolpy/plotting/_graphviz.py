#!/usr/bin/env python

import importlib
from typing import (
    Any,
    Mapping,
)


def _graphviz_attributes(attributes: Mapping[str, Any]) -> Mapping[str, str]:

    return {
        str(key): str(value) for key, value in attributes.items() if value is not None
    }


def _new_graphviz_digraph(program: str, **kwargs: Any) -> Any:

    try:
        graphviz = importlib.import_module("graphviz")

    except ImportError:
        raise ImportError(
            "to_graphviz() requires the 'graphviz' package. "
            "Install it with `pip install graphviz` or `bonesistools[graphviz]`."
        )

    graph = graphviz.Digraph(engine=program)

    if kwargs:
        graph.attr(**_graphviz_attributes(kwargs))

    return graph


def _networkx_to_graphviz(
    networkx_graph: Any,
    program: str = "dot",
    **kwargs: Any,
) -> Any:

    graphviz_graph = _new_graphviz_digraph(program=program, **kwargs)

    for node, data in networkx_graph.nodes(data=True):
        graphviz_graph.node(str(node), **_graphviz_attributes(data))

    for source, target, data in networkx_graph.edges(data=True):
        graphviz_graph.edge(
            str(source),
            str(target),
            **_graphviz_attributes(data),
        )

    return graphviz_graph
