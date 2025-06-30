#!/usr/bin/env python

import importlib

import networkx as nx

try:
    _mpbn_is_available = importlib.util.find_spec("mpbn") is not None
    _pydot_is_available = importlib.util.find_spec("pydot") is not None
except:
    _mpbn_is_available = importlib.find_loader("mpbn") is not None
    _pydot_is_available = importlib.find_loader("pydot") is not None

if _mpbn_is_available:
    from colomoto.minibn import BooleanNetwork # type: ignore
else:
    BooleanNetwork = type(NotImplemented)

if _pydot_is_available:
    from pydot import Dot # type: ignore
else:
    Dot = type(NotImplemented)

def bn_to_pydot(
    bn: BooleanNetwork, # type: ignore
    **kwargs
) -> Dot: # type: ignore
    """
    Convert Boolean network into pydot graph.

    Parameters
    ----------
    bn: colomoto.minibn.BooleanNetwork
        Boolean network.
    **kwargs: Mapping[str, Any]
        Keyword-arguments passed to update pydot graph.
    
    Returns
    -------
    Return pydot graph from Boolean network-based influence graph.
    """

    import boolean # type: ignore
    
    edges = set()
    ig = nx.MultiDiGraph()
    for target, f in bn.items():
        ig.add_node(target)
        for literal in f.simplify().literalize().get_literals():
            if isinstance(literal, boolean.NOT):
                source = literal.args[0].obj
                sign = -1
            else:
                source = literal.obj
                sign = 1
            edges.add((source,target,sign))
    for (source, target, sign) in edges:
        ig.add_edge(
            source,
            target,
            sign=sign,
            color="green4" if sign==1 else "red2",
            arrowhead="normal" if sign==1 else "tee",
            penwidth=2
        )
    dot = nx.drawing.nx_pydot.to_pydot(ig)
    for k, v in kwargs.items():
        dot.set(k, v)
    return dot
