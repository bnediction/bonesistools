#!/usr/bin/env python

import sys
from types import SimpleNamespace

import pytest


class FakeGraphvizDigraph:
    def __init__(self, engine="dot"):
        self.engine = engine
        self.graph_attr = {}
        self.node_attr = {}
        self.edge_attr = {}
        self.nodes = []
        self.edges = []

    def attr(self, kind=None, **attrs):
        if kind == "node":
            self.node_attr.update(attrs)
            return
        if kind == "edge":
            self.edge_attr.update(attrs)
            return
        self.graph_attr.update(attrs)

    def node(self, name, **attrs):
        self.nodes.append((name, attrs))

    def edge(self, tail_name, head_name, **attrs):
        self.edges.append((tail_name, head_name, attrs))


@pytest.fixture
def fake_graphviz(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "graphviz",
        SimpleNamespace(Digraph=FakeGraphvizDigraph),
    )

    return FakeGraphvizDigraph
