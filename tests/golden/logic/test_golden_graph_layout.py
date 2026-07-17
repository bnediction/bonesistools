#!/usr/bin/env python

import pytest

from ._graph_layout import EXPECTED_PATH, graphviz_input


def test_golden_aggregated_influence_graph_graphviz_input():
    pytest.importorskip("pydot")

    assert graphviz_input() == EXPECTED_PATH.read_text()
