#!/usr/bin/env python

import numpy as np
import pytest

import bonesistools as bt

from ._boolean_workflow import (
    DEATH_RECEPTOR_CELL_FATE_PATH,
    SYNTHETIC_30_NODE_DYNAMICS_PATH,
    available_expected_names,
    load_expected,
    run_boolean_workflow,
)


@pytest.fixture(scope="module")
def golden_boolean_outputs():

    return run_boolean_workflow()


@pytest.mark.parametrize("name", available_expected_names())
def test_golden_boolean_workflow_matches_expected_outputs(
    golden_boolean_outputs,
    name,
):
    result = golden_boolean_outputs[name]
    expected = load_expected(name)

    assert sorted(result) == sorted(expected)

    for key in sorted(expected):
        np.testing.assert_array_equal(result[key], expected[key])


@pytest.mark.parametrize(
    "path",
    [DEATH_RECEPTOR_CELL_FATE_PATH, SYNTHETIC_30_NODE_DYNAMICS_PATH],
)
def test_golden_imported_influence_graph_matches_boolean_network(path):
    model = bt.logic.io.read_zginml(path)
    imported = model.influence_graph
    expected = model.boolean_network.to_influence_graph()

    assert set(imported) == set(expected)
    assert {
        (source, target, attributes["sign"])
        for source, target, attributes in imported.edges(data=True)
    } == {
        (source, target, attributes["sign"])
        for source, target, attributes in expected.edges(data=True)
    }
