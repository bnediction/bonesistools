#!/usr/bin/env python

import numpy as np
import pytest

from ._boolean_workflow import (
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
