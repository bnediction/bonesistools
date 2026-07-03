#!/usr/bin/env python

import os

import pytest

_GOLDEN_ENV = "BONESISTOOLS_RUN_GOLDEN"


def pytest_collection_modifyitems(items):

    if os.environ.get(_GOLDEN_ENV) == "1":
        return

    skip = pytest.mark.skip(reason=f"requires golden CI mode ({_GOLDEN_ENV}=1)")
    for item in items:
        item.add_marker(skip)
