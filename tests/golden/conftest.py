#!/usr/bin/env python

import os
from pathlib import Path

import pytest

_GOLDEN_ENV = "BONESISTOOLS_RUN_GOLDEN"
_GOLDEN_DIR = Path(__file__).parent


def pytest_collection_modifyitems(items):

    if os.environ.get(_GOLDEN_ENV) == "1":
        return

    skip = pytest.mark.skip(reason=f"requires golden CI mode ({_GOLDEN_ENV}=1)")
    for item in items:
        try:
            Path(str(item.path)).relative_to(_GOLDEN_DIR)
        except ValueError:
            continue

        else:
            item.add_marker(skip)
