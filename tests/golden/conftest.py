#!/usr/bin/env python

import os
from pathlib import Path
from typing import Any, cast

import pytest
from typing_extensions import Literal

_GOLDEN_ENV = "BONESISTOOLS_RUN_GOLDEN"
_GOLDEN_DIR = Path(__file__).parent
GoldenMode = Literal["strict", "portable"]


def pytest_addoption(parser: Any) -> None:

    parser.getgroup("golden").addoption(
        "--golden-mode",
        choices=("strict", "portable"),
        default="strict",
        help="golden comparison contract (default: strict)",
    )


@pytest.fixture(scope="session")
def golden_mode(pytestconfig: Any) -> GoldenMode:

    return cast(GoldenMode, pytestconfig.getoption("--golden-mode"))


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
