#!/usr/bin/env python

"""
Dataset metadata loading.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, cast

_METADATA_PATH = Path(__file__).with_name("datasets.json")

with _METADATA_PATH.open(encoding="utf-8") as handle:
    _metadata = json.load(handle)

if not isinstance(_metadata, dict):
    raise ValueError("dataset metadata must be a JSON object")

_DATASETS = cast(Dict[str, Dict[str, Any]], _metadata)
