#!/usr/bin/env python

from typing import Any, cast

import pytest

from bonesistools._validation import _as_memory_size


def test_as_memory_size_parses_decimal_and_binary_units():
    assert _as_memory_size("1KB", "max_memory") == 1_000
    assert _as_memory_size("2MB", "max_memory") == 2_000_000
    assert _as_memory_size("3GB", "max_memory") == 3_000_000_000
    assert _as_memory_size("4TB", "max_memory") == 4_000_000_000_000

    assert _as_memory_size("1KiB", "max_memory") == 1024
    assert _as_memory_size("2MiB", "max_memory") == 2 * 1024**2
    assert _as_memory_size("3GiB", "max_memory") == 3 * 1024**3
    assert _as_memory_size("4TiB", "max_memory") == 4 * 1024**4


def test_as_memory_size_accepts_bytes_and_decimal_values():
    assert _as_memory_size(500_000_000, "max_memory") == 500_000_000
    assert _as_memory_size("0.5GB", "max_memory") == 500_000_000
    assert _as_memory_size("1.5MiB", "max_memory") == int(1.5 * 1024**2)


def test_as_memory_size_rejects_invalid_values():
    with pytest.raises(TypeError, match="unsupported argument type for 'max_memory'"):
        _as_memory_size(cast(Any, 1.0), "max_memory")

    with pytest.raises(ValueError, match="invalid argument value for 'max_memory'"):
        _as_memory_size(0, "max_memory")

    with pytest.raises(ValueError, match="invalid argument value for 'max_memory'"):
        _as_memory_size("1XB", "max_memory")

    with pytest.raises(ValueError, match="invalid argument value for 'max_memory'"):
        _as_memory_size("100", "max_memory")
