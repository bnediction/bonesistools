#!/usr/bin/env python

import os
from typing import Any, Tuple, cast

import pytest

import bonesistools as bt

pytestmark = pytest.mark.skipif(
    os.environ.get("BONESISTOOLS_RUN_REPRODUCIBILITY") != "1",
    reason="requires reproducibility CI mode",
)


def _canonical_attractors(
    attractors,
) -> Tuple[Tuple[Tuple[Tuple[str, int], ...], ...], ...]:

    canonical = []
    for attractor in attractors:
        states = tuple(
            sorted(tuple(sorted(state.items())) for state in attractor.enumerate())
        )
        canonical.append(states)

    return tuple(sorted(canonical))


@pytest.mark.parametrize("update", ["synchronous", "asynchronous", "general"])
def test_reachable_attractors_match_between_explicit_and_bdd_backends(update):
    pytest.importorskip("dd.autoref")

    bn = bt.bpy.bn.BooleanNetwork(
        {
            "A": "B",
            "B": "A",
            "C": "~C",
            "D": "(A & ~C) | D",
        }
    )
    initial_state = {
        "A": 0,
        "B": "*",
        "C": "*",
        "D": 0,
    }

    explicit = bn.reachable_attractors(
        initial_state,
        update=cast(Any, update),
        backend="explicit",
    )
    bdd = bn.reachable_attractors(
        initial_state,
        update=cast(Any, update),
        backend="bdd",
    )

    assert _canonical_attractors(bdd) == _canonical_attractors(explicit)
