#!/usr/bin/env python

from typing import Any, ItemsView

try:
    from typing import Protocol, runtime_checkable
except ImportError:
    from typing_extensions import Protocol, runtime_checkable

from ..boolean_algebra import BooleanRule, is_boolean_rule_like


@runtime_checkable
class BooleanNetworkLike(Protocol):
    """
    Protocol for Boolean network-like objects.

    A Boolean network-like object is expected to behave as a mapping from
    component names to Boolean rules.
    """

    def items(self) -> ItemsView[str, BooleanRule]: ...


def is_boolean_network_like(obj: Any) -> bool:
    if not isinstance(obj, BooleanNetworkLike):
        return False

    try:
        return all(
            isinstance(node, str) and is_boolean_rule_like(rule)
            for node, rule in obj.items()
        )
    except Exception:
        return False
