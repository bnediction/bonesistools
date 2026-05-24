#!/usr/bin/env python

"""
Typing aliases and runtime validators for BooleanNetworkLike objects.
"""

from collections.abc import Mapping as MappingABC
from typing import TYPE_CHECKING, Any, Mapping

from typing_extensions import TypeAlias, TypeGuard

from ..boolean_algebra import BooleanRule, is_boolean_rule_like

if TYPE_CHECKING:
    BooleanNetworkLike: TypeAlias = Mapping[str, BooleanRule]
else:
    BooleanNetworkLike = MappingABC


def is_boolean_network_like(obj: Any) -> TypeGuard[BooleanNetworkLike]:
    """
    Test whether an object behaves as a BooleanNetworkLike mapping.

    Parameters
    ----------
    obj: Any
        Object to test.

    Returns
    -------
    bool
        True if `obj` exposes valid string components and Boolean-rule-like
        values.
    """

    if not isinstance(obj, MappingABC):
        return False

    try:
        return all(
            isinstance(node, str) and is_boolean_rule_like(rule)
            for node, rule in obj.items()
        )
    except Exception:
        return False
