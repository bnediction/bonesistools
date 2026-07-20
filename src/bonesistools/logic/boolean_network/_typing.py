#!/usr/bin/env python

"""
Typing aliases and runtime validators for BooleanNetworkLike objects.
"""

from collections.abc import Mapping as _MappingABC
from typing import TYPE_CHECKING as _TYPE_CHECKING
from typing import Any as _Any
from typing import Mapping as _Mapping
from typing import Union as _Union

from typing_extensions import TypeAlias as _TypeAlias
from typing_extensions import TypeGuard as _TypeGuard

from ..._compat import Literal as _Literal
from ..boolean_algebra._typing import BooleanRule, is_boolean_rule_like

BooleanNetworkLike: _TypeAlias = _Mapping[str, BooleanRule]
BooleanNetworkMetric: _TypeAlias = _Literal["equivalence", "hamming"]

if _TYPE_CHECKING:
    from dd.autoref import BDD as _AutoRefBDD
    from dd.autoref import Function as _AutoRefBDDNode
    from dd.cudd import BDD as _CuddBDD
    from dd.cudd import Function as _CuddBDDNode

    _BDDManager: _TypeAlias = _Union[_CuddBDD, _AutoRefBDD]
    _BDDNode: _TypeAlias = _Union[_CuddBDDNode, _AutoRefBDDNode]
    _BDDConfigurationSetNode: _TypeAlias = _BDDNode
    _BDDTransitionRelationNode: _TypeAlias = _BDDNode


def is_boolean_network_like(obj: _Any) -> _TypeGuard[BooleanNetworkLike]:
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

    if not isinstance(obj, _MappingABC):
        return False

    try:
        return all(
            isinstance(node, str) and is_boolean_rule_like(rule)
            for node, rule in obj.items()
        )
    except Exception:
        return False
