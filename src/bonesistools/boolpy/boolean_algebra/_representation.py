#!/usr/bin/env python

from boolean import Expression
from boolean.boolean import _FALSE, _TRUE

def rule_to_string(rule: Expression) -> str:
    """
    Convert a Boolean rule into a readable string representation.

    Boolean constants are represented as `0` and `1`.
    """

    if isinstance(rule, _TRUE):
        return "1"

    if isinstance(rule, _FALSE):
        return "0"

    return (
        str(rule)
        .replace("~", "!")
        .replace("&", " & ")
        .replace("|", " | ")
    )