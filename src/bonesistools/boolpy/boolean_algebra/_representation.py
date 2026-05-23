#!/usr/bin/env python

"""
String representation helpers for Boolean algebra expressions.
"""

from boolean import Expression
from boolean.boolean import _FALSE, _TRUE


def rule_to_string(rule: Expression) -> str:
    """
    Convert a Boolean rule into a readable string representation.

    The function keeps the syntax used by the `boolean.py` package while adding
    spaces around binary Boolean operators. Negations are preserved as `~`.
    Boolean constants are represented as `0` and `1`, which matches the Boolean
    network string representation used elsewhere in boolpy.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> rule_to_string(ba.parse("A&B"))
    'A & B'

    Negations are kept compact:

    >>> rule_to_string(ba.parse("A|~B"))
    'A | ~B'

    Parentheses are preserved while binary operators are spaced:

    >>> rule_to_string(ba.parse("(A&B)|~C"))
    '(A & B) | ~C'

    Boolean constants are converted to numeric strings:

    >>> rule_to_string(ba.TRUE)
    '1'
    >>> rule_to_string(ba.FALSE)
    '0'

    Parameters
    ----------
    rule: boolean.Expression
        Boolean rule to convert.

    Returns
    -------
    str
        Readable string representation of `rule`.
    """

    if isinstance(rule, _TRUE):
        return "1"

    if isinstance(rule, _FALSE):
        return "0"

    return str(rule).replace("&", " & ").replace("|", " | ")
