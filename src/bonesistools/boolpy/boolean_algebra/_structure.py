#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Callable, Iterable, Tuple, Union

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from itertools import product

from boolean import (
    BooleanAlgebra,
    Expression,
)
from boolean.boolean import _TRUE, _FALSE

EquivalenceMethod = Literal["simplify", "truth_table"]

RegulatoryLiteral = Tuple[str, bool]
Clause = Iterable[RegulatoryLiteral]
DNFStructure = Iterable[Clause]
DNFValue = Union[bool, DNFStructure]

BA = BooleanAlgebra()


def _eval_as_bool(
    ba: BooleanAlgebra,
    expr: Expression,
    assignment: dict,
) -> bool:
    assignment = {
        symbol: ba.TRUE if value else ba.FALSE for symbol, value in assignment.items()
    }

    value = expr.subs(assignment).simplify()

    if isinstance(value, _TRUE):
        return True

    if isinstance(value, _FALSE):
        return False

    raise ValueError(f"expression did not evaluate to a Boolean constant: {value!r}")


def expressions_equivalent(
    expr1: Expression,
    expr2: Expression,
    method: EquivalenceMethod = "simplify",
    ba: BooleanAlgebra = BA,
) -> bool:
    """
    Test whether two Boolean expressions are logically equivalent.

    The `"simplify"` method is a fast structural comparison after
    `boolean.py` simplification. It can miss logical equivalences that require
    a full semantic check. The `"truth_table"` method evaluates both
    expressions on every assignment of their symbols and is therefore exact,
    but exponential in the number of symbols.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> expr1 = ba.parse("(B & C) | (~B & D) | (C & D)")
    >>> expr2 = ba.parse("(B & C) | (~B & D)")

    The simplification method compares simplified expression structure and may
    miss this equivalence:

    >>> expressions_equivalent(expr1, expr2, method="simplify", ba=ba)
    False

    Exhaustive truth-table comparison detects the logical equivalence:

    >>> expressions_equivalent(expr1, expr2, method="truth_table", ba=ba)
    True

    Non-equivalent expressions return False with truth-table comparison:

    >>> expressions_equivalent(
    ...     ba.parse("B & C"),
    ...     ba.parse("B | C"),
    ...     method="truth_table",
    ...     ba=ba,
    ... )
    False

    Parameters
    ----------
    expr1: boolean.Expression
        First Boolean expression.
    expr2: boolean.Expression
        Second Boolean expression.
    method: {"simplify", "truth_table"} (default: "simplify")
        Equivalence strategy.

        - `"simplify"` compares both expressions after `boolean.py`
          simplification. This is fast, but not guaranteed to detect all
          logical equivalences.
        - `"truth_table"` exhaustively evaluates both expressions on all
          assignments of their symbols. This is an exact logical equivalence
          check, but has exponential complexity in the number of symbols.
    ba: boolean.BooleanAlgebra (default: module-level BooleanAlgebra)
        Boolean algebra used to evaluate truth-table assignments. It should be
        compatible with the algebra used to build `expr1` and `expr2`.

    Returns
    -------
    bool
        Whether the two expressions are considered equivalent according to
        the selected method.

    Raises
    ------
    ValueError
        If `method` is not one of `"simplify"` or `"truth_table"`, or if an
        expression cannot be evaluated to a Boolean constant during
        truth-table comparison.
    """

    if method == "simplify":
        return expr1.simplify() == expr2.simplify()

    if method == "truth_table":
        expr1 = expr1.simplify()
        expr2 = expr2.simplify()

        symbols = sorted(expr1.symbols | expr2.symbols, key=str)

        for values in product([False, True], repeat=len(symbols)):
            assignment = dict(zip(symbols, values))

            value1 = _eval_as_bool(ba, expr1, assignment)
            value2 = _eval_as_bool(ba, expr2, assignment)

            if value1 != value2:
                return False

        return True

    raise ValueError(
        f"invalid argument value for 'method': "
        f"expected 'simplify' or 'truth_table' but received {method!r}"
    )


def _literal_to_pair(ba: BooleanAlgebra, literal: Expression) -> RegulatoryLiteral:
    literal_any: Any = literal

    if isinstance(literal_any, ba.NOT):
        operand = literal_any.args[0]

        if not operand.isliteral:
            raise ValueError(f"invalid DNF literal: {literal_any!r}")

        return str(operand.obj), False

    if literal_any.isliteral:
        return str(literal_any.obj), True

    raise ValueError(f"invalid DNF literal: {literal_any!r}")


def _clause_to_structure(
    ba: BooleanAlgebra,
    clause: Expression,
    container: Callable[[Iterable[Any]], Any],
    sort: bool,
):
    if isinstance(clause, ba.AND):
        literals = clause.args
    else:
        literals = (clause,)

    literals = [_literal_to_pair(ba, literal) for literal in literals]

    if sort:
        literals = sorted(literals)

    return container(literals)


def dnf_to_structure(
    ba: BooleanAlgebra,
    expr: Expression,
    container: Callable[[Iterable[Any]], Any] = frozenset,
    sort: bool = False,
) -> DNFValue:
    """
    Convert a Boolean expression in disjunctive normal form into a nested
    Python structure.

    The returned structure encodes a DNF as a collection of clauses, where each
    clause is a collection of signed literals. Each literal is represented as
    a pair `(symbol, sign)`, where `sign` is True for a positive literal and
    False for a negated literal.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> expr = ba.parse("(A & ~B) | C")
    >>> dnf_to_structure(ba, expr, container=list, sort=True)
    [[('A', True), ('B', False)], [('C', True)]]

    Parameters
    ----------
    ba: BooleanAlgebra
        Boolean algebra used to build the expression.
    expr: Expression
        Boolean expression assumed to be in disjunctive normal form.
    container: Callable[[Iterable[Any]], Any] (default: frozenset)
        Container used to store clauses and literals.
    sort: bool (default: False)
        Sort literals and clauses before storing them.

    Returns
    -------
    DNFValue
        Return True or False for constant expressions. Otherwise, return a
        nested structure representing the DNF.
    """

    if isinstance(expr, _TRUE):
        return True

    if isinstance(expr, _FALSE):
        return False

    if isinstance(expr, ba.OR):
        clauses = expr.args
    else:
        clauses = (expr,)

    clauses = [
        _clause_to_structure(ba, clause, container=container, sort=sort)
        for clause in clauses
    ]

    if sort:
        clauses = sorted(clauses, key=repr)

    return container(clauses)
