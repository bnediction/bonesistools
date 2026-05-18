#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Callable, Iterable, Literal, Tuple

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

    Parameters
    ----------
    expr1:
        First Boolean expression.
    expr2:
        Second Boolean expression.
    method:
        Equivalence strategy.

        - `"simplify"` compares both expressions after `boolean.py`
          simplification. This is fast, but not guaranteed to detect all
          logical equivalences.
        - `"truth_table"` exhaustively evaluates both expressions on all
          assignments of their symbols. This is an exact logical equivalence
          check, but has exponential complexity in the number of symbols.

    Returns
    -------
    bool
        Whether the two expressions are considered equivalent according to
        the selected method.
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
        f"unsupported equivalence method: {method}. "
        "Expected 'simplify' or 'truth_table'."
    )


def dnf_to_structure(
    ba: BooleanAlgebra,
    expr: Expression,
    container: Callable[[Iterable[Any]], Any] = frozenset,
    sort: bool = False,
) -> DNFStructure:
    """
    Convert a Boolean expression in disjunctive normal form into a nested
    Python structure.

    The returned structure encodes a DNF as a collection of clauses, where each
    clause is a collection of signed literals. Each literal is represented as
    a pair `(symbol, sign)`, where `sign` is True for a positive literal and
    False for a negated literal.

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
    DNFStructure
        Return True or False for constant expressions. Otherwise, return a
        nested structure representing the DNF.

    Examples
    --------
    The expression `(A & ~B) | C` is converted into:
    `{{('A', True), ('B', False)}, {('C', True)}}`
    """

    def make_literal(literal):
        if isinstance(literal, ba.NOT):
            return (literal.args[0].obj, False)
        return (literal.obj, True)

    def make_clause(clause):
        if isinstance(clause, ba.AND):
            literals = clause.args
        else:
            literals = [clause]

        literals = map(make_literal, literals)
        return container(sorted(literals) if sort else literals)

    if isinstance(expr, _TRUE):
        return True

    if isinstance(expr, _FALSE):
        return False

    if isinstance(expr, ba.OR):
        clauses = expr.args
    else:
        clauses = [expr]

    clauses = map(make_clause, clauses)
    return container(sorted(clauses) if sort else clauses)
