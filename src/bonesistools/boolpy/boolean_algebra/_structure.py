#!/usr/bin/env python

from __future__ import annotations

from typing import Union, Any, Callable, Iterable, Tuple

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


def _literal_to_pair(ba: BooleanAlgebra, literal: Expression) -> RegulatoryLiteral:
    if isinstance(literal, ba.NOT):
        if not literal.args[0].isliteral:
            raise ValueError(f"invalid DNF literal: {literal!r}")

        return literal.args[0].obj, False

    if literal.isliteral:
        return literal.obj, True

    raise ValueError(f"invalid DNF literal: {literal!r}")


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

    Examples
    --------
    The expression `(A & ~B) | C` is converted into:
    `{{('A', True), ('B', False)}, {('C', True)}}`
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
