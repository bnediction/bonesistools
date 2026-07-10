#!/usr/bin/env python

from __future__ import annotations

from importlib import import_module
from itertools import product
from typing import Any, FrozenSet, Mapping, Optional, Set, Tuple, cast

from boolean import (
    BooleanAlgebra,
    Expression,
)
from boolean.boolean import _FALSE, _TRUE

from ..._compat import Literal
from ..._validation import _as_literal
from ._hypercube import Hypercube
from ._typing import BooleanRule

EquivalenceMethod = Literal["simplify", "truth_table", "asp"]

_DNFLiteral = Tuple[str, int]
_DNFClause = FrozenSet[_DNFLiteral]
_BooleanCube = Tuple[Optional[int], ...]
Implicants = Tuple[Hypercube, ...]

BA = BooleanAlgebra()


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
    expressions on every assignment of their symbols. The `"asp"` method asks
    an ASP solver whether a counterexample assignment exists.

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

    Non-equivalent expressions return False with exact comparison:

    >>> expressions_equivalent(
    ...     ba.parse("B & C"),
    ...     ba.parse("B | C"),
    ...     method="asp",
    ...     ba=ba,
    ... )
    False

    Parameters
    ----------
    expr1: boolean.Expression
        First Boolean expression.
    expr2: boolean.Expression
        Second Boolean expression.
    method: {"simplify", "truth_table", "asp"} (default: "simplify")
        Equivalence strategy.

        - `"simplify"` compares both expressions after `boolean.py`
          simplification. This is fast, but not guaranteed to detect all
          logical equivalences.
        - `"truth_table"` exhaustively evaluates both expressions on all
          assignments of their symbols. This is an exact logical equivalence
          check, but has exponential complexity in the number of symbols.
        - `"asp"` uses ASP through `clingo` to search for a counterexample
          assignment. It is exact.
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
        If `method` is not one of `"simplify"`, `"truth_table"` or `"asp"`,
        or if an expression cannot be evaluated to a Boolean constant during
        exact comparison.
    """

    method = _as_literal(
        method,
        choices=("simplify", "truth_table", "asp"),
        name="method",
    )

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

    else:
        return _expressions_equivalent_with_asp(ba, expr1, expr2)


def dnf_implicants(
    rule: BooleanRule,
    value: Literal[0, 1] = 1,
    *,
    ba: BooleanAlgebra = BA,
) -> Implicants:
    """
    Decompose a Boolean rule into DNF implicants.

    DNF implicants are the conjunction clauses of a disjunctive normal form,
    represented as hypercubes. They are sufficient conditions for the rule to
    evaluate to `value`, but they are not necessarily prime implicants.

    Examples
    --------
    >>> dnf_implicants("(A & B) | (A & !B)")
    (Hypercube(A=1, B=0), Hypercube(A=1, B=1))

    Unlike `prime_implicants`, no absorption minimization is performed:

    >>> prime_implicants("(A & B) | (A & !B)")
    (Hypercube(A=1),)

    DNF implicants can also be computed for the false value:

    >>> dnf_implicants("A & !B", value=0)
    (Hypercube(A=0), Hypercube(B=1))

    Parameters
    ----------
    rule: str, bool, int or boolean.Expression
        Boolean rule to decompose.
    value: {0, 1} (default: 1)
        Target value forced by the returned DNF implicants.
    ba: BooleanAlgebra (default: module-level BooleanAlgebra)
        Boolean algebra used to parse string rules and evaluate expressions.

    Returns
    -------
    tuple of Hypercube
        DNF implicants forcing `rule` to `value`.

    Raises
    ------
    TypeError
        If `rule` has an unsupported type.
    ValueError
        If `value` is not 0 or 1, or if the DNF contains an invalid literal.
    """

    if value not in (0, 1):
        raise ValueError(
            "invalid argument value for 'value': "
            f"expected 0 or 1 but received {value!r}"
        )

    expr = _coerce_boolean_rule(ba, rule)
    if value == 0:
        expr = ba.NOT(expr)

    clauses = _dnf_clauses(ba, expr.literalize())

    return tuple(
        Hypercube(dict(sorted(clause)))
        for clause in sorted(clauses, key=_dnf_clause_sort_key)
    )


def prime_implicants(
    rule: BooleanRule,
    value: Literal[0, 1] = 1,
    *,
    backend: Literal["truth_table", "asp"] = "truth_table",
    ba: BooleanAlgebra = BA,
) -> Implicants:
    """
    Compute prime implicants of a Boolean rule.

    Prime implicants are minimal partial configurations that force a Boolean
    rule to a target value. Each returned `Hypercube` contains only fixed
    components; free components are omitted.

    Examples
    --------
    >>> prime_implicants("B | C")
    (Hypercube(B=1), Hypercube(C=1))

    >>> prime_implicants("B & !C")
    (Hypercube(B=1, C=0),)

    Prime implicants can also be computed for the false value:

    >>> prime_implicants("B & !C", value=0)
    (Hypercube(B=0), Hypercube(C=1))

    Parameters
    ----------
    rule: str, bool, int or boolean.Expression
        Boolean rule to analyze.
    value: {0, 1} (default: 1)
        Target value forced by the returned prime implicants.
    backend: {"truth_table", "asp"} (default: "truth_table")
        Backend used to enumerate minimal implicants. The `"truth_table"`
        backend uses an exact truth-table minimization. The `"asp"` backend
        uses ASP through `clingo` to enumerate minimal cubes.
    ba: BooleanAlgebra (default: module-level BooleanAlgebra)
        Boolean algebra used to parse string rules and evaluate expressions.

    Returns
    -------
    tuple of Hypercube
        Minimal partial configurations forcing `rule` to `value`.

    Raises
    ------
    TypeError
        If `rule` has an unsupported type.
    ValueError
        If `value` is not 0 or 1, or if `rule` cannot be evaluated to a
        Boolean constant for some assignment.
    """

    if value not in (0, 1):
        raise ValueError(
            "invalid argument value for 'value': "
            f"expected 0 or 1 but received {value!r}"
        )

    backend = _as_literal(
        backend,
        choices=("truth_table", "asp"),
        name="backend",
    )

    expr = _coerce_boolean_rule(ba, rule).simplify()
    symbols, minterms = _matching_minterms(ba, expr, target=bool(value))

    if not minterms:
        return ()

    if backend == "truth_table":
        prime_cubes = _compute_prime_cubes(tuple(minterms))

    else:
        prime_cubes = _compute_prime_cubes_with_clingo(
            tuple(minterms),
            n_variables=len(symbols),
        )

    prime_cubes = sorted(prime_cubes, key=_boolean_cube_sort_key)

    return tuple(
        Hypercube(
            {str(symbol): bit for symbol, bit in zip(symbols, cube) if bit is not None}
        )
        for cube in prime_cubes
    )


def _eval_as_bool(
    ba: BooleanAlgebra,
    expr: Expression,
    assignment: Mapping[Any, bool],
) -> bool:

    substitutions = {
        symbol: ba.TRUE if value else ba.FALSE for symbol, value in assignment.items()
    }

    value = expr.subs(substitutions).simplify()
    evaluated = _constant_expression_as_bool(ba, value)

    if evaluated is not None:
        return evaluated

    raise ValueError(f"expression did not evaluate to a Boolean constant: {value!r}")


def _constant_expression_as_bool(
    ba: BooleanAlgebra,
    expr: Expression,
) -> Optional[bool]:

    if isinstance(expr, _TRUE):
        return True

    if isinstance(expr, _FALSE):
        return False

    expr_any: Any = expr

    if isinstance(expr_any, ba.NOT):
        child = _constant_expression_as_bool(ba, expr_any.args[0].simplify())
        return None if child is None else not child

    if isinstance(expr_any, ba.AND):
        result = True
        for operand in expr_any.args:
            value = _constant_expression_as_bool(ba, operand.simplify())
            if value is None:
                return None
            result = result and value

        return result

    if isinstance(expr_any, ba.OR):
        result = False
        for operand in expr_any.args:
            value = _constant_expression_as_bool(ba, operand.simplify())
            if value is None:
                return None
            result = result or value

        return result

    return None


def _expressions_equivalent_with_asp(
    ba: BooleanAlgebra,
    expr1: Expression,
    expr2: Expression,
) -> bool:

    clingo = _import_clingo()
    expr1 = expr1.simplify()
    expr2 = expr2.simplify()
    symbols = tuple(sorted(expr1.symbols | expr2.symbols, key=str))
    symbol_ids = {str(symbol): index for index, symbol in enumerate(symbols)}

    root1, facts1, next_expression_id = _clingo_expression_tree_facts(
        ba,
        expr1,
        symbol_ids=symbol_ids,
        next_expression_id=0,
    )
    root2, facts2, _ = _clingo_expression_tree_facts(
        ba,
        expr2,
        symbol_ids=symbol_ids,
        next_expression_id=next_expression_id,
    )

    program = "\n".join(
        [
            *[f"symbol({index})." for index in range(len(symbols))],
            *facts1,
            *facts2,
            f"root(1, {root1}).",
            f"root(2, {root2}).",
            """
        1 { assign(S, 0); assign(S, 1) } 1 :- symbol(S).

        value(E, V) :- const_expr(E, V).
        value(E, V) :- symbol_expr(E, S), assign(S, V).
        value(E, 0) :- not_expr(E, C), value(C, 1).
        value(E, 1) :- not_expr(E, C), value(C, 0).
        value(E, 0) :- and_expr(E), arg(E, C), value(C, 0).
        value(E, 1) :- and_expr(E), not value(E, 0).
        value(E, 1) :- or_expr(E), arg(E, C), value(C, 1).
        value(E, 0) :- or_expr(E), not value(E, 1).

        different :- root(1, E1), root(2, E2), value(E1, V1), value(E2, V2),
                     V1 != V2.
        :- not different.
        #show assign/2.
        """,
        ]
    )

    control = clingo.Control(["--models=1", "--warn=none"])
    control.add("base", [], program)
    control.ground([("base", [])])

    with control.solve(yield_=True) as handle:
        return next(iter(handle), None) is None


def _clingo_expression_tree_facts(
    ba: BooleanAlgebra,
    expr: Expression,
    *,
    symbol_ids: Mapping[str, int],
    next_expression_id: int,
) -> Tuple[int, Tuple[str, ...], int]:

    expression_id = next_expression_id
    next_expression_id += 1

    if isinstance(expr, _TRUE):
        return expression_id, (f"const_expr({expression_id}, 1).",), next_expression_id

    if isinstance(expr, _FALSE):
        return expression_id, (f"const_expr({expression_id}, 0).",), next_expression_id

    expr_any: Any = expr

    if isinstance(expr_any, ba.NOT):
        child_id, child_facts, next_expression_id = _clingo_expression_tree_facts(
            ba,
            expr_any.args[0],
            symbol_ids=symbol_ids,
            next_expression_id=next_expression_id,
        )
        return (
            expression_id,
            (f"not_expr({expression_id}, {child_id}).", *child_facts),
            next_expression_id,
        )

    if isinstance(expr_any, ba.AND) or isinstance(expr_any, ba.OR):
        expression_facts = [
            f"{'and' if isinstance(expr_any, ba.AND) else 'or'}_expr({expression_id})."
        ]

        for operand in expr_any.args:
            child_id, child_facts, next_expression_id = _clingo_expression_tree_facts(
                ba,
                operand,
                symbol_ids=symbol_ids,
                next_expression_id=next_expression_id,
            )
            expression_facts.append(f"arg({expression_id}, {child_id}).")
            expression_facts.extend(child_facts)

        return expression_id, tuple(expression_facts), next_expression_id

    if expr_any.isliteral:
        return (
            expression_id,
            (f"symbol_expr({expression_id}, {symbol_ids[str(expr_any.obj)]}).",),
            next_expression_id,
        )

    raise ValueError(f"invalid Boolean expression for ASP equivalence: {expr!r}")


def _matching_minterms(
    ba: BooleanAlgebra,
    expr: Expression,
    target: bool,
) -> Tuple[Tuple[Expression, ...], Tuple[Tuple[int, ...], ...]]:

    symbols = tuple(sorted(expr.symbols, key=str))
    minterms = []

    for values in product([False, True], repeat=len(symbols)):
        assignment = dict(zip(symbols, values))
        if _eval_as_bool(ba, expr, assignment) == target:
            minterms.append(tuple(int(value) for value in values))

    return symbols, tuple(minterms)


def _coerce_boolean_rule(ba: BooleanAlgebra, rule: BooleanRule) -> Expression:

    if isinstance(rule, Expression):
        return ba.parse(str(rule))

    if isinstance(rule, bool):
        return cast(Expression, ba.TRUE if rule else ba.FALSE)

    if isinstance(rule, int):
        if rule in (0, 1):
            return cast(Expression, ba.TRUE if rule else ba.FALSE)

        raise ValueError(
            "invalid argument value for 'rule': "
            f"expected 0 or 1 for an integer rule but received {rule!r}"
        )

    if isinstance(rule, str):
        expr = rule.strip()

        if expr == "0":
            return cast(Expression, ba.FALSE)

        if expr == "1":
            return cast(Expression, ba.TRUE)

        return ba.parse(expr)

    raise TypeError(
        "unsupported argument type for 'rule': "
        f"expected str, bool, int or Expression but received {type(rule)}"
    )


def _compute_prime_cubes(
    minterms: Tuple[Tuple[int, ...], ...],
) -> FrozenSet[_BooleanCube]:

    current: Set[_BooleanCube] = set(minterms)
    prime_cubes: Set[_BooleanCube] = set()

    while current:
        combined: Set[_BooleanCube] = set()
        used: Set[_BooleanCube] = set()
        current_cubes = sorted(current, key=_boolean_cube_sort_key)

        for index, cube1 in enumerate(current_cubes):
            for cube2 in current_cubes[index + 1 :]:
                cube = _combine_boolean_cubes(cube1, cube2)

                if cube is None:
                    continue

                used.add(cube1)
                used.add(cube2)
                combined.add(cube)

        prime_cubes.update(cube for cube in current if cube not in used)
        current = combined

    return frozenset(
        cube
        for cube in prime_cubes
        if not any(
            cube != other and _boolean_cube_covers(other, cube) for other in prime_cubes
        )
    )


def _compute_prime_cubes_with_clingo(
    minterms: Tuple[Tuple[int, ...], ...],
    n_variables: int,
) -> FrozenSet[_BooleanCube]:

    clingo = _import_clingo()
    target_minterms = frozenset(minterms)
    opposite_minterms = tuple(
        minterm
        for minterm in product((0, 1), repeat=n_variables)
        if minterm not in target_minterms
    )

    program = "\n".join(
        [
            *_clingo_truth_table_facts(
                minterms,
                opposite_minterms,
                n_variables=n_variables,
            ),
            """
            { fixed(I, 0); fixed(I, 1) } :- var(I).
            :- fixed(I, 0), fixed(I, 1).

            mismatch(M) :- fixed(I, V), bit(M, I, W), V != W.
            covered_target(M) :- target(M), not mismatch(M).
            covered_opposite(M) :- opposite(M), not mismatch(M).

            :- not 1 { covered_target(M) : target(M) }.
            :- covered_opposite(M).

            other_mismatch(M, I) :-
                var(I), fixed(J, V), bit(M, J, W), I != J, V != W.
            drop_witness(I) :-
                fixed(I, V), opposite(M), bit(M, I, W), V != W,
                not other_mismatch(M, I).
            :- fixed(I, V), not drop_witness(I).

            #show fixed/2.
            """,
        ]
    )

    control = clingo.Control(
        ["--models=0", "--opt-mode=opt", "--opt-strategy=usc", "--warn=none"]
    )
    control.add("base", [], program)
    control.ground([("base", [])])

    cubes: Set[_BooleanCube] = set()
    with control.solve(yield_=True) as handle:
        for model in handle:
            cube = [None] * n_variables

            for atom in model.symbols(shown=True):
                variable, fixed_value = atom.arguments
                cube[variable.number] = fixed_value.number

            cubes.add(tuple(cube))

    return frozenset(cubes)


def _import_clingo() -> Any:

    try:
        return import_module("clingo")

    except ImportError as error:
        raise ImportError(
            "The ASP backend requires `clingo` to be installed."
        ) from error


def _clingo_truth_table_facts(
    minterms: Tuple[Tuple[int, ...], ...],
    opposite_minterms: Tuple[Tuple[int, ...], ...],
    n_variables: int,
) -> Tuple[str, ...]:

    facts = [f"var({index})." for index in range(n_variables)]

    for minterm_id, minterm in enumerate(minterms):
        facts.append(f"target({minterm_id}).")
        facts.extend(_clingo_minterm_facts(minterm_id, minterm))

    offset = len(minterms)
    for index, minterm in enumerate(opposite_minterms, start=offset):
        facts.append(f"opposite({index}).")
        facts.extend(_clingo_minterm_facts(index, minterm))

    return tuple(facts)


def _clingo_minterm_facts(minterm_id: int, minterm: Tuple[int, ...]) -> Tuple[str, ...]:

    return tuple(
        f"bit({minterm_id}, {variable}, {value})."
        for variable, value in enumerate(minterm)
    )


def _combine_boolean_cubes(
    cube1: _BooleanCube,
    cube2: _BooleanCube,
) -> Optional[_BooleanCube]:

    differing_index = None

    for index, (value1, value2) in enumerate(zip(cube1, cube2)):
        if value1 == value2:
            continue

        if value1 is None or value2 is None or differing_index is not None:
            return None

        differing_index = index

    if differing_index is None:
        return None

    combined = list(cube1)
    combined[differing_index] = None

    return tuple(combined)


def _boolean_cube_covers(cube1: _BooleanCube, cube2: _BooleanCube) -> bool:

    return all(
        value1 is None or value1 == value2 for value1, value2 in zip(cube1, cube2)
    )


def _boolean_cube_sort_key(cube: _BooleanCube) -> Tuple[int, Tuple[int, ...]]:

    return (
        sum(value is not None for value in cube),
        tuple(2 if value is None else value for value in cube),
    )


def _dnf_clauses(ba: BooleanAlgebra, expr: Expression) -> Tuple[_DNFClause, ...]:

    if isinstance(expr, _TRUE):
        return (frozenset(),)

    if isinstance(expr, _FALSE):
        return ()

    expr_any: Any = expr

    if isinstance(expr_any, ba.OR):
        return tuple(
            clause for operand in expr_any.args for clause in _dnf_clauses(ba, operand)
        )

    if isinstance(expr_any, ba.AND):
        clauses = (frozenset(),)

        for operand in expr_any.args:
            operand_clauses = _dnf_clauses(ba, operand)
            clauses = tuple(
                combined
                for clause in clauses
                for operand_clause in operand_clauses
                for combined in [_combine_dnf_clauses(clause, operand_clause)]
                if combined is not None
            )

        return clauses

    return (_literal_to_dnf_clause(ba, expr),)


def _combine_dnf_clauses(
    clause1: _DNFClause,
    clause2: _DNFClause,
) -> Optional[_DNFClause]:

    literals = dict(clause1)

    for component, value in clause2:
        if component in literals and literals[component] != value:
            return None

        literals[component] = value

    return frozenset(literals.items())


def _literal_to_dnf_clause(ba: BooleanAlgebra, literal: Expression) -> _DNFClause:

    literal_any: Any = literal

    if isinstance(literal_any, ba.NOT):
        operand = literal_any.args[0]

        if not operand.isliteral:
            raise ValueError(f"invalid DNF literal: {literal_any!r}")

        return frozenset({(str(operand.obj), 0)})

    if literal_any.isliteral:
        return frozenset({(str(literal_any.obj), 1)})

    raise ValueError(f"invalid DNF literal: {literal_any!r}")


def _dnf_clause_sort_key(clause: _DNFClause) -> Tuple[int, Tuple[_DNFLiteral, ...]]:

    return len(clause), tuple(sorted(clause))
