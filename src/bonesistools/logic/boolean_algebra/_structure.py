#!/usr/bin/env python

from __future__ import annotations

from itertools import product
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    FrozenSet,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    cast,
)

import clingo
from boolean import (
    BooleanAlgebra,
    Expression,
)
from boolean.boolean import _FALSE, _TRUE

from ..._compat import Literal
from ._hypercube import Hypercube
from ._typing import BooleanRule

if TYPE_CHECKING:
    from ._robdd import ROBDD

_DNFLiteral = Tuple[str, int]
_DNFClause = FrozenSet[_DNFLiteral]
_BooleanCube = Tuple[Optional[int], ...]
Implicants = Tuple[Hypercube, ...]

BA = BooleanAlgebra()


def equivalence(
    expr1: Expression,
    expr2: Expression,
    ba: BooleanAlgebra = BA,
) -> bool:
    """
    Test whether two Boolean expressions are logically equivalent.

    Examples
    --------
    >>> from boolean import BooleanAlgebra
    >>> ba = BooleanAlgebra()
    >>> expr1 = ba.parse("(B & C) | (~B & D) | (C & D)")
    >>> expr2 = ba.parse("(B & C) | (~B & D)")

    >>> equivalence(expr1, expr2, ba=ba)
    True

    Non-equivalent expressions return False:

    >>> equivalence(
    ...     ba.parse("B & C"),
    ...     ba.parse("B | C"),
    ...     ba=ba,
    ... )
    False

    Parameters
    ----------
    expr1: boolean.Expression
        First Boolean expression.
    expr2: boolean.Expression
        Second Boolean expression.
    ba: boolean.BooleanAlgebra (default: module-level BooleanAlgebra)
        Boolean algebra compatible with the one used to build `expr1` and
        `expr2`.

    Returns
    -------
    bool
        Whether the two expressions are logically equivalent.
    """

    if expr1 == expr2:
        return True

    expr1 = expr1.simplify()
    expr2 = expr2.simplify()
    if expr1 == expr2:
        return True

    symbols = tuple(sorted(expr1.symbols | expr2.symbols, key=str))

    if len(symbols) <= 15:
        return _equivalence_with_bitsets(
            ba,
            expr1,
            expr2,
            symbols=symbols,
        )

    return _equivalence_with_asp(
        ba,
        expr1,
        expr2,
        symbols=symbols,
    )


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
    DNF implicants preserve clauses that prime implicants can merge:

    >>> rule = "(A & B) | (A & !B)"
    >>> dnf_implicants(rule)
    (Hypercube(A=1, B=0), Hypercube(A=1, B=1))
    >>> prime_implicants(rule)
    (Hypercube(A=1),)

    Prime implicants can also expose a consensus condition that is absent from
    the supplied DNF:

    >>> rule = "(A & B) | (!A & C)"
    >>> dnf_implicants(rule)
    (Hypercube(A=0, C=1), Hypercube(A=1, B=1))
    >>> prime_implicants(rule)
    (Hypercube(A=0, C=1), Hypercube(A=1, B=1), Hypercube(B=1, C=1))

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
    ba: BooleanAlgebra = BA,
) -> Implicants:
    """
    Compute prime implicants of a Boolean rule.

    Prime implicants are minimal partial configurations that force a Boolean
    rule to a target value. Each returned `Hypercube` contains only fixed
    components; free components are omitted.

    Examples
    --------
    DNF implicants preserve clauses that prime implicants can merge:

    >>> rule = "(A & B) | (A & !B)"
    >>> dnf_implicants(rule)
    (Hypercube(A=1, B=0), Hypercube(A=1, B=1))
    >>> prime_implicants(rule)
    (Hypercube(A=1),)

    Prime implicants can also expose a consensus condition that is absent from
    the supplied DNF:

    >>> rule = "(A & B) | (!A & C)"
    >>> dnf_implicants(rule)
    (Hypercube(A=0, C=1), Hypercube(A=1, B=1))
    >>> prime_implicants(rule)
    (Hypercube(A=0, C=1), Hypercube(A=1, B=1), Hypercube(B=1, C=1))

    Prime implicants can also be computed for the false value:

    >>> prime_implicants("A & !B", value=0)
    (Hypercube(A=0), Hypercube(B=1))

    Parameters
    ----------
    rule: str, bool, int or boolean.Expression
        Boolean rule to analyze.
    value: {0, 1} (default: 1)
        Target value forced by the returned prime implicants.
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

    expr = _coerce_boolean_rule(ba, rule).simplify()
    constant = _constant_expression_as_bool(ba, expr)
    if constant is not None:
        return (Hypercube({}),) if constant == bool(value) else ()

    symbols = tuple(sorted(expr.symbols, key=str))

    if len(symbols) <= 4:
        _, minterms = _matching_minterms(ba, expr, target=bool(value))
        if not minterms:
            return ()

        prime_cubes = _compute_prime_cubes(tuple(minterms))

    else:
        prime_cubes = _compute_prime_cubes_with_clingo(
            ba,
            expr,
            variables=tuple(str(symbol) for symbol in symbols),
            target=value,
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


def _equivalence_with_bitsets(
    ba: BooleanAlgebra,
    expr1: Expression,
    expr2: Expression,
    *,
    symbols: Tuple[Expression, ...],
) -> bool:

    full_mask, symbol_bits = _boolean_symbol_bitsets(symbols)
    cache: Dict[Expression, int] = {}

    return _expression_truth_bits(
        ba,
        expr1,
        symbol_bits=symbol_bits,
        full_mask=full_mask,
        cache=cache,
    ) == _expression_truth_bits(
        ba,
        expr2,
        symbol_bits=symbol_bits,
        full_mask=full_mask,
        cache=cache,
    )


def _boolean_symbol_bitsets(
    symbols: Tuple[Expression, ...],
) -> Tuple[int, Dict[str, int]]:

    n_symbols = len(symbols)
    n_assignments = 1 << n_symbols
    full_mask = (1 << n_assignments) - 1
    symbol_bits = {}

    for index, symbol in enumerate(symbols):
        block_size = 1 << (n_symbols - index - 1)
        period = 2 * block_size
        unit = ((1 << block_size) - 1) << block_size
        repetitions = full_mask // ((1 << period) - 1)
        symbol_bits[str(symbol)] = unit * repetitions

    return full_mask, symbol_bits


def _expression_truth_bits(
    ba: BooleanAlgebra,
    expr: Expression,
    *,
    symbol_bits: Mapping[str, int],
    full_mask: int,
    cache: Dict[Expression, int],
) -> int:

    if expr in cache:
        return cache[expr]

    if isinstance(expr, _TRUE):
        bits = full_mask

    elif isinstance(expr, _FALSE):
        bits = 0

    else:
        expr_any: Any = expr

        if isinstance(expr_any, ba.NOT):
            bits = full_mask ^ _expression_truth_bits(
                ba,
                expr_any.args[0],
                symbol_bits=symbol_bits,
                full_mask=full_mask,
                cache=cache,
            )

        elif isinstance(expr_any, ba.AND):
            bits = full_mask
            for operand in expr_any.args:
                bits &= _expression_truth_bits(
                    ba,
                    operand,
                    symbol_bits=symbol_bits,
                    full_mask=full_mask,
                    cache=cache,
                )

        elif isinstance(expr_any, ba.OR):
            bits = 0
            for operand in expr_any.args:
                bits |= _expression_truth_bits(
                    ba,
                    operand,
                    symbol_bits=symbol_bits,
                    full_mask=full_mask,
                    cache=cache,
                )

        elif expr_any.isliteral:
            bits = symbol_bits[str(expr_any.obj)]

        else:
            raise ValueError(f"invalid Boolean expression: {expr!r}")

    cache[expr] = bits
    return bits


def _equivalence_with_asp(
    ba: BooleanAlgebra,
    expr1: Expression,
    expr2: Expression,
    *,
    symbols: Tuple[Expression, ...],
) -> bool:

    symbol_ids = {str(symbol): index for index, symbol in enumerate(symbols)}
    expression_ids: Dict[Expression, int] = {}

    root1, facts1, next_expression_id = _clingo_expression_tree_facts(
        ba,
        expr1,
        symbol_ids=symbol_ids,
        expression_ids=expression_ids,
        next_expression_id=0,
    )
    root2, facts2, _ = _clingo_expression_tree_facts(
        ba,
        expr2,
        symbol_ids=symbol_ids,
        expression_ids=expression_ids,
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
        """,
        ]
    )

    control = clingo.Control(["--models=1", "--warn=none"])
    control.add("base", [], program)
    control.ground([("base", [])])

    return bool(control.solve().unsatisfiable)


def _clingo_expression_tree_facts(
    ba: BooleanAlgebra,
    expr: Expression,
    *,
    symbol_ids: Mapping[str, int],
    expression_ids: Dict[Expression, int],
    next_expression_id: int,
) -> Tuple[int, Tuple[str, ...], int]:

    if expr in expression_ids:
        return expression_ids[expr], (), next_expression_id

    expression_id = next_expression_id
    expression_ids[expr] = expression_id
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
            expression_ids=expression_ids,
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
                expression_ids=expression_ids,
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
        for cube in current:
            for index, value in enumerate(cube):
                if value != 0:
                    continue

                neighbor = cube[:index] + (1,) + cube[index + 1 :]
                if neighbor not in current:
                    continue

                used.add(cube)
                used.add(neighbor)
                combined.add(cube[:index] + (None,) + cube[index + 1 :])

        prime_cubes.update(cube for cube in current if cube not in used)
        current = combined

    return frozenset(prime_cubes)


def _compute_prime_cubes_with_clingo(
    ba: BooleanAlgebra,
    expr: Expression,
    *,
    variables: Tuple[str, ...],
    target: Literal[0, 1],
) -> FrozenSet[_BooleanCube]:

    from ._robdd import ROBDD

    robdd = ROBDD._from_rule(
        expr,
        order=variables,
        ba=ba,
    )
    opposite = 1 - target

    program = "\n".join(
        [
            *_clingo_robdd_facts(robdd, variables=variables),
            f"""
            {{ fixed(V, 0); fixed(V, 1) }} :- variable(V).
            :- fixed(V, 0), fixed(V, 1).

            reachable(R) :- root(R).
            reachable(L) :-
                reachable(N), decision(N, V, L, H), not fixed(V, 1).
            reachable(H) :-
                reachable(N), decision(N, V, L, H), not fixed(V, 0).

            :- not reachable({target}).
            :- reachable({opposite}).

            drop_reachable(V, R) :- fixed(V, _), root(R).
            drop_reachable(V, L) :-
                drop_reachable(V, N), decision(N, V, L, H).
            drop_reachable(V, H) :-
                drop_reachable(V, N), decision(N, V, L, H).
            drop_reachable(V, L) :-
                drop_reachable(V, N), decision(N, W, L, H), V != W,
                not fixed(W, 1).
            drop_reachable(V, H) :-
                drop_reachable(V, N), decision(N, W, L, H), V != W,
                not fixed(W, 0).

            :- fixed(V, _), not drop_reachable(V, {opposite}).

            #show fixed/2.
            """,
        ]
    )

    control = clingo.Control(["--models=0", "--warn=none"])
    control.add("base", [], program)
    control.ground([("base", [])])

    cubes: Set[_BooleanCube] = set()
    with control.solve(yield_=True) as handle:
        for model in handle:
            cube: List[Optional[int]] = [None] * len(variables)

            for atom in model.symbols(shown=True):
                variable, fixed_value = atom.arguments
                cube[variable.number] = fixed_value.number

            cubes.add(tuple(cube))

    return frozenset(cubes)


def _clingo_robdd_facts(
    robdd: "ROBDD",
    *,
    variables: Tuple[str, ...],
) -> Tuple[str, ...]:

    variable_indices = {variable: index for index, variable in enumerate(variables)}
    facts = [
        *(f"variable({index})." for index in range(len(variables))),
        f"root({robdd._root_id}).",
    ]
    facts.extend(
        "decision(" f"{node}, {variable_indices[variable]}, {low}, {high}" ")."
        for node, variable, low, high in robdd._iter_nodes()
    )

    return tuple(facts)


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
