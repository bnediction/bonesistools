#!/usr/bin/env python

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Set, cast

from ..boolean_algebra import dnf_implicants
from ..boolean_algebra._robdd import ROBDD

if TYPE_CHECKING:
    from ._network import BooleanNetwork


def _boolean_function_asp_facts(
    network: "BooleanNetwork",
    *,
    clingo: object,
) -> List[str]:
    """Encode unate rules as DNF and non-unate rules as exact ROBDDs."""

    facts = []

    for target, rule in network.items():
        target_symbol = _asp_string(clingo, target)
        if _expression_is_unate(network, rule):
            facts.append(f"unate({target_symbol}).")
            facts.extend(
                _positive_dnf_asp_facts(
                    network,
                    target=target,
                    clingo=clingo,
                )
            )
        else:
            facts.extend(
                _non_unate_bdd_asp_facts(
                    network,
                    target=target,
                    clingo=clingo,
                )
            )

    return facts


def _boolean_function_evaluation_asp_rules() -> str:
    """Return shared MP evaluation rules for DNF and ROBDD encodings."""

    return """
        has_positive_clause(N) :- positive_clause(N,C).

        positive_clause_satisfied(T,N,C) :-
            timepoint(T),
            positive_clause(N,C),
            mp_reach(T,L,W) : positive_clause_literal(N,C,L,W).
        positive_clause_falsified(T,N,C) :-
            timepoint(T),
            positive_clause_literal(N,C,L,W),
            opposite(W,V),
            mp_reach(T,L,V).

        mp_eval(T,N,1) :- positive_clause_satisfied(T,N,C).
        mp_eval(T,N,0) :-
            timepoint(T), unate(N), not has_positive_clause(N).
        mp_eval(T,N,0) :-
            timepoint(T),
            unate(N),
            has_positive_clause(N),
            positive_clause_falsified(T,N,C) : positive_clause(N,C).

        bdd_eval(T,N,0,0) :- timepoint(T), non_unate(N).
        bdd_eval(T,N,1,1) :- timepoint(T), non_unate(N).
        bdd_eval(T,N,B,V) :-
            bdd_node(N,B,L,LO,HI),
            mp_reach(T,L,0),
            bdd_eval(T,N,LO,V).
        bdd_eval(T,N,B,V) :-
            bdd_node(N,B,L,LO,HI),
            mp_reach(T,L,1),
            bdd_eval(T,N,HI,V).
        mp_eval(T,N,V) :-
            non_unate(N),
            bdd_root(N,B),
            bdd_eval(T,N,B,V).
    """


def _asp_string(clingo: object, value: str) -> str:
    """Return a clingo string literal."""

    return str(clingo.String(value))  # pyright: ignore[reportAttributeAccessIssue]


def _expression_is_unate(network: "BooleanNetwork", rule: Any) -> bool:
    """Return whether every symbol occurs with only one polarity."""

    polarities: Dict[str, Set[int]] = {}

    for literal in rule.literalize().get_literals():
        if isinstance(literal, network.ba.NOT):
            component = str(literal.args[0].obj)
            value = 0
        else:
            component = str(literal.obj)
            value = 1

        polarities.setdefault(component, set()).add(value)

    return all(len(values) == 1 for values in polarities.values())


def _positive_dnf_asp_facts(
    network: "BooleanNetwork",
    *,
    target: str,
    clingo: object,
) -> List[str]:
    """Encode positive DNF clauses for one unate function."""

    target_symbol = _asp_string(clingo, target)
    facts = []

    for clause_id, implicant in enumerate(
        dnf_implicants(network[target], value=1, ba=network.ba)
    ):
        facts.append(f"positive_clause({target_symbol}, {clause_id}).")
        for source, source_value in implicant.items():
            if not source_value.is_fixed:
                raise ValueError(
                    "invalid DNF implicant with free component " f"{source!r}"
                )

            facts.append(
                "positive_clause_literal("
                f"{target_symbol}, {clause_id}, "
                f"{_asp_string(clingo, source)}, "
                f"{cast(int, source_value.value)}"
                ")."
            )

    return facts


def _non_unate_bdd_asp_facts(
    network: "BooleanNetwork",
    *,
    target: str,
    clingo: object,
) -> List[str]:
    """Encode one non-unate Boolean function as a reduced ordered BDD."""

    rule = network[target]
    robdd = ROBDD._from_rule(rule, ba=network.ba)
    target_symbol = _asp_string(clingo, target)
    facts = [
        f"non_unate({target_symbol}).",
        f"bdd_root({target_symbol}, {robdd._root_id}).",
    ]

    for node, source, low, high in robdd._iter_nodes():
        facts.append(
            "bdd_node("
            f"{target_symbol}, {node}, {_asp_string(clingo, source)}, "
            f"{low}, {high}"
            ")."
        )

    return facts
