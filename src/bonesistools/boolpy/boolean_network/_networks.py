#!/usr/bin/env python

from typing import List, Dict
from ._typing import BooleanNetwork, is_boolean_network_like


def dnf_to_structure(ba, expr, container=frozenset, sort=False):
    """
    Convert a Boolean expression in disjunctive normal form into a nested
    Python structure.

    The returned structure encodes a DNF as a collection of clauses, where each
    clause is a collection of signed literals. Each literal is represented as
    a pair `(symbol, sign)`, where `sign` is True for a positive literal and
    False for a negated literal.

    Parameters
    ----------
    ba: boolean.boolean.BooleanAlgebra
        Boolean algebra used to build the expression.
    expr: boolean.boolean.Expression
        Boolean expression assumed to be in disjunctive normal form.
    container: type (default: frozenset)
        Container used to store clauses and literals.
    sort: bool (default: False)
        Sort literals and clauses before storing them.

    Returns
    -------
    bool | container
        Return True or False for constant expressions. Otherwise, return a
        nested structure representing the DNF.

    Examples
    --------
    The expression `(A & ~B) | C` is converted into:
    `{{('A', True), ('B', False)}, {('C', True)}}`
    """

    from boolean.boolean import _TRUE, _FALSE

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


class BooleanNetworkEnsemble(list):
    """
    Store and analyse an ensemble of Boolean networks sharing the same components.

    Parameters
    ----------
    components: list[str] (optional, default: None)
        Component names expected in each Boolean network. If specified without
        `bns`, initialise an empty ensemble with this component set.
    bns: list[BooleanNetwork] (optional, default: None)
        Boolean networks used to initialise the ensemble. All networks must
        contain the same components.

    Notes
    -----
    Boolean networks are expected to behave as mappings from component names to
    Boolean rules.
    """

    def __init__(
        self,
        components: List[str] = None,
        bns: List[BooleanNetwork] = None,
    ) -> None:

        if components is None and bns is None:
            raise TypeError(
                "cannot instantiate BooleanNetworkEnsemble: "
                "either 'components' or 'bns' must be provided"
            )

        if components is not None and bns is not None:
            raise TypeError(
                "cannot instantiate BooleanNetworkEnsemble: "
                "'components' and 'bns' are mutually exclusive"
            )

        if components is not None:
            super().__init__()
            self.__components = set(components)
            return

        if not isinstance(bns, list):
            raise TypeError(
                f"unsupported argument type for 'bns': expected {list}, "
                f"but received {type(bns)}"
            )

        if len(bns) == 0:
            raise ValueError(
                "cannot infer components from an empty Boolean network list"
            )

        if not all(is_boolean_network_like(bn) for bn in bns):
            raise TypeError(
                "unsupported argument type for 'bns': "
                "all elements must be Boolean network-like objects"
            )

        self.__components = set(bns[0])

        if not all(set(bn) == self.__components for bn in bns[1:]):
            raise ValueError("invalid value: different components between networks")

        super().__init__(bns)

    def __setitem__(self, key, value) -> None:
        if not is_boolean_network_like(value):
            raise TypeError(
                "unsupported argument type: expected Boolean network-like object, "
                f"but received {type(value)}"
            )

        if set(value) != self.__components:
            raise ValueError("invalid value: missing or additional components")

        super().__setitem__(key, value)

    def append(self, other) -> None:
        """
        Append a Boolean network to the ensemble.

        Parameters
        ----------
        other: BooleanNetwork
            Boolean network-like object to append.

        Raises
        ------
        TypeError
            If `other` is not a Boolean network-like object.
        ValueError
            If `other` does not contain exactly the expected components.
        """

        if not is_boolean_network_like(other):
            raise TypeError(
                "unsupported argument type: expected Boolean network-like object, "
                f"but received {type(other)}"
            )

        if set(other) != self.__components:
            raise ValueError("invalid value: missing or additional components")

        super().append(other)

    def insert(self, index, other) -> None:
        """
        Insert a Boolean network into the ensemble.

        Parameters
        ----------
        index: int
            Position where the Boolean network is inserted.
        other: BooleanNetwork
            Boolean network-like object to insert.

        Raises
        ------
        TypeError
            If `other` is not a Boolean network-like object.
        ValueError
            If `other` does not contain exactly the expected components.
        """

        if not is_boolean_network_like(other):
            raise TypeError(
                "unsupported argument type: expected Boolean network-like object, "
                f"but received {type(other)}"
            )

        if set(other) != self.__components:
            raise ValueError("invalid value: missing or additional components")

        super().insert(index, other)

    def get_components(self) -> set:
        """
        Return the set of components shared by all Boolean networks.

        Returns
        -------
        set
            Copy of the component set.
        """

        return self.__components.copy()

    def get_clauses(self) -> Dict[str, list]:
        """
        Return Boolean rules encoded as DNF-like structures for each component.

        Returns
        -------
        dict[str, list]
            Dictionary mapping each component to the list of rules observed
            across the ensemble. Constant rules are kept as booleans; non-constant
            rules are converted into nested DNF structures.
        """

        clauses = {component: [] for component in self.__components}

        for bn in self:
            for component, rule in bn.items():
                if rule is True or rule is False or rule in [0, 1]:
                    clauses[component].append(bool(rule))
                else:
                    clauses[component].append(dnf_to_structure(bn.ba, rule))

        return clauses

    def get_transcription_factors(self) -> dict:
        """
        Count regulators associated with each target across the ensemble.

        Returns
        -------
        dict
            Nested dictionary of the form target -> source -> sign -> count,
            where sign is True for positive literals and False for negative
            literals.
        """

        def get_transcription_factors_from_clause(clause: frozenset) -> dict:
            transcription_factors = {}

            for conjunction in clause:
                for factor, sign in conjunction:
                    if factor not in transcription_factors:
                        transcription_factors[factor] = sign

            return transcription_factors

        clauses = self.get_clauses()
        transcription_factors = {component: {} for component in self.__components}

        for target, clause_set in clauses.items():
            for clause in clause_set:
                if clause is True or clause is False:
                    continue

                for factor, sign in get_transcription_factors_from_clause(
                    clause
                ).items():
                    if factor not in transcription_factors[target]:
                        transcription_factors[target][factor] = {}

                    if sign not in transcription_factors[target][factor]:
                        transcription_factors[target][factor][sign] = 0

                    transcription_factors[target][factor][sign] += 1

        return transcription_factors

    def get_influences(self) -> dict:
        """
        Return regulator-target influences counted across the ensemble.

        Returns
        -------
        dict
            Nested dictionary of the form source -> target -> sign -> count,
            where sign is True for positive influences and False for negative
            influences.
        """

        influences = {component: {} for component in self.__components}
        transcription_factors = self.get_transcription_factors()

        for target, sources in transcription_factors.items():
            for source, influence in sources.items():
                influences[source][target] = influence

        return influences
