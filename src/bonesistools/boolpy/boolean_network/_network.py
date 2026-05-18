#!/usr/bin/env python

from typing import List, Dict
from ._typing import BooleanNetworkLike, is_boolean_network_like

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional, Union

from boolean.boolean import (
    BooleanAlgebra,
    Expression,
    _TRUE,
    _FALSE,
)
from ..boolean_algebra import rule_to_string, expressions_equivalent, dnf_to_structure

EquivalenceMethod = Literal["simplify", "truth_table"]


class BooleanNetwork(dict):
    """
    Dictionary-like representation of a Boolean network.

    A Boolean network maps each component to a Boolean rule represented
    as a `boolean.py` expression.

    Examples
    --------
    >>> bn = BooleanNetwork(
    ...     {
    ...         "A": "B & !C",
    ...         "B": 1,
    ...         "C": 0,
    ...     }
    ... )
    >>> print(bn)
    A <- B & !C
    B <- 1
    C <- 0
    """

    def __init__(
        self,
        rules: Optional[Mapping[str, Any]] = None,
        ba: Optional[BooleanAlgebra] = None,
        check: bool = True,
    ) -> None:
        """
        Initialize a Boolean network.

        Parameters
        ----------
        rules:
            Mapping associating components to Boolean rules.
        ba:
            Boolean algebra used to parse Boolean expressions.
        check:
            Whether to validate that all symbols referenced by Boolean
            rules are defined as network components.
        """

        self.ba = BooleanAlgebra() if ba is None else ba
        super().__init__()

        if rules is not None:
            for component, rule in rules.items():
                self[component] = self._coerce_rule(rule)

        if check:
            self.validate()

    def __str__(self) -> str:
        """
        Return a readable rule-based representation of the Boolean network.
        """

        return "\n".join(
            f"{component} <- {rule_to_string(rule)}" for component, rule in self.items()
        )

    __repr__ = __str__

    def __eq__(self, other: object) -> bool:
        """
        Test structural equality between Boolean networks.

        Two Boolean networks are equal when they have the same components
        and each associated rule has the same `boolean.py` expression
        structure.

        This does not test logical equivalence. For example, `x` and
        `x | (y & ~y)` may be logically equivalent but not structurally equal.
        Use `equivalent()` for simplified logical comparison.
        """

        if not is_boolean_network_like(other):
            return NotImplemented

        if set(self.keys()) != set(other.keys()):
            return False

        return all(self[component] == other[component] for component in self)

    def __ne__(self, other: object) -> bool:
        eq = self.__eq__(other)

        if eq is NotImplemented:
            return NotImplemented

        return not eq

    @property
    def components(self) -> set[str]:
        """
        Return the set of components defined by the Boolean network.
        """

        return set(self.keys())

    @property
    def symbols(self) -> set[str]:
        """
        Return the set of symbols referenced by Boolean rules.
        """

        return {str(symbol) for rule in self.values() for symbol in rule.symbols}

    @property
    def undefined_symbols(self) -> set[str]:
        """
        Return symbols referenced by rules but not defined as components.
        """

        return self.symbols - self.components

    @property
    def is_closed(self) -> bool:
        """
        Whether all symbols referenced by rules are defined as components.
        """

        return self.symbols <= self.components

    @property
    def rules(self) -> dict[str, str]:
        """
        Return readable string representations of Boolean rules.
        """

        return {component: rule_to_string(rule) for component, rule in self.items()}

    def validate(self) -> None:
        """
        Validate the Boolean network structure.

        Raises
        ------
        ValueError
            If Boolean rules reference undefined components.
        """

        if not self.is_closed:
            raise ValueError(
                "Boolean network rules reference undefined components: "
                f"{sorted(self.undefined_symbols)}"
            )

    def rule(self, component: str) -> str:
        """
        Return a readable string representation of a Boolean rule.
        """

        return rule_to_string(self[component])

    def equivalent(
        self,
        other: object,
        method: EquivalenceMethod = "simplify",
    ) -> bool:
        """
        Test logical equivalence between Boolean networks.

        Two Boolean networks are equivalent when they have the same components
        and each component has an equivalent Boolean rule.

        Parameters
        ----------
        other:
            Boolean-network-like object to compare with.
        method:
            Equivalence strategy used to compare component rules.

            - `"simplify"` compares rules after `boolean.py` simplification.
              This is fast, but not guaranteed to detect all logical
              equivalences.
            - `"truth_table"` exhaustively evaluates each pair of rules on all
              assignments of their symbols. This is exact, but exponential in
              the number of symbols per rule.

        Returns
        -------
        bool
            Whether both Boolean networks are equivalent according to the
            selected method.
        """

        if not is_boolean_network_like(other):
            return NotImplemented

        if set(self.keys()) != set(other.keys()):
            return False

        for component in self:
            rule1 = self._coerce_rule(self[component])
            rule2 = self._coerce_rule(other[component])

            if not expressions_equivalent(
                rule1,
                rule2,
                method=method,
                ba=self.ba,
            ):
                return False

        return True

    def to_bnet(self, file: Optional[Union[str, Path]] = None) -> Optional[str]:
        """
        Export the Boolean network in .bnet format.
        """

        lines = []

        for component, rule in self.items():
            if isinstance(rule, _TRUE):
                rule = "1"
            elif isinstance(rule, _FALSE):
                rule = "0"
            else:
                rule = str(rule).replace("~", "!")

            lines.append(f"{component}, {rule}")

        content = "\n".join(lines) + "\n"

        if file is None:
            return content

        file = Path(file)
        file.write_text(content)

        return None

    def _coerce_rule(self, rule: Any) -> Expression:
        """
        Coerce a Boolean rule into an internal `boolean.py` expression.

        Accepted rule representations include:
        - `boolean.py.Expression`
        - Python booleans (`True`, `False`)
        - integer Boolean constants (`0`, `1`)
        - Boolean expressions encoded as strings

        All rules are converted to `boolean.py` expressions. In particular,
        Boolean constants are mapped to `BooleanAlgebra.TRUE` and
        `BooleanAlgebra.FALSE` rather than Python booleans.
        """
        if isinstance(rule, Expression):
            return rule

        if isinstance(rule, bool):
            return self.ba.TRUE if rule else self.ba.FALSE

        if rule in [0, 1]:
            return self.ba.TRUE if rule else self.ba.FALSE

        if isinstance(rule, str):
            expr = rule.strip()

            if expr == "0":
                return self.ba.FALSE

            if expr == "1":
                return self.ba.TRUE

            return self.ba.parse(expr)

        raise TypeError(
            f"unsupported Boolean rule type: expected str, bool, int or "
            f"Expression, but received {type(rule)}"
        )


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
        bns: List[BooleanNetworkLike] = None,
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
                if isinstance(rule, (_TRUE, _FALSE)):
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
