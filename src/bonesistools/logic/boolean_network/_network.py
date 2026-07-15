#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from importlib import import_module
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Mapping,
    MutableSequence,
    Optional,
    Set,
    Tuple,
    Union,
    cast,
    overload,
)

import clingo
from boolean.boolean import (
    _FALSE,
    _TRUE,
    BooleanAlgebra,
    Expression,
)

from ..._compat import Literal
from ..._validation import _as_literal, _as_non_negative_integer
from ..boolean_algebra import (
    ConfigurationSet,
    Hypercube,
    PartialBoolean,
    dnf_implicants,
    equivalence,
    prime_implicants,
    rule_to_string,
)
from ..boolean_algebra._structure import (
    Implicants,
    _constant_expression_as_bool,
)
from ..boolean_algebra._typing import (
    BooleanRule,
    Configuration,
    ConfigurationLike,
    HypercubeLike,
    is_configuration_like,
    is_hypercube_like,
)
from ..influence_graph._influence_graph import (
    AggregatedInfluenceGraph,
    InfluenceGraph,
)
from ._dynamics import (
    _automatic_reachable_attractors_backend,
    _automatic_reachable_configurations_backend,
    _bdd_reachability,
    _bdd_reachable_attractors,
    _bdd_reachable_configurations,
    _explicit_reachable_attractors,
    _explicit_reachable_configurations,
    _finite_state_transition,
    _synchronous_reachability,
    _synchronous_reachable_attractors,
    _synchronous_reachable_configurations,
    _validate_complete_hypercube,
)
from ._most_permissive import (
    _most_permissive_reachability,
    _most_permissive_reachable_attractors,
    _most_permissive_reachable_configurations,
    _smallest_closed_hypercube,
)
from ._typing import BooleanNetworkLike, is_boolean_network_like

NodeStyle = Literal["count", "stability"]
_PyBoolNetPrimes = Dict[str, List[List[Dict[str, int]]]]


class BooleanNetwork(Dict[str, Expression]):
    """
    Dictionary-like representation of a Boolean network.

    A Boolean network maps each component to a Boolean rule represented
    as a `boolean.py` expression. Input rules may be provided as strings,
    Boolean constants, integer constants, or existing `boolean.py`
    expressions. They are coerced to internal Boolean algebra expressions on
    assignment.

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
    A <- B & ~C
    B <- 1
    C <- 0

    Rules are stored internally as Boolean algebra expressions, while the
    `rules` property exposes readable strings:

    >>> bn.rules
    {'A': 'B & ~C', 'B': '1', 'C': '0'}

    Parameters
    ----------
    rules: Mapping[str, BooleanRule] (optional, default: None)
        Mapping associating component names to Boolean rules.
    ba: BooleanAlgebra (optional, default: None)
        Boolean algebra used to parse and store Boolean expressions. If `None`,
        a new BooleanAlgebra instance is created.
    check: bool (default: True)
        If `True`, validate that all symbols referenced by rules are defined as
        network components.

    Raises
    ------
    TypeError
        If a component is not a string or a Boolean rule has an unsupported
        type.
    ValueError
        If `check=True` and rules reference undefined components.
    """

    def __init__(
        self,
        rules: Optional[Mapping[str, BooleanRule]] = None,
        ba: Optional[BooleanAlgebra] = None,
        check: bool = True,
    ) -> None:
        """
        Initialize a Boolean network.

        Examples
        --------
        >>> BooleanNetwork({"A": "B", "B": 1}).rules
        {'A': 'B', 'B': '1'}

        Set `check=False` to allow external or temporarily undefined
        regulators:

        >>> bn = BooleanNetwork({"A": "B"}, check=False)
        >>> sorted(bn.undefined_symbols)
        ['B']

        Parameters
        ----------
        rules: Mapping[str, BooleanRule] (optional, default: None)
            Mapping associating component names to Boolean rules.
        ba: BooleanAlgebra (optional, default: None)
            Boolean algebra used to parse and store Boolean expressions. If
            None, a new BooleanAlgebra instance is created.
        check: bool (default: True)
            If `True`, validate that all symbols referenced by rules are defined
            as network components.

        Raises
        ------
        TypeError
            If a component is not a string or a rule has an unsupported type.
        ValueError
            If `check=True` and a rule references an undefined component.
        """

        self.ba = BooleanAlgebra() if ba is None else ba
        super().__init__()

        if rules is not None:
            for component, rule in rules.items():
                self[component] = self._coerce_rule(rule)

        if check:
            self.validate()

    def __setitem__(self, component: str, rule: BooleanRule) -> None:
        """
        Set the Boolean rule associated with a component.

        Components must be strings. Rules are coerced into `boolean.py`
        expressions using `_coerce_rule`.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": 1})
        >>> bn["B"] = "A"
        >>> bn.rules
        {'A': '1', 'B': 'A'}

        Parameters
        ----------
        component: str
            Component name.
        rule: BooleanRule
            Boolean rule to associate with `component`.

        Raises
        ------
        TypeError
            If `component` is not a string or `rule` has an unsupported type.
        """

        if not isinstance(component, str):
            raise TypeError(
                f"unsupported argument type for 'component': "
                f"expected {str} but received {type(component)}"
            )

        super().__setitem__(component, self._coerce_rule(rule))

    def __str__(self) -> str:
        """
        Return a readable rule-based representation of the Boolean network.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": 0, "C": 1})
        >>> print(bn)
        A <- B & ~C
        B <- 0
        C <- 1

        Returns
        -------
        str
            Multiline string representation of network rules.
        """

        return "\n".join(
            f"{component} <- {rule_to_string(rule)}" for component, rule in self.items()
        )

    def __repr__(self) -> str:
        """
        Return a compact representation of the Boolean network.

        Unlike `__str__`, this representation does not dump all rules, which
        keeps notebook displays readable for large networks.

        Examples
        --------
        >>> BooleanNetwork({"A": "B", "B": 1})
        BooleanNetwork(components=2)

        Returns
        -------
        str
            Compact representation with the number of components.
        """

        return f"{type(self).__name__}(components={len(self.components)})"

    def __eq__(self, other: object) -> bool:
        """
        Test exact logical equivalence between Boolean networks.

        Two Boolean networks are equal when they have the same components
        and each associated Boolean rule is logically equivalent.

        Examples
        --------
        >>> bn1 = BooleanNetwork({"A": "B | (B & C)", "B": 0, "C": 1})
        >>> bn2 = BooleanNetwork({"A": "B", "B": 0, "C": 1})
        >>> bn1 == bn2
        True

        Parameters
        ----------
        other: object
            BooleanNetworkLike object to compare against.

        Returns
        -------
        bool or NotImplemented
            True if both networks have the same components and logically
            equivalent rules. Returns NotImplemented for unsupported objects.
        """

        if not is_boolean_network_like(other):
            return NotImplemented

        if set(self.keys()) != set(other.keys()):
            return False

        for component in self:
            rule1 = self._coerce_rule(self[component])
            rule2 = self._coerce_rule(other[component])

            if not equivalence(
                rule1,
                rule2,
                ba=self.ba,
            ):
                return False

        return True

    def __ne__(self, other: object) -> bool:
        """
        Test logical inequality between Boolean networks.

        This is the logical negation of `__eq__` when the other object can be
        interpreted as a BooleanNetworkLike object.

        Parameters
        ----------
        other: object
            BooleanNetworkLike object to compare against.

        Returns
        -------
        bool or NotImplemented
            True if the networks are not logically equivalent. Returns
            NotImplemented for unsupported objects.
        """

        eq = self.__eq__(other)

        if eq is NotImplemented:
            return NotImplemented

        return not eq

    def copy(self) -> "BooleanNetwork":
        """
        Return a shallow BooleanNetwork copy.

        The copied network keeps the same Boolean algebra instance and preserves
        the current rule expressions without re-validating closure. This mirrors
        dictionary copy semantics while preserving the BooleanNetwork type.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B"}, check=False)
        >>> copied = bn.copy()
        >>> copied.rules
        {'A': 'B'}
        >>> copied is bn
        False

        Returns
        -------
        BooleanNetwork
            Shallow copy of the Boolean network.
        """

        return type(self)(self, ba=self.ba, check=False)

    def to_biolqm(self, **kwargs: Any) -> Any:
        """
        Convert the Boolean network to a bioLQM logical model.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `biolqm.load`.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "~B", "B": 1})
        >>> # model = bn.to_biolqm()

        Returns
        -------
        Any
            bioLQM logical model.

        Raises
        ------
        ImportError
            If the bioLQM Python interface is not installed.
        RuntimeError
            If the GINsim executable is unavailable.
        """

        try:
            biolqm = import_module("biolqm")

        except ImportError as error:
            raise ImportError(
                "BooleanNetwork.to_biolqm() requires the bioLQM Python "
                "interface to be installed."
            ) from error

        except FileNotFoundError as error:
            raise RuntimeError(
                "BooleanNetwork.to_biolqm() requires the GINsim executable "
                "to be available."
            ) from error

        with TemporaryDirectory() as directory:
            file = Path(directory) / "network.bnet"
            self.save(file)
            return biolqm.load(str(file), **kwargs)

    def to_minibn(self, **kwargs: Any) -> Any:
        """
        Convert the Boolean network to a colomoto MiniBN network.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `minibn.BooleanNetwork`.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "~B", "B": 1})
        >>> # network = bn.to_minibn()

        Returns
        -------
        Any
            MiniBN Boolean network.

        Raises
        ------
        ImportError
            If `colomoto.minibn` is not installed.
        """

        try:
            from colomoto import minibn

        except ImportError as error:
            raise ImportError(
                "BooleanNetwork.to_minibn() requires colomoto.minibn to be installed."
            ) from error

        return minibn.BooleanNetwork(self._bnet_string(), **kwargs)

    def to_mpbn(self, **kwargs: Any) -> Any:
        """
        Convert the Boolean network to an mpbn network.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `mpbn.MPBooleanNetwork`.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "~B", "B": 1})
        >>> # network = bn.to_mpbn()

        Returns
        -------
        Any
            mpbn Boolean network.

        Raises
        ------
        ImportError
            If `mpbn` is not installed.
        """

        try:
            import mpbn

        except ImportError as error:
            raise ImportError(
                "BooleanNetwork.to_mpbn() requires mpbn to be installed."
            ) from error

        mpbn_constructor = cast(Callable[..., Any], mpbn.MPBooleanNetwork)
        return mpbn_constructor(self._bnet_string(), **kwargs)

    def to_pyboolnet(self) -> _PyBoolNetPrimes:
        """
        Convert the Boolean network to PyBoolNet prime implicants.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B", "B": 1})
        >>> bn.to_pyboolnet()
        {'A': [[{'B': 0}], [{'B': 1}]], 'B': [[], [{}]]}

        Returns
        -------
        dict
            PyBoolNet prime-implicant representation.
        """

        return {
            component: [
                [
                    {
                        source: cast(int, source_value.value)
                        for source, source_value in implicant.items()
                    }
                    for implicant in prime_implicants(
                        rule,
                        value=value,
                        ba=self.ba,
                    )
                ]
                for value in (0, 1)
            ]
            for component, rule in self.items()
        }

    @property
    def components(self) -> Set[str]:
        """
        Return the set of components defined by the Boolean network.

        Examples
        --------
        >>> sorted(BooleanNetwork({"A": "B", "B": 1}).components)
        ['A', 'B']

        Returns
        -------
        Set[str]
            Component names defined by the network.
        """

        return set(self.keys())

    @property
    def symbols(self) -> Set[str]:
        """
        Return the set of symbols referenced by Boolean rules.

        Examples
        --------
        >>> sorted(BooleanNetwork({"A": "B & ~C"}, check=False).symbols)
        ['B', 'C']

        Returns
        -------
        Set[str]
            Symbols referenced by all Boolean rules.
        """

        return {str(symbol) for rule in self.values() for symbol in rule.symbols}

    @property
    def undefined_symbols(self) -> Set[str]:
        """
        Return symbols referenced by rules but not defined as components.

        Examples
        --------
        >>> sorted(BooleanNetwork({"A": "B"}, check=False).undefined_symbols)
        ['B']

        Returns
        -------
        Set[str]
            Referenced symbols that are not network components.
        """

        return self.symbols - self.components

    @property
    def is_closed(self) -> bool:
        """
        Test whether all symbols referenced by rules are defined as components.

        Examples
        --------
        >>> BooleanNetwork({"A": "B", "B": 1}).is_closed
        True
        >>> BooleanNetwork({"A": "B"}, check=False).is_closed
        False

        Returns
        -------
        bool
            True if the network has no undefined symbols.
        """

        return self.symbols <= self.components

    @property
    def rules(self) -> Dict[str, str]:
        """
        Return readable string representations of Boolean rules.

        Examples
        --------
        >>> BooleanNetwork({"A": "B & !C", "B": 0, "C": 1}).rules
        {'A': 'B & ~C', 'B': '0', 'C': '1'}

        Returns
        -------
        Dict[str, str]
            Mapping from components to readable Boolean rule strings.
        """

        return {component: rule_to_string(rule) for component, rule in self.items()}

    def validate(self) -> None:
        """
        Validate the Boolean network structure.

        A valid closed Boolean network only references symbols that are also
        defined as components.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B"}, check=False)
        >>> bn.validate()
        Traceback (most recent call last):
        ...
        ValueError: Boolean network rules reference undefined components: ['B']

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

        Examples
        --------
        >>> BooleanNetwork({"A": "B & !C", "B": 0, "C": 1}).rule("A")
        'B & ~C'

        Parameters
        ----------
        component: str
            Component whose rule should be returned.

        Returns
        -------
        str
            Readable Boolean rule associated with `component`.

        Raises
        ------
        KeyError
            If `component` is not defined in the network.
        """

        return rule_to_string(self[component])

    def simplify(self) -> None:
        """
        Simplify every Boolean rule in place.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B | (B & C)", "B": 0, "C": 1})
        >>> bn.simplify()
        >>> bn.rules
        {'A': 'B', 'B': '0', 'C': '1'}

        Returns
        -------
        None
            Simplifies the Boolean network in place.
        """

        for component, rule in self.items():
            super().__setitem__(component, rule.simplify())

    def rename(self, old: str, new: str) -> None:
        """
        Rename a Boolean network component and update all rule references.

        The component key is replaced by `new`, and every occurrence of `old`
        in Boolean rules is rewritten to `new`. Rule
        rewriting is performed without simplification so that renaming does not
        otherwise alter rule structure.

        Examples
        --------
        >>> bn = BooleanNetwork(
        ...     {
        ...         "A": "B",
        ...         "B": "A & C",
        ...         "C": 1,
        ...     }
        ... )
        >>> bn.rename("A", "X")
        >>> bn.rules
        {'X': 'B', 'B': 'X & C', 'C': '1'}

        Parameters
        ----------
        old: str
            Component to rename.
        new: str
            New component name.

        Raises
        ------
        TypeError
            If `old` or `new` is not a string.
        KeyError
            If `old` is not defined in the Boolean network.
        ValueError
            If `new` already exists or if the renamed Boolean network
            references undefined components.
        """

        if not isinstance(old, str):
            raise TypeError(
                f"unsupported argument type for 'old': "
                f"expected {str} but received {type(old)}"
            )

        if not isinstance(new, str):
            raise TypeError(
                f"unsupported argument type for 'new': "
                f"expected {str} but received {type(new)}"
            )

        if new == old:
            return None

        if old not in self:
            raise KeyError(f"component {old!r} not found")

        if new in self:
            raise ValueError(
                f"invalid argument value for 'new': component {new!r} already exists"
            )

        self.relabel({old: new})

    def relabel(self, mapping: Mapping[str, str]) -> None:
        """
        Relabel several Boolean network components in one operation.

        Components absent from the network are ignored. Rule references to
        relabeled components are updated without simplification so that
        relabeling does not otherwise alter rule structure.

        Examples
        --------
        >>> bn = BooleanNetwork({"Trp53": "Sox2", "Sox2": 1})
        >>> bn.relabel({"Trp53": "TP53", "Sox2": "SOX2"})
        >>> bn.rules
        {'TP53': 'SOX2', 'SOX2': '1'}

        Parameters
        ----------
        mapping: Mapping[str, str]
            Component rename mapping.

        Raises
        ------
        TypeError
            If `mapping` is not a mapping from strings to strings.
        ValueError
            If relabeling would merge two Boolean network components or if the
            relabeled Boolean network references undefined components.
        """

        if not isinstance(mapping, MappingABC):
            raise TypeError(
                f"unsupported argument type for 'mapping': "
                f"expected {Mapping} but received {type(mapping)}"
            )

        for old, new in mapping.items():
            if not isinstance(old, str):
                raise TypeError(
                    "unsupported mapping key type: "
                    f"expected {str} but received {type(old)}"
                )
            if not isinstance(new, str):
                raise TypeError(
                    "unsupported mapping value type: "
                    f"expected {str} but received {type(new)}"
                )

        active_mapping = {
            old: new for old, new in mapping.items() if old in self and old != new
        }
        if not active_mapping:
            return None

        renamed_components = [
            active_mapping.get(component, component) for component in self
        ]
        if len(set(renamed_components)) != len(renamed_components):
            raise ValueError("relabeling would merge Boolean network components")

        substitutions = {
            self.ba.Symbol(old): self.ba.Symbol(new)
            for old, new in active_mapping.items()
        }

        renamed_rules = {}
        for current_component, rule in self.items():
            renamed_component = active_mapping.get(current_component, current_component)
            renamed_rules[renamed_component] = rule.subs(substitutions, simplify=False)

        renamed_network = type(self)(renamed_rules, ba=self.ba, check=True)

        self.clear()

        for renamed_component, rule in renamed_network.items():
            self[renamed_component] = rule

        self.validate()

    def next_state(self, component: str, configuration: ConfigurationLike) -> int:
        """
        Return the next state of one component from a Boolean configuration.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": 1, "C": 0})
        >>> bn.next_state("A", {"A": 0, "B": 1, "C": 0})
        1
        >>> bn.next_state("A", {"A": 0, "B": 1, "C": 1})
        0

        Parameters
        ----------
        component: str
            Component whose Boolean rule is evaluated.
        configuration: ConfigurationLike
            Complete Boolean configuration used to evaluate the rule.

        Returns
        -------
        int
            Next Boolean state of `component`, either 0 or 1.

        Raises
        ------
        ValueError
            If `configuration` is incomplete, contains extra components,
            contains non-Boolean values or if the rule cannot be fully
            evaluated.
        """

        configuration = self._normalize_configuration(configuration)
        substitutions = {
            self.ba.Symbol(name): self.ba.TRUE if value else self.ba.FALSE
            for name, value in configuration.items()
        }
        value = self[component].subs(substitutions).simplify()
        evaluated = _constant_expression_as_bool(self.ba, value)

        if evaluated is not None:
            return int(evaluated)

        raise ValueError(
            "Boolean rule could not be fully evaluated: "
            f"component {component!r} still depends on "
            f"{sorted(map(str, value.symbols))}"
        )

    def next_configuration(
        self,
        configuration: ConfigurationLike,
    ) -> Configuration:
        """
        Return the next Boolean configuration.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": "A", "C": 0})
        >>> bn.next_configuration({"A": 1, "B": 1, "C": 0})
        {'A': 1, 'B': 1, 'C': 0}
        >>> bn.next_configuration({"A": 0, "B": 1, "C": 1})
        {'A': 0, 'B': 0, 'C': 0}

        Parameters
        ----------
        configuration: ConfigurationLike
            Complete Boolean configuration used to evaluate every rule.

        Returns
        -------
        Configuration
            Next Boolean configuration.

        Raises
        ------
        ValueError
            If `configuration` is incomplete, contains extra components,
            contains non-Boolean values or if one rule cannot be fully
            evaluated.
        """

        configuration = self._normalize_configuration(configuration)
        substitutions = {
            self.ba.Symbol(name): self.ba.TRUE if value else self.ba.FALSE
            for name, value in configuration.items()
        }
        next_configuration = {}

        for component, rule in self.items():
            value = rule.subs(substitutions).simplify()
            evaluated = _constant_expression_as_bool(self.ba, value)

            if evaluated is not None:
                next_configuration[component] = int(evaluated)

            else:
                raise ValueError(
                    "Boolean rule could not be fully evaluated: "
                    f"component {component!r} still depends on "
                    f"{sorted(map(str, value.symbols))}"
                )

        return next_configuration

    def is_fixed_point(self, configuration: HypercubeLike) -> bool:
        """
        Test whether a configuration is a fixed point.

        A fixed point is a Boolean configuration `x` such that `f_i(x) = x_i`
        for every component `i`. Configurations containing free Boolean values
        such as `"*"` are accepted and return False.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> bn.is_fixed_point({"A": 1, "B": 1})
        True
        >>> bn.is_fixed_point({"A": 1, "B": 0})
        False

        Parameters
        ----------
        configuration: HypercubeLike
            Boolean configuration to test.

        Returns
        -------
        bool
            Whether `configuration` is a fixed point.

        Raises
        ------
        ValueError
            If the network is not closed or `configuration` does not define
            exactly the network components.
        """

        self.validate()
        if not isinstance(configuration, MappingABC):
            raise TypeError(
                "unsupported argument type for 'configuration': "
                f"expected Mapping but received {type(configuration)}"
            )

        if not all(isinstance(component, str) for component in configuration):
            raise ValueError(
                "invalid argument value for 'configuration': "
                "expected string component names"
            )

        if set(configuration) != self.components:
            raise ValueError(
                "invalid argument value for 'configuration': "
                f"expected components {sorted(self.components)} but received "
                f"{sorted(configuration)}"
            )

        if not is_configuration_like(configuration):
            if is_hypercube_like(configuration):
                hypercube = Hypercube(configuration)
                if not hypercube.is_fully_specified:
                    return False

                configuration = {
                    component: cast(int, hypercube[component].value)
                    for component in self
                }

            else:
                for component in self:
                    self._normalize_boolean_value(
                        configuration[component],
                        name="configuration",
                        component=component,
                    )

        configuration = self._normalize_configuration(
            cast(ConfigurationLike, configuration)
        )

        return self.next_configuration(configuration) == configuration

    def fixed_points(self, limit: int = 0) -> List[Configuration]:
        """
        Return fixed points of the Boolean network.

        Fixed points are fully specified configurations `x` such that
        `f_i(x) = x_i` for every component `i`. The computation is exhaustive
        and therefore exponential in the number of components.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> bn.fixed_points()
        [{'A': 0, 'B': 0}, {'A': 1, 'B': 1}]

        Parameters
        ----------
        limit: int (default: 0)
            Maximum number of fixed points to return. If 0, return all fixed
            points.

        Returns
        -------
        List[Configuration]
            Fixed points as component-value mappings.

        Raises
        ------
        ValueError
            If the network is not closed or `limit` is negative.
        """

        self.validate()

        limit = _as_non_negative_integer(limit, "limit")

        components = list(self.keys())
        fixed_points = []

        for values in product([0, 1], repeat=len(components)):
            configuration = dict(zip(components, values))

            if self.next_configuration(configuration) == configuration:
                fixed_points.append(configuration)

                if limit and len(fixed_points) >= limit:
                    break

        return fixed_points

    def trap_spaces(
        self,
        *,
        kind: Literal["minimal"] = "minimal",
    ) -> Tuple[Hypercube, ...]:
        """
        Enumerate trap spaces of the Boolean network.

        A trap space is a partial Boolean configuration closed under the
        network dynamics: every component fixed in the hypercube is forced by
        its update rule to keep the same value throughout the hypercube.

        Currently only minimal trap spaces are supported. Minimality is with
        respect to hypercube inclusion, so returned trap spaces are
        inclusion-minimal closed hypercubes.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> bn.trap_spaces()
        (Hypercube(A=0, B=0), Hypercube(A=1, B=1))

        Parameters
        ----------
        kind: {"minimal"} (default: "minimal")
            Trap-space family to enumerate.

        Returns
        -------
        tuple of Hypercube
            Trap spaces represented as partial Boolean configurations.

        Raises
        ------
        ValueError
            If `kind` is unsupported or the network is not closed.
        """

        _as_literal(kind, choices=("minimal",), name="kind")

        self.validate()

        return self._minimal_trap_spaces_with_asp()

    def principal_trap_space(
        self,
        configuration: ConfigurationLike,
    ) -> Hypercube:
        """
        Return the smallest trap space containing a configuration.

        The principal trap space of a configuration is the smallest closed
        hypercube, with respect to inclusion, that contains it.

        Examples
        --------
        A fixed point generates a fully specified principal trap space:

        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> bn.principal_trap_space({"A": 0, "B": 0})
        Hypercube(A=0, B=0)

        Components are freed when their fixed values are not preserved
        throughout the enclosing hypercube:

        >>> bn.principal_trap_space({"A": 0, "B": 1})
        Hypercube()

        Parameters
        ----------
        configuration: ConfigurationLike
            Complete Boolean configuration enclosed by the returned trap
            space.

        Returns
        -------
        Hypercube
            Smallest trap space containing `configuration`.

        Raises
        ------
        ValueError
            If the network is not closed or `configuration` is invalid or
            incomplete.
        """

        self.validate()
        resolved_configuration = _validate_complete_hypercube(
            tuple(self.keys()),
            configuration,
            name="configuration",
        )

        return _smallest_closed_hypercube(
            self,
            resolved_configuration,
            relaxed_components=self.keys(),
        )

    def attractors(
        self,
        initial: Optional[HypercubeLike] = None,
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "asynchronous",
    ) -> Tuple[ConfigurationSet, ...]:
        """
        Return attractors, optionally restricted by an initial configuration.

        Attractors are returned as exact `ConfigurationSet` objects. Iterating
        over a returned set yields complete Boolean configurations, while
        `len(...)` returns the number of configurations in the attractor. The
        order of attractors and configurations is unspecified.

        Examples
        --------
        Synchronous dynamics are deterministic, so one initial configuration
        reaches one cycle:

        >>> bn = BooleanNetwork({"A": "~A"})
        >>> attractors = bn.attractors({"A": 0}, update="synchronous")
        >>> attractors[0].enumerate()
        ({'A': 0}, {'A': 1})

        Asynchronous dynamics can branch and reach several terminal strongly
        connected components:

        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> attractors = bn.attractors(
        ...     {"A": 0, "B": 1},
        ...     update="asynchronous",
        ... )
        >>> [attractor.enumerate() for attractor in attractors]
        [({'A': 0, 'B': 0},), ({'A': 1, 'B': 1},)]

        A partial initial configuration is interpreted as the set of all
        compatible complete configurations:

        >>> bn = BooleanNetwork({"A": "B", "B": "A"})
        >>> attractors = bn.attractors({"A": 0})
        >>> [attractor.enumerate() for attractor in attractors]
        [({'A': 0, 'B': 0},), ({'A': 1, 'B': 1},)]

        If no initial configuration is provided, all configurations are
        considered initial:

        >>> attractors = bn.attractors()
        >>> [attractor.enumerate() for attractor in attractors]
        [({'A': 0, 'B': 0},), ({'A': 1, 'B': 1},)]

        Parameters
        ----------
        initial: HypercubeLike, optional
            Initial Boolean configuration. It may be partially specified;
            omitted or free components define multiple compatible
            configurations. If `None`, all configurations are considered
            initial.
        update: {"asynchronous", "synchronous", "general", "most-permissive"}
            (default: "asynchronous")
            Update semantics.

            - `"synchronous"` updates every unstable component at once.
            - `"asynchronous"` updates one unstable component at a time and
              uses transition-guided reduction for symbolic attractor
              detection [1, 2].
            - `"general"` updates any non-empty subset of unstable components.
            - `"most-permissive"` returns reachable minimal trap spaces using
              most-permissive reachability [3].

        Returns
        -------
        tuple of ConfigurationSet
            Reachable attractors. Each `ConfigurationSet` is an exact set of
            complete Boolean configurations.

        Raises
        ------
        ValueError
            If the network is not closed, `initial` is invalid,
            or `update` is invalid.

        References
        ----------
        [1] Benes et al. (2021). Computing bottom SCCs symbolically using
        transition guided reduction. International Conference on Computer
        Aided Verification, 505-528.

        [2] Xie and Beerel (2000). Implicit enumeration of strongly connected
        components and an application to formal verification. IEEE
        Transactions on Computer-Aided Design of Integrated Circuits and
        Systems, 19(10), 1225-1230.

        [3] Pauleve et al. (2020). Reconciling qualitative, abstract, and
        scalable modeling of biological networks. Nature Communications,
        11(1), 4256.
        """

        update = _as_literal(
            update,
            choices=("asynchronous", "synchronous", "general", "most-permissive"),
            name="update",
        )

        self.validate()
        resolved_initial: HypercubeLike
        if initial is None:
            resolved_initial = {}
        else:
            resolved_initial = initial

        initial_configurations = ConfigurationSet(
            tuple(self.keys()),
            [resolved_initial],
        )

        if update == "most-permissive":
            return _most_permissive_reachable_attractors(
                self,
                resolved_initial,
            )

        selected_backend = _automatic_reachable_attractors_backend(
            self,
            resolved_initial,
            initial_configurations,
            update=update,
        )

        if selected_backend == "bdd":
            return _bdd_reachable_attractors(
                self,
                initial_configurations,
                update=update,
            )
        elif update == "synchronous":
            return _synchronous_reachable_attractors(
                self,
                initial_configurations,
            )
        else:
            return _explicit_reachable_attractors(
                self,
                initial_configurations,
                update=update,
            )

    def reachable_configurations(
        self,
        initial: HypercubeLike,
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "most-permissive",
    ) -> Iterator[Configuration]:
        """
        Iterate over configurations reachable from one initial configuration.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": 1, "B": "A"})
        >>> list(bn.reachable_configurations({"A": 0, "B": 0}))
        [{'A': 1, 'B': 0}, {'A': 1, 'B': 1}, {'A': 0, 'B': 0}]

        Parameters
        ----------
        initial: HypercubeLike
            Initial Boolean configuration. It must be fully specified.
        update: {"asynchronous", "synchronous", "general", "most-permissive"}
            (default: "most-permissive")
            Update semantics. Finite-state dynamics use an automatically
            selected explicit or symbolic traversal. Most-permissive dynamics
            use a compact decomposition into closed hypercubes and irreversible
            components [1].

        Returns
        -------
        Iterator[Configuration]
            Complete Boolean configurations reachable from `initial`.
            Iteration order is unspecified.

        Raises
        ------
        ValueError
            If the network is not closed, `initial` is invalid,
            or `update` is invalid.

        References
        ----------
        [1] Pauleve et al. (2020). Reconciling qualitative, abstract, and
        scalable modeling of biological networks. Nature Communications,
        11(1), 4256.
        """

        update = _as_literal(
            update,
            choices=("asynchronous", "synchronous", "general", "most-permissive"),
            name="update",
        )

        self.validate()
        initial_hypercube = _validate_complete_hypercube(
            tuple(self.keys()),
            initial,
            name="initial",
        )

        if update == "synchronous":
            return _synchronous_reachable_configurations(
                self,
                initial_hypercube,
            )

        if update == "most-permissive":
            return _most_permissive_reachable_configurations(
                self,
                initial_hypercube,
            )

        selected_backend = _automatic_reachable_configurations_backend(
            self,
            initial_hypercube,
            update=update,
        )
        if selected_backend == "bdd":
            return _bdd_reachable_configurations(
                self,
                initial_hypercube,
                update=update,
            )

        return _explicit_reachable_configurations(
            self,
            initial_hypercube,
            update=update,
        )

    def transition(
        self,
        source: HypercubeLike,
        target: HypercubeLike,
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "most-permissive",
        quantifier: Literal["exists", "robust", "universal"] = "robust",
    ) -> bool:
        """
        Test a direct transition between two Boolean configurations.

        For most-permissive dynamics [1], a transition is equivalent to
        most-permissive reachability. Consequently, `transition(...)` and
        `reachability(...)` return the same result when
        `update="most-permissive"`.

        Examples
        --------
        Asynchronous dynamics update exactly one unstable component:

        >>> bn = BooleanNetwork({"A": "~A", "B": "~B"})
        >>> bn.transition(
        ...     {"A": 0, "B": 0},
        ...     {"A": 1, "B": 0},
        ...     update="asynchronous",
        ... )
        True
        >>> bn.transition(
        ...     {"A": 0, "B": 0},
        ...     {"A": 1, "B": 1},
        ...     update="asynchronous",
        ... )
        False

        General dynamics may update both components at once:

        >>> bn.transition(
        ...     {"A": 0, "B": 0},
        ...     {"A": 1, "B": 1},
        ...     update="general",
        ... )
        True

        Parameters
        ----------
        source: HypercubeLike
            Source Boolean configuration. It may be partially specified;
            omitted or free components represent all compatible complete
            configurations.
        target: HypercubeLike
            Target Boolean configuration. It may be partially specified;
            omitted or free components represent all compatible complete
            configurations.
        update: {"asynchronous", "synchronous", "general", "most-permissive"}
            (default: "most-permissive")
            Update semantics.
        quantifier: {"exists", "robust", "universal"} (default: "robust")
            Quantification used for partial configurations. If `"exists"`, at
            least one source configuration must have a transition to one target
            configuration. If `"robust"`, every source configuration must have
            a transition to at least one target configuration. If
            `"universal"`, every source configuration must have a transition to
            every target configuration.

        Returns
        -------
        bool
            Whether the direct transition satisfies the selected quantifier.

        Raises
        ------
        ImportError
            If partial-configuration dynamics require the optional dependency
            `dd` and it is not installed.
        ValueError
            If the network is not closed, `source` or `target` is invalid,
            `update` is invalid or `quantifier` is invalid.

        References
        ----------
        [1] Pauleve et al. (2020). Reconciling qualitative, abstract, and
        scalable modeling of biological networks. Nature Communications,
        11(1), 4256.
        """

        update = _as_literal(
            update,
            choices=(
                "asynchronous",
                "synchronous",
                "general",
                "most-permissive",
            ),
            name="update",
        )
        quantifier = _as_literal(
            quantifier,
            choices=("exists", "robust", "universal"),
            name="quantifier",
        )

        self.validate()

        if update == "most-permissive":
            return _most_permissive_reachability(
                self,
                source,
                target,
                quantifier=quantifier,
            )

        return _finite_state_transition(
            self,
            source,
            target,
            update=update,
            quantifier=quantifier,
        )

    def reachability(
        self,
        source: HypercubeLike,
        target: HypercubeLike,
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "most-permissive",
        quantifier: Literal["exists", "robust", "universal"] = "robust",
    ) -> bool:
        """
        Test reachability between two Boolean configurations.

        Most-permissive reachability follows Pauleve et al. [1].

        Examples
        --------
        Most-permissive reachability keeps irreversible components explicit.
        In this example, `A` can irreversibly become 1 from the initial
        configuration, but `B` cannot become 1 while `A` remains 0:

        >>> bn = BooleanNetwork({"A": 1, "B": "A"})
        >>> bn.reachability({"A": 0, "B": 0}, {"A": 1, "B": 0})
        True
        >>> bn.reachability({"A": 0, "B": 0}, {"A": 1, "B": 1})
        True
        >>> bn.reachability({"A": 0, "B": 0}, {"A": 0, "B": 1})
        False

        Asynchronous reachability uses a symbolic BDD backward closure:

        >>> bn.reachability(
        ...     {"A": 0, "B": 0},
        ...     {"A": 1, "B": 1},
        ...     update="asynchronous",
        ... )
        True

        Parameters
        ----------
        source: HypercubeLike
            Source Boolean configuration. It may be partially specified;
            omitted or free components represent all compatible complete
            configurations.
        target: HypercubeLike
            Target Boolean configuration. It may be partially specified;
            omitted or free components represent all compatible complete
            configurations.
        update: {"asynchronous", "synchronous", "general", "most-permissive"}
            (default: "most-permissive")
            Update semantics.
        quantifier: {"exists", "robust", "universal"} (default: "robust")
            Quantification used for partial configurations. If `"exists"`, at
            least one configuration compatible with `source` must reach one
            configuration compatible with `target`. If `"robust"`, every
            configuration compatible with `source` must reach at least one
            configuration compatible with `target`. If `"universal"`, every
            configuration compatible with `source` must reach every
            configuration compatible with `target`.

        Returns
        -------
        bool
            Whether `target` is reachable from `source` under the selected
            update semantics.

        Raises
        ------
        ImportError
            If the selected dynamics require the optional dependency `dd` and
            it is not installed.
        ValueError
            If the network is not closed, `source` or `target` is invalid,
            `update` is invalid or `quantifier` is invalid.

        References
        ----------
        [1] Pauleve et al. (2020). Reconciling qualitative, abstract, and
        scalable modeling of biological networks. Nature Communications,
        11(1), 4256.
        """

        update = _as_literal(
            update,
            choices=(
                "asynchronous",
                "synchronous",
                "general",
                "most-permissive",
            ),
            name="update",
        )
        quantifier = _as_literal(
            quantifier,
            choices=("exists", "robust", "universal"),
            name="quantifier",
        )

        self.validate()

        if update == "synchronous":
            return _synchronous_reachability(
                self,
                source,
                target,
                quantifier=quantifier,
            )

        if update == "asynchronous" or update == "general":
            return _bdd_reachability(
                self,
                source,
                target,
                update=update,
                quantifier=quantifier,
            )

        return _most_permissive_reachability(
            self,
            source,
            target,
            quantifier=quantifier,
        )

    def influences(self) -> Set[Tuple[str, str, int]]:
        """
        Return signed syntactic influences induced by literalized rules.

        Each influence is represented as `(source, target, sign)`, where
        `sign` is `1` for positive literals and `-1` for negative literals.
        Influences are extracted after `boolean.py` simplification and
        literalization with `rule.simplify().literalize().get_literals()`.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": 0, "C": 1})
        >>> sorted(bn.influences())
        [('B', 'A', 1), ('C', 'A', -1)]

        Returns
        -------
        Set[Tuple[str, str, int]]
            Signed regulator-target influences induced by network rules.

        Notes
        -----
        This method reports syntactic/literalized influences. It should not be
        interpreted as a general causal influence analysis: logically
        equivalent rule formulations can expose literals differently, and the
        same regulator may appear with both signs in more complex rules.
        """

        influences = set()

        for target, rule in self.items():
            literals = cast(Iterable[Any], rule.simplify().literalize().get_literals())

            for literal in literals:
                if isinstance(literal, self.ba.NOT):
                    operand = literal.args[0]
                    source = operand.obj
                    sign = -1
                else:
                    source = literal.obj
                    sign = 1

                influences.add((source, target, sign))

        return influences

    def to_influence_graph(self) -> InfluenceGraph:
        """
        Convert the Boolean network into a signed influence graph.

        Nodes correspond to Boolean network components. Edges correspond to
        signed regulatory influences extracted from Boolean rules.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": 0, "C": 1})
        >>> graph = bn.to_influence_graph()
        >>> sorted(graph.nodes)
        ['A', 'B', 'C']
        >>> graph["B"]["A"][0]["sign"]
        1

        Returns
        -------
        InfluenceGraph
            Signed influence graph induced by the Boolean network.
        """

        graph = InfluenceGraph()

        for component in self.components:
            graph.add_node(component)

        for source, target, sign in self.influences():
            graph.add_edge(source, target, sign=1 if sign == 1 else -1)

        return graph

    def save(
        self,
        file: Union[str, Path],
        *,
        format: Literal["auto", "bnet"] = "auto",
    ) -> None:
        """
        Save the Boolean network to a file.

        The output format is inferred from the file extension by default.
        Currently, `.bnet` is supported. Its text representation stores one
        rule per line as `component, rule`, using `!` for negation and `0` or
        `1` for Boolean constants.

        Parameters
        ----------
        file: str or Path
            Output file path.
        format: {"auto", "bnet"} (default: "auto")
            Output format. If `"auto"`, infer it from the file extension.

        Examples
        --------
        >>> bn = BooleanNetwork({"A": "B & ~C", "B": 0, "C": 1})
        >>> bn.save("network.bnet")

        Returns
        -------
        None
            Writes the Boolean network to `file`.

        Raises
        ------
        ValueError
            If `format` is unsupported or cannot be inferred from the file
            extension.
        """

        format = _as_literal(
            format,
            choices=("auto", "bnet"),
            name="format",
        )
        file = Path(file)

        if format == "auto" and file.suffix.lower() != ".bnet":
            raise ValueError(
                "cannot infer Boolean network format from file extension "
                f"{file.suffix!r}; specify 'format' explicitly"
            )

        file.write_text(self._bnet_string())

    def _bnet_string(self) -> str:
        """Serialize the Boolean network in `.bnet` format."""

        lines = []

        for component, rule in self.items():
            if isinstance(rule, _TRUE):
                rule = "1"
            elif isinstance(rule, _FALSE):
                rule = "0"
            else:
                rule = str(rule).replace("~", "!")

            lines.append(f"{component}, {rule}")

        return "\n".join(lines) + "\n"

    def _minimal_trap_spaces_with_asp(self) -> Tuple[Hypercube, ...]:

        components = tuple(sorted(self.components))
        component_indices = {
            component: index for index, component in enumerate(components)
        }
        base_program = self._trap_space_asp_base_program(
            components=components,
            component_indices=component_indices,
        )

        blockers: List[FrozenSet[Tuple[int, int]]] = []
        trap_spaces: Set[FrozenSet[Tuple[int, int]]] = set()

        while True:
            literals = self._solve_maximal_trap_space_literals(
                clingo,
                base_program=base_program,
                blockers=blockers,
                n_components=len(components),
            )

            if literals is None:
                break

            if literals in trap_spaces:
                break

            trap_spaces.add(literals)
            blockers.append(literals)

        return tuple(
            Hypercube({components[index]: value for index, value in sorted(literals)})
            for literals in sorted(trap_spaces, key=self._trap_space_sort_key)
        )

    def _trap_space_asp_base_program(
        self,
        components: Tuple[str, ...],
        component_indices: Mapping[str, int],
    ) -> str:

        facts = [f"component({index})." for index in range(len(components))]
        implicant_id = 0

        for target, rule in self.items():
            target_index = component_indices[target]

            for value in (0, 1):
                for implicant in prime_implicants(
                    rule,
                    value=value,
                    ba=self.ba,
                ):
                    facts.append(f"implicant({target_index}, {value}, {implicant_id}).")

                    for source, source_value in implicant.items():
                        if not source_value.is_fixed:
                            raise ValueError(
                                "invalid prime implicant with free component "
                                f"{source!r}"
                            )

                        facts.append(
                            "implicant_literal("
                            f"{target_index}, {value}, {implicant_id}, "
                            f"{component_indices[source]}, {source_value.value}"
                            ")."
                        )

                    implicant_id += 1

        return "\n".join(
            [
                *facts,
                """
                { fixed(I, 0); fixed(I, 1) } :- component(I).
                :- fixed(I, 0), fixed(I, 1).

                missing_literal(I, V, P) :-
                    implicant_literal(I, V, P, J, W), not fixed(J, W).
                supported(I, V) :-
                    implicant(I, V, P), not missing_literal(I, V, P).
                :- fixed(I, V), not supported(I, V).

                #maximize { 1,I,V : fixed(I, V) }.
                #show fixed/2.
                """,
            ]
        )

    def _solve_maximal_trap_space_literals(
        self,
        clingo: Any,
        *,
        base_program: str,
        blockers: List[FrozenSet[Tuple[int, int]]],
        n_components: int,
    ) -> Optional[FrozenSet[Tuple[int, int]]]:

        program = "\n".join(
            [
                base_program,
                self._trap_space_blocker_program(
                    blockers,
                    n_components=n_components,
                ),
            ]
        )

        control = clingo.Control(
            ["--opt-mode=opt", "--opt-strategy=usc", "--warn=none"]
        )
        control.add("base", [], program)
        control.ground([("base", [])])

        optimal_literals = None
        with control.solve(yield_=True) as handle:
            for model in handle:
                optimal_literals = frozenset(
                    (atom.arguments[0].number, atom.arguments[1].number)
                    for atom in model.symbols(shown=True)
                )

        if optimal_literals is None:
            return None

        return optimal_literals

    def _trap_space_blocker_program(
        self,
        blockers: List[FrozenSet[Tuple[int, int]]],
        *,
        n_components: int,
    ) -> str:

        lines = []

        for blocker_id, literals in enumerate(blockers):
            lines.append(f"blocker({blocker_id}).")

            for component in range(n_components):
                for value in (0, 1):
                    if (component, value) not in literals:
                        lines.append(
                            f"outside_literal({blocker_id}, {component}, {value})."
                        )

        if blockers:
            lines.extend(
                [
                    "outside(B) :- outside_literal(B, I, V), fixed(I, V).",
                    ":- blocker(B), not outside(B).",
                ]
            )

        return "\n".join(lines)

    def _normalize_configuration(
        self,
        configuration: ConfigurationLike,
    ) -> Configuration:

        if not isinstance(configuration, MappingABC):
            raise TypeError(
                "unsupported argument type for 'configuration': "
                f"expected Mapping but received {type(configuration)}"
            )

        if not all(isinstance(component, str) for component in configuration):
            raise ValueError(
                "invalid argument value for 'configuration': "
                "expected string component names"
            )

        if set(configuration) != self.components:
            raise ValueError(
                "invalid argument value for 'configuration': "
                f"expected components {sorted(self.components)} but received "
                f"{sorted(configuration)}"
            )

        if not is_configuration_like(configuration):
            for component in self:
                self._normalize_boolean_value(
                    configuration[component],
                    name="configuration",
                    component=component,
                )

            raise ValueError(
                "invalid argument value for 'configuration': "
                "expected a concrete Boolean configuration"
            )

        return {
            component: self._normalize_boolean_value(
                configuration[component],
                name="configuration",
                component=component,
            )
            for component in self
        }

    def _coerce_rule(self, rule: BooleanRule) -> Expression:
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

        Examples
        --------
        >>> bn = BooleanNetwork({"A": 1})
        >>> bn._coerce_rule("A")
        Symbol('A')

        Parameters
        ----------
        rule: BooleanRule
            Boolean rule to coerce.

        Returns
        -------
        Expression
            Boolean algebra expression representing `rule`.

        Raises
        ------
        TypeError
            If `rule` has an unsupported type.
        """
        if isinstance(rule, Expression):
            return self.ba.parse(rule_to_string(rule))

        if isinstance(rule, bool):
            return cast(Expression, self.ba.TRUE if rule else self.ba.FALSE)

        if rule in [0, 1]:
            return cast(Expression, self.ba.TRUE if rule else self.ba.FALSE)

        if isinstance(rule, str):
            expr = rule.strip()

            if expr == "0":
                return cast(Expression, self.ba.FALSE)

            if expr == "1":
                return cast(Expression, self.ba.TRUE)

            return self.ba.parse(expr)

        raise TypeError(
            f"unsupported argument type for 'rule': "
            f"expected str, bool, int or Expression but received {type(rule)}"
        )

    @staticmethod
    def _normalize_boolean_value(
        value: Any,
        name: str,
        component: str,
    ) -> int:

        if isinstance(value, PartialBoolean):
            value = value.value

        if value not in [0, 1]:
            raise ValueError(
                f"invalid argument value for '{name}': "
                f"expected 0 or 1 for component {component!r} "
                f"but received {value!r}"
            )

        return int(value)

    @staticmethod
    def _trap_space_sort_key(
        literals: FrozenSet[Tuple[int, int]],
    ) -> Tuple[int, Tuple[Tuple[int, int], ...]]:

        return len(literals), tuple(sorted(literals))


class BooleanNetworkEnsemble(MutableSequence[BooleanNetwork]):
    """
    Mutable sequence of Boolean networks sharing the same components.

    A Boolean network ensemble stores BooleanNetworkLike objects while
    enforcing that all networks are defined over the same component set.

    Examples
    --------
    >>> ensemble = BooleanNetworkEnsemble(
    ...     {"A": "B", "B": 1},
    ...     {"A": 0, "B": "A"},
    ... )
    >>> len(ensemble)
    2
    >>> sorted(ensemble.components)
    ['A', 'B']

    Parameters
    ----------
    *networks: BooleanNetworkLike
        Boolean networks used to initialise the ensemble. All networks must
        contain exactly the same components.
    components: Iterable[str] (optional, default: None)
        Components expected in each Boolean network. If specified without
        `networks`, initialise an empty ensemble with this component set.

    Notes
    -----
    This class behaves like a mutable sequence: networks can be accessed by
    index, appended, inserted, replaced and deleted. Inserted or replaced
    networks are validated before being stored.

    Raises
    ------
    TypeError
        If initialisation arguments are inconsistent or if provided networks
        are not BooleanNetworkLike objects.
    ValueError
        If a network has missing or additional components.
    """

    def __init__(
        self,
        *networks: BooleanNetworkLike,
        components: Optional[Iterable[str]] = None,
        ba: Optional[BooleanAlgebra] = None,
    ) -> None:

        if components is None and not networks:
            raise TypeError(
                "cannot instantiate BooleanNetworkEnsemble: "
                "either 'components' or at least one network must be provided"
            )

        if components is not None and networks:
            raise TypeError(
                "cannot instantiate BooleanNetworkEnsemble: "
                "'components' and positional networks are mutually exclusive"
            )

        self.__ba = BooleanAlgebra() if ba is None else ba
        self._networks: List[BooleanNetwork] = []

        if components is not None:
            self._components: Set[str] = set(components)
            return

        if not all(is_boolean_network_like(network) for network in networks):
            raise TypeError(
                "unsupported positional argument type: "
                "expected Boolean network-like objects"
            )

        self._components = set(networks[0])

        for network in networks:
            self._networks.append(self._coerce_network(network))

    def __len__(self) -> int:
        """
        Return the number of Boolean networks in the ensemble.

        Examples
        --------
        >>> len(BooleanNetworkEnsemble({"A": 1}, {"A": 0}))
        2

        Returns
        -------
        int
            Number of Boolean networks stored in the ensemble.
        """

        return len(self._networks)

    @overload
    def __getitem__(self, index: int) -> BooleanNetwork: ...

    @overload
    def __getitem__(self, index: slice) -> List[BooleanNetwork]: ...

    def __getitem__(
        self,
        index: Union[int, slice],
    ) -> Union[BooleanNetwork, List[BooleanNetwork]]:
        """
        Return one or more Boolean networks from the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": 1}, {"A": 0})
        >>> ensemble[0].rules
        {'A': '1'}

        Parameters
        ----------
        index: int or slice
            Index or slice selecting Boolean networks.

        Returns
        -------
        BooleanNetworkLike or List[BooleanNetworkLike]
            Selected Boolean network or List of Boolean networks.
        """

        return self._networks[index]

    @overload
    def __setitem__(self, index: int, value: BooleanNetworkLike) -> None: ...

    @overload
    def __setitem__(
        self,
        index: slice,
        value: Iterable[BooleanNetworkLike],
    ) -> None: ...

    def __setitem__(
        self,
        index: Union[int, slice],
        value: Union[BooleanNetworkLike, Iterable[BooleanNetworkLike]],
    ) -> None:
        """
        Replace one or more Boolean networks in the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": 1}, {"A": 0})
        >>> ensemble[0] = {"A": "0"}
        >>> ensemble[0].rules
        {'A': '0'}

        Parameters
        ----------
        index: int or slice
            Index or slice selecting positions to replace.
        value: BooleanNetworkLike or Iterable[BooleanNetworkLike]
            BooleanNetworkLike object, or iterable of BooleanNetworkLike
            objects when `index` is a slice.

        Raises
        ------
        TypeError
            If a provided object is not a BooleanNetworkLike object.
        ValueError
            If a provided object does not contain exactly the expected
            components.
        """

        if isinstance(index, slice):
            networks = cast(Iterable[BooleanNetworkLike], value)
            self._networks[index] = [self._coerce_network(bn) for bn in networks]
            return

        self._networks[index] = self._coerce_network(cast(BooleanNetworkLike, value))

    def __delitem__(self, index: Union[int, slice]) -> None:
        """
        Delete one or more Boolean networks from the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": 1}, {"A": 0})
        >>> del ensemble[0]
        >>> len(ensemble)
        1

        Parameters
        ----------
        index: int or slice
            Index or slice selecting Boolean networks to delete.
        """

        del self._networks[index]

    def to_biolqm(self, **kwargs: Any) -> List[Any]:
        """
        Convert all Boolean networks to bioLQM logical models.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `BooleanNetwork.to_biolqm`.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": "~B", "B": 1})
        >>> # models = ensemble.to_biolqm()

        Returns
        -------
        list
            bioLQM logical models.
        """

        return [network.to_biolqm(**kwargs) for network in self]

    def to_minibn(self, **kwargs: Any) -> List[Any]:
        """
        Convert all Boolean networks to colomoto MiniBN networks.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `BooleanNetwork.to_minibn`.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": "~B", "B": 1})
        >>> # networks = ensemble.to_minibn()

        Returns
        -------
        list
            MiniBN Boolean networks.
        """

        return [network.to_minibn(**kwargs) for network in self]

    def to_mpbn(self, **kwargs: Any) -> List[Any]:
        """
        Convert all Boolean networks to mpbn networks.

        Parameters
        ----------
        **kwargs: Any
            Keyword arguments forwarded to `BooleanNetwork.to_mpbn`.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": "~B", "B": 1})
        >>> # networks = ensemble.to_mpbn()

        Returns
        -------
        list
            mpbn Boolean networks.
        """

        return [network.to_mpbn(**kwargs) for network in self]

    def to_pyboolnet(self) -> List[_PyBoolNetPrimes]:
        """
        Convert all Boolean networks to PyBoolNet prime implicants.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": "B", "B": 1})
        >>> ensemble.to_pyboolnet()
        [{'A': [[{'B': 0}], [{'B': 1}]], 'B': [[], [{}]]}]

        Returns
        -------
        list
            PyBoolNet prime-implicant representations.
        """

        return [network.to_pyboolnet() for network in self]

    def to_influence_graph(
        self,
        drop_isolates: bool = False,
    ) -> AggregatedInfluenceGraph:
        """
        Convert the Boolean network ensemble into an aggregated influence graph.

        Nodes correspond to Boolean network components. Each node stores the
        number of distinct DNF implicant structures observed for the corresponding
        component and the stability of the most frequent structure. Edges correspond
        to signed influences aggregated across the ensemble and store their occurrence
        count.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B", "B": 1},
        ...     {"A": "B", "B": 1},
        ... )
        >>> graph = ensemble.to_influence_graph()
        >>> graph.edge_count("B", "A", sign=1)
        2
        >>> graph.edge_frequency("B", "A", sign=1)
        1.0

        Parameters
        ----------
        drop_isolates: bool (default: False)
            If `True`, remove components with no incoming or outgoing influence.

        Returns
        -------
        AggregatedInfluenceGraph
            Aggregated signed influence graph.

        Raises
        ------
        ValueError
            If the ensemble is empty.
        """

        if len(self) == 0:
            raise ValueError("cannot aggregate an empty Boolean network ensemble")

        graph = AggregatedInfluenceGraph(total=len(self))

        component_implicants = self.dnf_implicants()

        function_counts = {}
        function_stabilities = {}

        for component, implicants_by_network in component_implicants.items():
            implicant_counts = {}

            for implicants in implicants_by_network:
                if implicants not in implicant_counts:
                    implicant_counts[implicants] = 0

                implicant_counts[implicants] += 1

            function_counts[component] = len(implicant_counts)
            function_stabilities[component] = (
                max(implicant_counts.values()) / len(self) if implicant_counts else 0
            )

        influence_counts = self.influence_counts()

        for component in self.components:
            graph.add_node(
                component,
                function_count=function_counts[component],
                function_stability=function_stabilities[component],
            )

        for source, targets in influence_counts.items():
            for target, signs in targets.items():
                for sign, count in signs.items():
                    graph.add_edge(
                        source,
                        target,
                        sign=1 if sign else -1,
                        count=count,
                    )

        if drop_isolates:
            isolated = [node for node, degree in graph.degree() if degree == 0]
            graph.remove_nodes_from(isolated)

        return graph

    def insert(self, index: int, value: BooleanNetworkLike) -> None:
        """
        Insert a Boolean network into the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(components=["A"])
        >>> ensemble.insert(0, {"A": 1})
        >>> ensemble[0].rules
        {'A': '1'}

        Parameters
        ----------
        index: int
            Position where the Boolean network is inserted.
        value: BooleanNetworkLike
            BooleanNetworkLike object to insert.

        Raises
        ------
        TypeError
            If `value` is not a BooleanNetworkLike object.
        ValueError
            If `value` does not contain exactly the expected components.
        """

        self._check_network(value)
        self._networks.insert(index, self._coerce_network(value))

    def simplify(self) -> None:
        """
        Simplify every Boolean rule in every network in place.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B | (B & C)", "B": 0, "C": 1},
        ...     {"A": "B & 1", "B": 1, "C": 0},
        ... )
        >>> ensemble.simplify()
        >>> [bn.rule("A") for bn in ensemble]
        ['B', 'B']

        Returns
        -------
        None
            Simplifies all Boolean networks in place.
        """

        for network in self:
            network.simplify()

    @property
    def components(self) -> FrozenSet[str]:
        """
        Return the set of components shared by all Boolean networks.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble({"A": 1, "B": 0})
        >>> sorted(ensemble.components)
        ['A', 'B']

        Returns
        -------
        FrozenSet[str]
            Components expected in every Boolean network of the ensemble.
        """

        return frozenset(self._components)

    @property
    def ba(self) -> BooleanAlgebra:
        """
        Boolean algebra shared by all Boolean networks in the ensemble.

        Returns
        -------
        BooleanAlgebra
            Boolean algebra used to parse and store network rules.
        """

        return self.__ba

    def dnf_implicants(
        self,
        value: Literal[0, 1] = 1,
    ) -> Dict[str, List[Implicants]]:
        """
        Return DNF implicants for each component across the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B", "B": 1, "C": 0},
        ...     {"A": "B & C", "B": 1, "C": 0},
        ... )
        >>> implicants = ensemble.dnf_implicants()
        >>> implicants == {
        ...     "A": [(Hypercube(B=1),), (Hypercube(B=1, C=1),)],
        ...     "B": [(Hypercube(),), (Hypercube(),)],
        ...     "C": [(), ()],
        ... }
        True

        Parameters
        ----------
        value: {0, 1} (default: 1)
            Target value forced by the returned DNF implicants.

        Returns
        -------
        Dict[str, List[Implicants]]
            Dictionary mapping each component to one tuple of DNF implicants per
            Boolean network.

        See Also
        --------
        BooleanNetworkEnsemble.prime_implicants
            Return prime implicants across the ensemble.
        """

        if value not in (0, 1):
            raise ValueError(
                "invalid argument value for 'value': "
                f"expected 0 or 1 but received {value!r}"
            )

        component_implicants: Dict[str, List[Implicants]] = {
            component: [] for component in self._components
        }

        for bn in self:
            for component, rule in bn.items():
                component_implicants[component].append(
                    dnf_implicants(rule, value=value, ba=bn.ba)
                )

        return component_implicants

    def prime_implicants(
        self,
        value: Literal[0, 1] = 1,
    ) -> Dict[str, List[Implicants]]:
        """
        Return prime implicants for each component across the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B | C", "B": 1, "C": 0},
        ...     {"A": "B & C", "B": 1, "C": 0},
        ... )
        >>> ensemble.prime_implicants()["A"]
        [(Hypercube(B=1), Hypercube(C=1)), (Hypercube(B=1, C=1),)]

        Parameters
        ----------
        value: {0, 1} (default: 1)
            Target value forced by the returned prime implicants.

        Returns
        -------
        Dict[str, List[Implicants]]
            Dictionary mapping each component to one tuple of prime implicants
            per Boolean network.

        See Also
        --------
        BooleanNetworkEnsemble.dnf_implicants
            Return DNF implicants across the ensemble.
        """

        if value not in (0, 1):
            raise ValueError(
                "invalid argument value for 'value': "
                f"expected 0 or 1 but received {value!r}"
            )

        component_implicants: Dict[str, List[Implicants]] = {
            component: [] for component in self._components
        }

        for bn in self:
            for component, rule in bn.items():
                component_implicants[component].append(
                    prime_implicants(
                        rule,
                        value=value,
                        ba=bn.ba,
                    )
                )

        return component_implicants

    def regulator_counts(self) -> Dict[str, Dict[str, Dict[bool, int]]]:
        """
        Count signed regulators associated with each target across the ensemble.

        Each signed regulator is counted at most once per Boolean network,
        regardless of how many times it appears in the rule structure.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B", "B": 1, "C": 0},
        ...     {"A": "B & ~C", "B": 1, "C": 0},
        ... )
        >>> ensemble.regulator_counts()["A"]["B"][True]
        2
        >>> ensemble.regulator_counts()["A"]["C"][False]
        1

        Returns
        -------
        Dict
            Nested dictionary of the form target -> regulator -> sign -> count,
            where sign is True for positive regulation and False for negative
            regulation.
        """

        counts: Dict[str, Dict[str, Dict[bool, int]]] = {
            component: {} for component in self._components
        }

        for bn in self:
            for target, rule in bn.items():
                influences: Set[Tuple[str, bool]] = set()

                literals = cast(
                    Iterable[Any],
                    rule.simplify().literalize().get_literals(),
                )

                for literal in literals:
                    if isinstance(literal, bn.ba.NOT):
                        operand = cast(Any, literal.args[0])
                        regulator = cast(str, operand.obj)
                        sign = False
                    else:
                        regulator = cast(str, literal.obj)
                        sign = True

                    influences.add((regulator, sign))

                for regulator, sign in influences:
                    if regulator not in counts[target]:
                        counts[target][regulator] = {}

                    if sign not in counts[target][regulator]:
                        counts[target][regulator][sign] = 0

                    counts[target][regulator][sign] += 1

        return counts

    def influence_counts(self) -> Dict[str, Dict[str, Dict[bool, int]]]:
        """
        Count signed regulator-target influences across the ensemble.

        Examples
        --------
        >>> ensemble = BooleanNetworkEnsemble(
        ...     {"A": "B", "B": 1, "C": 0},
        ...     {"A": "B & ~C", "B": 1, "C": 0},
        ... )
        >>> ensemble.influence_counts()["B"]["A"][True]
        2

        Returns
        -------
        Dict
            Nested dictionary of the form source -> target -> sign -> count,
            where sign is True for positive influences and False for negative
            influences.
        """

        influences: Dict[str, Dict[str, Dict[bool, int]]] = {
            component: {} for component in self._components
        }
        regulator_counts = self.regulator_counts()

        for target, regulators in regulator_counts.items():
            for source, signs in regulators.items():
                if source not in influences:
                    influences[source] = {}

                influences[source][target] = signs

        return influences

    def _check_network(self, bn: BooleanNetworkLike) -> None:
        """
        Validate a BooleanNetworkLike object before insertion.

        Parameters
        ----------
        bn: BooleanNetworkLike
            BooleanNetworkLike object to validate.

        Raises
        ------
        TypeError
            If `bn` is not a BooleanNetworkLike object.
        ValueError
            If `bn` does not contain exactly the expected components.
        """

        if not is_boolean_network_like(bn):
            raise TypeError(
                "unsupported argument type for 'bn': "
                "expected Boolean network-like object "
                f"but received {type(bn)}"
            )

        if set(bn) != self._components:
            raise ValueError(
                "invalid argument value for 'bn': missing or additional components"
            )

    def _coerce_network(self, bn: BooleanNetworkLike) -> BooleanNetwork:
        """
        Convert a BooleanNetworkLike object to a BooleanNetwork.

        Parameters
        ----------
        bn: BooleanNetworkLike
            BooleanNetworkLike object to coerce.

        Returns
        -------
        BooleanNetwork
            Coerced Boolean network using the ensemble Boolean algebra.

        Raises
        ------
        TypeError
            If `bn` is not BooleanNetworkLike.
        ValueError
            If `bn` does not contain exactly the expected components.
        """

        self._check_network(bn)
        return BooleanNetwork(bn, ba=self.__ba, check=False)
