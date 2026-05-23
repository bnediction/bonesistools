#!/usr/bin/env python

from collections.abc import MutableSequence
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Iterable,
    Mapping,
    Optional,
    Union,
    Dict,
    FrozenSet,
    List,
    Set,
    Tuple,
)
from ._typing import BooleanNetworkLike, is_boolean_network_like

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from pathlib import Path

import networkx as nx

from boolean.boolean import (
    BooleanAlgebra,
    Expression,
    _TRUE,
    _FALSE,
)
from ..boolean_algebra import rule_to_string, expressions_equivalent, dnf_to_structure
from ..interaction_graph import InfluenceGraph

if TYPE_CHECKING:
    from pydot import Dot

EquivalenceMethod = Literal["simplify", "truth_table"]
NodeStyle = Literal["count", "stability"]

from ..plotting import ratio_edge_style, count_node_style, stability_node_style


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

    def __setitem__(self, component: str, rule: Any) -> None:
        """
        Set the Boolean rule associated with a component.

        Components must be strings. Rules are coerced into `boolean.py`
        expressions using `_coerce_rule`.
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
        """

        return "\n".join(
            f"{component} <- {rule_to_string(rule)}" for component, rule in self.items()
        )

    __repr__ = __str__

    def copy(self) -> "BooleanNetwork":
        """
        Return a shallow BooleanNetwork copy.

        The copied network keeps the same Boolean algebra instance and preserves
        the current rule expressions without re-validating closure. This mirrors
        dictionary copy semantics while preserving the BooleanNetwork type.
        """

        return type(self)(self, ba=self.ba, check=False)

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
    def components(self) -> Set[str]:
        """
        Return the set of components defined by the Boolean network.
        """

        return set(self.keys())

    @property
    def symbols(self) -> Set[str]:
        """
        Return the set of symbols referenced by Boolean rules.
        """

        return {str(symbol) for rule in self.values() for symbol in rule.symbols}

    @property
    def undefined_symbols(self) -> Set[str]:
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
    def rules(self) -> Dict[str, str]:
        """
        Return readable string representations of Boolean rules.
        """

        return {component: rule_to_string(rule) for component, rule in self.items()}

    @classmethod
    def from_bnet(
        cls,
        file: Union[str, Path],
        ba: Optional[BooleanAlgebra] = None,
        check: bool = True,
    ) -> "BooleanNetwork":
        from ._parser import read_bnet

        return read_bnet(file, ba=ba, check=check)

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

    def rename(self, old_name: str, new_name: str) -> None:
        """
        Rename a Boolean network component and update all rule references.

        The component key is replaced by `new_name`, and every occurrence of
        `old_name` in Boolean rules is rewritten to `new_name`. Rule
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
        old_name: str
            Component to rename.
        new_name: str
            New component name.

        Raises
        ------
        TypeError
            If `old_name` or `new_name` is not a string.
        KeyError
            If `old_name` is not defined in the Boolean network.
        ValueError
            If `new_name` already exists.
        """

        if not isinstance(old_name, str):
            raise TypeError(
                f"unsupported argument type for 'old_name': "
                f"expected {str} but received {type(old_name)}"
            )

        if not isinstance(new_name, str):
            raise TypeError(
                f"unsupported argument type for 'new_name': "
                f"expected {str} but received {type(new_name)}"
            )

        if new_name == old_name:
            return None

        if old_name not in self:
            raise KeyError(f"component {old_name!r} not found")

        if new_name in self:
            raise ValueError(
                f"invalid argument value for 'new_name': "
                f"component {new_name!r} already exists"
            )

        old_symbol = self.ba.Symbol(old_name)
        new_symbol = self.ba.Symbol(new_name)

        renamed_rules = {}

        for current_component, rule in self.items():
            renamed_component = new_name if current_component == old_name else current_component
            renamed_rules[renamed_component] = rule.subs(
                {old_symbol: new_symbol},
                simplify=False,
            )

        self.clear()

        for renamed_component, rule in renamed_rules.items():
            self[renamed_component] = rule

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

    def influences(self) -> Set[Tuple[str, str, int]]:
        """
        Return the signed influences induced by Boolean rules.

        Each influence is represented as `(source, target, sign)`, where
        `sign` is `1` for positive literals and `-1` for negative literals.
        """

        influences = set()

        for target, rule in self.items():
            for literal in rule.simplify().literalize().get_literals():
                if isinstance(literal, self.ba.NOT):
                    source = literal.args[0].obj
                    sign = -1
                else:
                    source = literal.obj
                    sign = 1

                influences.add((source, target, sign))

        return influences

    def to_networkx(self) -> nx.MultiDiGraph:
        """
        Convert the Boolean network into a signed NetworkX influence graph.

        Nodes correspond to Boolean network components. Edges correspond to
        signed regulatory influences extracted from Boolean rules.
        """

        graph = nx.MultiDiGraph()

        for component in self.components:
            graph.add_node(component)

        for source, target, sign in self.influences():
            graph.add_edge(source, target, sign=sign)

        return graph

    def to_pydot(
        self,
        program: str = "dot",
        edge_style: Optional[Callable[[Mapping[str, Any]], Mapping[str, Any]]] = None,
        **kwargs: Mapping[str, Any],
    ) -> "Dot":
        """
        Convert the influence graph to a styled pydot graph.

        Positive influences are represented as green activating edges, while negative
        influences are represented as red inhibitory edges. Edge signs are inferred
        from the `sign` edge attribute and normalized to -1 or 1.

        Examples
        --------
        >>> dot = ig.to_pydot(rankdir="LR")
        >>> dot.write_pdf("graph.pdf")

        Parameters
        ----------
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting pydot graph.
        edge_style: Callable, optional
            Optional callable used to update edge attributes. The callable receives
            edge attribute dictionaries and must return a mapping of pydot edge
            attributes.
        **kwargs: Mapping[str, Any]
            Keyword arguments passed to the resulting pydot graph using
            `dot.set(key, value)`.

        Returns
        -------
        Dot
            Styled pydot influence graph.
        """

        return InfluenceGraph(self.to_networkx()).to_pydot(
            program=program, edge_style=edge_style, **kwargs
        )

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
            return self.ba.parse(rule_to_string(rule))

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
            f"unsupported argument type for 'rule': "
            f"expected str, bool, int or Expression but received {type(rule)}"
        )


class BooleanNetworkEnsemble(MutableSequence):
    """
    Mutable sequence of Boolean networks sharing the same components.

    A Boolean network ensemble stores Boolean-network-like objects while
    enforcing that all networks are defined over the same component set.

    Parameters
    ----------
    components: Iterable[str] (optional, default: None)
        Components expected in each Boolean network. If specified without
        `bns`, initialise an empty ensemble with this component set.
    bns: Iterable[BooleanNetworkLike] (optional, default: None)
        Boolean-network-like objects used to initialise the ensemble. All
        networks must contain exactly the same components.

    Notes
    -----
    This class behaves like a mutable sequence: networks can be accessed by
    index, appended, inserted, replaced and deleted. Inserted or replaced
    networks are validated before being stored.

    Raises
    ------
    TypeError
        If initialisation arguments are inconsistent or if provided networks
        are not Boolean network-like objects.
    ValueError
        If `bns` is empty or if a network has missing or additional
        components.
    """

    def __init__(
        self,
        components: Optional[Iterable[str]] = None,
        bns: Optional[Iterable[BooleanNetworkLike]] = None,
        ba: Optional[BooleanAlgebra] = None,
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

        self.__ba = BooleanAlgebra() if ba is None else ba
        self._networks = []

        if components is not None:
            self._components = set(components)
            return

        bns = list(bns)

        if len(bns) == 0:
            raise ValueError(
                "cannot infer components from an empty Boolean network collection"
            )

        if not all(is_boolean_network_like(bn) for bn in bns):
            raise TypeError(
                "unsupported argument type for 'bns': "
                "expected iterable of Boolean network-like objects"
            )

        self._components = set(bns[0])

        for bn in bns:
            self._networks.append(self._coerce_network(bn))

    def __len__(self) -> int:
        """
        Return the number of Boolean networks in the ensemble.
        """

        return len(self._networks)

    def __getitem__(self, index):
        """
        Return one or more Boolean networks from the ensemble.

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

    def __setitem__(self, index, value) -> None:
        """
        Replace one or more Boolean networks in the ensemble.

        Parameters
        ----------
        index: int or slice
            Index or slice selecting positions to replace.
        value: BooleanNetworkLike or Iterable[BooleanNetworkLike]
            Boolean-network-like object, or iterable of Boolean-network-like
            objects when `index` is a slice.

        Raises
        ------
        TypeError
            If a provided object is not a Boolean-network-like object.
        ValueError
            If a provided object does not contain exactly the expected
            components.
        """

        if isinstance(index, slice):
            self._networks[index] = [self._coerce_network(bn) for bn in value]
            return

        self._networks[index] = self._coerce_network(value)

    def __delitem__(self, index) -> None:
        """
        Delete one or more Boolean networks from the ensemble.

        Parameters
        ----------
        index: int or slice
            Index or slice selecting Boolean networks to delete.
        """

        del self._networks[index]

    def insert(self, index: int, value: BooleanNetworkLike) -> None:
        """
        Insert a Boolean network into the ensemble.

        Parameters
        ----------
        index: int
            Position where the Boolean network is inserted.
        value: BooleanNetworkLike
            Boolean-network-like object to insert.

        Raises
        ------
        TypeError
            If `value` is not a Boolean-network-like object.
        ValueError
            If `value` does not contain exactly the expected components.
        """

        self._check_network(value)
        self._networks.insert(index, self._coerce_network(value))

    @property
    def components(self) -> FrozenSet[str]:
        """
        Return the set of components shared by all Boolean networks.
        """

        return frozenset(self._components)

    @property
    def ba(self) -> BooleanAlgebra:
        """
        Boolean algebra shared by all Boolean networks in the ensemble.
        """

        return self.__ba

    def rule_structures(self) -> Dict[str, List]:
        """
        Return Boolean rules encoded as DNF-like structures.

        Returns
        -------
        Dict[str, List]
            Dictionary mapping each component to the list of rules observed
            across the ensemble. Constant rules are stored as Python booleans;
            non-constant rules are converted into nested DNF structures.
        """

        rule_structures = {component: [] for component in self._components}

        for bn in self:
            for component, rule in bn.items():
                if isinstance(rule, _TRUE):
                    rule_structures[component].append(True)

                elif isinstance(rule, _FALSE):
                    rule_structures[component].append(False)

                else:
                    rule_structures[component].append(dnf_to_structure(bn.ba, rule))

        return rule_structures

    def regulator_counts(self) -> Dict:
        """
        Count signed regulators associated with each target across the ensemble.

        Each signed regulator is counted at most once per Boolean network,
        regardless of how many times it appears in the rule structure.
        """

        counts = {component: {} for component in self._components}

        for bn in self:
            for target, rule in bn.items():
                influences = set()

                for literal in rule.simplify().literalize().get_literals():
                    if isinstance(literal, bn.ba.NOT):
                        regulator = literal.args[0].obj
                        sign = False
                    else:
                        regulator = literal.obj
                        sign = True

                    influences.add((regulator, sign))

                for regulator, sign in influences:
                    if regulator not in counts[target]:
                        counts[target][regulator] = {}

                    if sign not in counts[target][regulator]:
                        counts[target][regulator][sign] = 0

                    counts[target][regulator][sign] += 1

        return counts

    def influence_counts(self) -> Dict:
        """
        Count signed regulator-target influences across the ensemble.

        Returns
        -------
        Dict
            Nested dictionary of the form source -> target -> sign -> count,
            where sign is True for positive influences and False for negative
            influences.
        """

        influences = {component: {} for component in self._components}
        regulator_counts = self.regulator_counts()

        for target, regulators in regulator_counts.items():
            for source, signs in regulators.items():
                if source not in influences:
                    influences[source] = {}

                influences[source][target] = signs

        return influences

    def to_networkx(self, remove_isolated_nodes: bool = False) -> nx.MultiDiGraph:
        """
        Convert the Boolean network ensemble into an aggregated signed influence graph.

        Nodes correspond to Boolean network components. Each node stores the number
        of distinct Boolean rule structures observed for the corresponding component
        and the stability of the most frequent rule structure. Edges correspond to
        signed influences aggregated across the ensemble and store their occurrence
        count.

        Parameters
        ----------
        remove_isolated_nodes: bool (default: False)
            If True, remove components with no incoming or outgoing influence.

        Returns
        -------
        nx.MultiDiGraph
            Aggregated signed influence graph.
        """

        graph = nx.MultiDiGraph()

        rule_structures = self.rule_structures()

        function_counts = {}
        function_stabilities = {}

        for component, structures in rule_structures.items():
            structure_counts = {}

            for structure in structures:
                if structure not in structure_counts:
                    structure_counts[structure] = 0

                structure_counts[structure] += 1

            function_counts[component] = len(structure_counts)
            function_stabilities[component] = (
                max(structure_counts.values()) / len(self) if structure_counts else 0
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
                        sign=sign,
                        count=count,
                        ratio=count / len(self),
                    )

        if remove_isolated_nodes:
            isolated = list(nx.isolates(graph))
            graph.remove_nodes_from(isolated)

        return graph

    def to_pydot(
        self,
        remove_isolated_nodes: bool = False,
        node_style: Union[
            Literal["count", "stability"],
            Callable[[Mapping[str, Any]], Mapping[str, Any]],
            bool,
            None,
        ] = None,
        min_ratio: float = 0.0,
        show_edge_labels: bool = True,
        edge_style: Union[
            Callable[[float], Mapping[str, Any]],
            bool,
            None,
        ] = ratio_edge_style,
        program: str = "dot",
        **kwargs,
    ) -> "Dot":
        """
        Convert the Boolean network ensemble into an aggregated pydot graph.

        The resulting graph represents signed influences aggregated across the
        ensemble. Edge occurrence counts correspond to the number of Boolean
        networks in which a signed influence is observed. Node attributes include
        both the number of distinct rule structures and the stability of the most
        frequent rule structure.

        Parameters
        ----------
        remove_isolated_nodes: bool (default: False)
            If True, remove components with no incoming or outgoing influence.
        node_style: Literal["count", "stability"] or Callable[[Mapping[str, Any]], Mapping[str, Any]] or bool or None
            Node styling strategy.
            If assigned to `"count"`, nodes are styled according to the number of
            distinct rule structures observed for the component.
            If assigned to `"stability"`, nodes are styled according to the
            frequency of the most common rule structure.
            If assigned to a callable, the callable receives the node attributes
            and must return a mapping of pydot node attributes.
            If assigned to False or None, no additional node styling is applied.
        show_edge_labels: bool (default: True)
            If True, display edge occurrence counts as edge labels.
        min_ratio: float (default: 0.0)
            Minimum edge occurrence ratio required for an influence to be displayed.
            Edges with occurrence ratios strictly smaller than `min_ratio` are removed
            from the aggregated graph before rendering.
        edge_style: Callable[[float], Mapping[str, Any]] or bool or None
            Function used to style edges according to their occurrence ratio in
            the ensemble. The callable receives a ratio between 0 and 1 and must
            return a mapping of pydot edge attributes.
            If assigned to True, use `ratio_edge_style`.
            If assigned to False or None, no additional edge styling is applied.
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting pydot graph.
        **kwargs: Mapping[str, Any]
            Keyword arguments passed to the resulting pydot graph using
            `dot.set(key, value)`.

        Returns
        -------
        Dot
            Aggregated signed influence graph as a pydot object.
        """

        graph = self.to_networkx(remove_isolated_nodes=False)

        if not 0 <= min_ratio <= 1:
            raise ValueError(
                f"invalid argument value for 'min_ratio': "
                f"expected value between 0 and 1 but received {min_ratio!r}"
            )

        edges_to_remove = [
            (u, v) for u, v, data in graph.edges(data=True) if data["ratio"] < min_ratio
        ]

        graph.remove_edges_from(edges_to_remove)

        if remove_isolated_nodes:

            isolated = list(nx.isolates(graph))
            graph.remove_nodes_from(isolated)

        if edge_style is True:
            edge_style = ratio_edge_style

        if node_style == "count":
            node_style = count_node_style

        elif node_style == "stability":
            node_style = stability_node_style

        for _, data in graph.nodes(data=True):

            if node_style not in (None, False):
                data.update(node_style(data))

        for _, _, data in graph.edges(data=True):
            sign = data["sign"]
            count = data["count"]
            ratio = data["ratio"]

            data.update(
                color="green4" if sign is True else "red2",
                arrowhead="normal" if sign is True else "tee",
                penwidth="2",
            )

            if show_edge_labels:
                data["label"] = str(count)

            if edge_style not in (None, False):
                data.update(edge_style(ratio))

        dot = nx.drawing.nx_pydot.to_pydot(graph)

        dot.set_prog(program)

        for key, value in kwargs.items():
            dot.set(key, value)

        return dot

    def _check_network(self, bn: BooleanNetworkLike) -> None:
        """
        Validate a Boolean-network-like object before insertion.

        Parameters
        ----------
        bn: BooleanNetworkLike
            Boolean-network-like object to validate.

        Raises
        ------
        TypeError
            If `bn` is not a Boolean-network-like object.
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
                "invalid argument value for 'bn': " "missing or additional components"
            )

    def _coerce_network(self, bn: BooleanNetworkLike) -> BooleanNetwork:
        self._check_network(bn)
        return BooleanNetwork(bn, ba=self.__ba, check=False)
