#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, Dict, Optional, Tuple, Union

import networkx as nx
from boolean import BooleanAlgebra, Expression

from ._configuration import ConfigurationSet
from ._hypercube import Hypercube
from ._structure import BA, _coerce_boolean_rule
from ._typing import HypercubeLike


class ROBDD:
    """
    Reduced ordered binary decision diagram for one Boolean function.

    A ROBDD represents one scalar function
    ``f: {0, 1}^n -> {0, 1}``, such as one local update rule ``f_i`` of a
    Boolean network. It does not represent the complete vector-valued network
    update function.

    Variables follow the specified order, or lexicographic order by default.
    Nodes with identical low and high branches are removed, and nodes sharing
    the same variable and branches are merged. The resulting diagram is
    canonical for its variable order.

    Examples
    --------
    Build the ROBDD for ``f(A, B) = A | B``:

    >>> robdd = ROBDD("A | B")

    With the variable order ``A < B``, the reduced diagram is::

                 A
              0 / \\ 1
               B   1
            0 / \\ 1
             0   1

    Evaluate the function on a complete configuration:

    >>> robdd.evaluate({"A": 0, "B": 1})
    True

    `can_take_value()` tests whether at least one assignment of the free
    variables produces the requested value. With only ``A = 0`` fixed, both
    terminal values are possible:

    >>> robdd.can_take_value({"A": 0}, 0)
    True
    >>> robdd.can_take_value({"A": 0}, 1)
    True

    Parameters
    ----------
    rule: boolean.Expression, bool, int or str
        Scalar Boolean rule represented by the diagram.
    order: Iterable[str], optional
        Variable order. It must contain every variable exactly once. If None,
        use lexicographic order.
    """

    def __init__(
        self,
        rule: Union[Expression, bool, int, str],
        *,
        order: Optional[Iterable[str]] = None,
    ) -> None:

        self._initialize(rule, order=order, ba=BA)

    def __repr__(self) -> str:
        """Return a compact representation."""

        return f"ROBDD(variables={self.variables!r}, n_nodes={self.n_nodes})"

    def configurations(
        self,
        value: Union[bool, int] = 1,
    ) -> ConfigurationSet:
        """
        Return assignments producing a terminal value.

        Root-to-terminal paths are stored as Boolean subspaces in a
        `ConfigurationSet`; complete configurations are not eagerly
        enumerated.

        Examples
        --------
        >>> configurations = ROBDD("A | B").configurations()
        >>> configurations.enumerate()
        ({'A': 0, 'B': 1}, {'A': 1, 'B': 0}, {'A': 1, 'B': 1})

        Parameters
        ----------
        value: bool or int (default: 1)
            Requested terminal value, 0 or 1.

        Returns
        -------
        ConfigurationSet
            Exact set of complete assignments producing `value`.

        Raises
        ------
        ValueError
            If `value` is not 0 or 1.
        """

        target = _coerce_terminal_value(value)
        subspaces = []
        waiting = [(self._root, 0, 0)]

        while waiting:
            node, fixed_mask, value_mask = waiting.pop()
            if node == target:
                subspaces.append((fixed_mask, value_mask))
                continue
            if node <= 1:
                continue

            variable, low, high = self._nodes[node]
            bit = 1 << variable
            waiting.append((high, fixed_mask | bit, value_mask | bit))
            waiting.append((low, fixed_mask | bit, value_mask))

        return ConfigurationSet._from_encoded_hypercubes(
            self._variables,
            subspaces,
        )

    def to_networkx(self) -> nx.DiGraph:
        """
        Convert the ROBDD into a directed NetworkX graph.

        Decision nodes expose `variable` and `label` attributes. Terminal nodes
        expose `terminal`, `value` and `label`. Edges expose the corresponding
        Boolean branch through `value` and `label`. The root identifier is
        stored in `graph.graph["root"]`.

        Examples
        --------
        >>> graph = ROBDD("A | B").to_networkx()
        >>> graph.graph["root"] == next(
        ...     node for node, data in graph.nodes(data=True)
        ...     if data.get("variable") == "A"
        ... )
        True

        Returns
        -------
        networkx.DiGraph
            Directed acyclic decision graph.
        """

        graph = nx.DiGraph()
        graph.graph["root"] = self._root

        if self._root <= 1:
            graph.add_node(
                self._root,
                terminal=True,
                value=self._root,
                label=str(self._root),
            )
            return graph

        for node, variable, low, high in self._iter_nodes():
            graph.add_node(
                node,
                terminal=False,
                variable=variable,
                label=variable,
            )
            graph.add_edge(node, low, value=0, label="0")
            graph.add_edge(node, high, value=1, label="1")

        for terminal in (0, 1):
            if terminal in graph:
                graph.nodes[terminal].update(
                    terminal=True,
                    value=terminal,
                    label=str(terminal),
                )

        return graph

    @property
    def variables(self) -> Tuple[str, ...]:
        """Return variables in decision order."""

        return self._variables

    @property
    def n_nodes(self) -> int:
        """Return the number of non-terminal decision nodes."""

        return len(self._iter_nodes())

    def evaluate(self, configuration: Mapping[str, Union[bool, int]]) -> bool:
        """
        Evaluate the function on a complete variable assignment.

        Examples
        --------
        >>> ROBDD("A & ~B").evaluate({"A": 1, "B": 0})
        True

        Parameters
        ----------
        configuration: Mapping[str, bool or int]
            Assignment fixing every function variable to 0 or 1. Additional
            keys are ignored.

        Returns
        -------
        bool
            Function value for the assignment.

        Raises
        ------
        ValueError
            If a function variable is missing, free or not Boolean.
        """

        hypercube = Hypercube(
            {
                variable: configuration[variable]
                for variable in self._variables
                if variable in configuration
            }
        )
        missing = [
            variable
            for variable in self._variables
            if variable not in hypercube or not hypercube[variable].is_fixed
        ]
        if missing:
            raise ValueError(
                "configuration must define fixed values for all ROBDD variables: "
                f"{missing}"
            )

        node = self._root
        while node > 1:
            variable, low, high = self._nodes[node]
            value = hypercube[self._variables[variable]].value
            node = high if value else low

        return bool(node)

    def can_take_value(
        self,
        configuration: HypercubeLike,
        value: Union[bool, int],
    ) -> bool:
        """
        Test whether the function can take a value inside a Boolean subspace.

        Missing variables and variables assigned to ``"*"`` are free. The
        method returns True when at least one compatible complete assignment
        produces `value`.

        Examples
        --------
        >>> robdd = ROBDD("A | B")
        >>> robdd.can_take_value({"A": 0}, 1)
        True
        >>> robdd.can_take_value({"A": 0, "B": 0}, 1)
        False

        Parameters
        ----------
        configuration: HypercubeLike
            Boolean configuration.
        value: bool or int
            Requested terminal value, 0 or 1.

        Returns
        -------
        bool
            Whether the requested value is possible in the subspace.

        Raises
        ------
        ValueError
            If `value` is not Boolean or a hypercube value is invalid.
        """

        target = _coerce_terminal_value(value)
        subspace = Hypercube(
            {
                variable: configuration[variable]
                for variable in self._variables
                if variable in configuration
            }
        )
        fixed_values = {
            variable: int(subspace[variable].value)
            for variable in self._variables
            if variable in subspace and subspace[variable].is_fixed
        }

        return self._can_reach_terminal(fixed_values, target)

    def count(self, value: Union[bool, int] = 1) -> int:
        """
        Count complete assignments producing a terminal value.

        Counting uses the ROBDD structure directly and does not enumerate
        complete configurations.

        Examples
        --------
        >>> ROBDD("A | B").count(1)
        3
        >>> ROBDD("A | B").count(0)
        1

        Parameters
        ----------
        value: bool or int (default: 1)
            Requested terminal value, 0 or 1.

        Returns
        -------
        int
            Number of complete variable assignments producing `value`.

        Raises
        ------
        ValueError
            If `value` is not 0 or 1.
        """

        target = _coerce_terminal_value(value)
        cache: Dict[Tuple[int, int], int] = {}
        return self._count_assignments(
            self._root,
            next_variable=0,
            target=target,
            cache=cache,
        )

    def _can_reach_terminal(
        self,
        fixed_values: Mapping[str, int],
        target: int,
    ) -> bool:
        """Test whether a terminal is reachable under fixed variable values."""

        waiting = [self._root]
        visited = set()

        while waiting:
            node = waiting.pop()
            if node == target:
                return True
            if node <= 1 or node in visited:
                continue

            visited.add(node)
            variable, low, high = self._nodes[node]
            variable_name = self._variables[variable]
            if variable_name in fixed_values:
                waiting.append(high if fixed_values[variable_name] else low)
            else:
                waiting.extend((low, high))

        return False

    def _initialize(
        self,
        rule: Union[Expression, bool, int, str],
        *,
        order: Optional[Iterable[str]],
        ba: BooleanAlgebra,
    ) -> None:

        expression = _coerce_boolean_rule(ba, rule)
        self._initialize_expression(expression, order=order, ba=ba)

    def _initialize_expression(
        self,
        expression: Expression,
        *,
        order: Optional[Iterable[str]],
        ba: BooleanAlgebra,
    ) -> None:
        """Initialize directly from a trusted compatible expression."""

        variables = tuple(sorted(str(symbol) for symbol in expression.symbols))
        self._variables = _validate_variable_order(variables, order)
        self._variable_indices = {
            variable: index for index, variable in enumerate(self._variables)
        }
        self._nodes: Dict[int, Tuple[int, int, int]] = {}
        self._unique_nodes: Dict[Tuple[int, int, int], int] = {}
        self._expression_cache: Dict[Expression, int] = {}
        self._apply_cache: Dict[Tuple[str, int, int], int] = {}
        self._negation_cache: Dict[int, int] = {0: 1, 1: 0}
        self._root = self._build(expression, ba)

    @classmethod
    def _from_rule(
        cls,
        rule: Union[Expression, bool, int, str],
        *,
        order: Optional[Iterable[str]] = None,
        ba: BooleanAlgebra,
    ) -> "ROBDD":
        """Build a ROBDD with a specific internal Boolean algebra."""

        robdd = cls.__new__(cls)
        robdd._initialize(rule, order=order, ba=ba)
        return robdd

    @classmethod
    def _from_expression(
        cls,
        expression: Expression,
        *,
        order: Optional[Iterable[str]] = None,
        ba: BooleanAlgebra,
    ) -> "ROBDD":
        """Build a ROBDD from a trusted compatible Boolean expression."""

        robdd = cls.__new__(cls)
        robdd._initialize_expression(expression, order=order, ba=ba)
        return robdd

    def _iter_nodes(self) -> Tuple[Tuple[int, str, int, int], ...]:
        """Return non-terminal nodes reachable from the root."""

        waiting = [self._root]
        visited = set()

        while waiting:
            node = waiting.pop()
            if node <= 1 or node in visited:
                continue

            visited.add(node)
            _, low, high = self._nodes[node]
            waiting.extend((low, high))

        return tuple(
            (
                node,
                self._variables[self._nodes[node][0]],
                self._nodes[node][1],
                self._nodes[node][2],
            )
            for node in sorted(visited)
        )

    @property
    def _root_id(self) -> int:
        """Return the root identifier for internal encoders."""

        return self._root

    def _build(self, expression: Expression, ba: BooleanAlgebra) -> int:

        if expression in self._expression_cache:
            return self._expression_cache[expression]

        expression_any: Any = expression

        if expression is ba.FALSE:
            root = 0
        elif expression is ba.TRUE:
            root = 1
        elif isinstance(expression_any, ba.Symbol):
            variable = str(expression_any.obj)
            root = self._make_node(self._variable_indices[variable], 0, 1)
        elif isinstance(expression_any, ba.NOT):
            root = self._negate(self._build(expression_any.args[0], ba))
        elif isinstance(expression_any, ba.AND):
            root = 1
            for operand in expression_any.args:
                root = self._apply("and", root, self._build(operand, ba))
        elif isinstance(expression_any, ba.OR):
            root = 0
            for operand in expression_any.args:
                root = self._apply("or", root, self._build(operand, ba))
        else:
            raise TypeError(f"unsupported Boolean expression type: {type(expression)}")

        self._expression_cache[expression] = root
        return root

    def _count_assignments(
        self,
        node: int,
        *,
        next_variable: int,
        target: int,
        cache: Dict[Tuple[int, int], int],
    ) -> int:

        key = (node, next_variable)
        if key in cache:
            return cache[key]

        if node <= 1:
            count = 2 ** (len(self._variables) - next_variable) if node == target else 0
        else:
            variable, low, high = self._nodes[node]
            skipped_assignments = 2 ** (variable - next_variable)
            count = skipped_assignments * (
                self._count_assignments(
                    low,
                    next_variable=variable + 1,
                    target=target,
                    cache=cache,
                )
                + self._count_assignments(
                    high,
                    next_variable=variable + 1,
                    target=target,
                    cache=cache,
                )
            )

        cache[key] = count
        return count

    def _make_node(self, variable: int, low: int, high: int) -> int:

        if low == high:
            return low

        key = (variable, low, high)
        if key not in self._unique_nodes:
            node = len(self._nodes) + 2
            self._nodes[node] = key
            self._unique_nodes[key] = node

        return self._unique_nodes[key]

    def _negate(self, node: int) -> int:

        if node in self._negation_cache:
            return self._negation_cache[node]

        variable, low, high = self._nodes[node]
        negated = self._make_node(
            variable,
            self._negate(low),
            self._negate(high),
        )
        self._negation_cache[node] = negated
        self._negation_cache[negated] = node
        return negated

    def _apply(self, operation: str, left: int, right: int) -> int:

        if left > right:
            left, right = right, left

        key = (operation, left, right)
        if key in self._apply_cache:
            return self._apply_cache[key]

        if operation == "and":
            if left == 0:
                result = 0
            elif left == 1:
                result = right
            elif left == right:
                result = left
            else:
                result = self._apply_nodes(operation, left, right)
        elif operation == "or":
            if right == 1:
                result = 1
            elif left == 0:
                result = right
            elif left == right:
                result = left
            else:
                result = self._apply_nodes(operation, left, right)
        else:
            raise ValueError(f"unsupported ROBDD operation: {operation!r}")

        self._apply_cache[key] = result
        return result

    def _apply_nodes(self, operation: str, left: int, right: int) -> int:

        terminal_index = len(self._variables)
        left_index = self._nodes[left][0] if left > 1 else terminal_index
        right_index = self._nodes[right][0] if right > 1 else terminal_index
        variable = min(left_index, right_index)

        if left_index == variable:
            _, left_low, left_high = self._nodes[left]
        else:
            left_low = left_high = left

        if right_index == variable:
            _, right_low, right_high = self._nodes[right]
        else:
            right_low = right_high = right

        return self._make_node(
            variable,
            self._apply(operation, left_low, right_low),
            self._apply(operation, left_high, right_high),
        )


def _validate_variable_order(
    variables: Tuple[str, ...],
    order: Optional[Iterable[str]],
) -> Tuple[str, ...]:
    """Validate and return a complete ROBDD variable order."""

    if order is None:
        return variables

    ordered = tuple(order)
    if any(not isinstance(variable, str) for variable in ordered):
        raise TypeError("order must contain only strings")
    if len(set(ordered)) != len(ordered):
        raise ValueError("order contains duplicated variables")

    missing = sorted(set(variables) - set(ordered))
    unknown = sorted(set(ordered) - set(variables))
    if missing or unknown:
        raise ValueError(
            "order must contain every rule variable exactly once; "
            f"missing={missing}, unknown={unknown}"
        )

    return ordered


def _coerce_terminal_value(value: Union[bool, int]) -> int:
    """Coerce a terminal value to 0 or 1."""

    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int) and value in (0, 1):
        return value

    raise ValueError(f"value must be 0 or 1, received {value!r}")
