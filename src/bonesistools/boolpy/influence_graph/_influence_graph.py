#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingInstance
from itertools import product
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Mapping,
    MutableMapping,
    NoReturn,
    Optional,
    Set,
    Tuple,
    Union,
    cast,
    overload,
)

import networkx as nx

from ..._compat import Literal
from ..._validation import (
    _as_non_negative_integer,
    _as_positive_integer,
    _as_probability,
)
from ..._warnings import _warn_deprecated
from ..plotting._graphviz import _networkx_to_graphviz
from ..plotting._styles import ratio_edge_style
from ..plotting._svg import SvgLength, scale_svg

if TYPE_CHECKING:
    from pydot import Dot

    from ..boolean_network import BooleanNetwork, BooleanNetworkEnsemble

CircuitSign = Literal[-1, 1]
InfluenceSign = Literal[-1, 1, "+", "-", "positive", "negative"]
Direction = Literal["upstream", "downstream", "both"]
CollapseMode = Literal["family", "feedback", "both"]

StructuralSignature = Tuple[
    FrozenSet[Tuple[str, int]],
    FrozenSet[Tuple[str, int]],
]

if TYPE_CHECKING:
    _MultiDiGraphBase = nx.MultiDiGraph[Any]
else:
    _MultiDiGraphBase = nx.MultiDiGraph


class InfluenceGraph(_MultiDiGraphBase):
    """
    Signed influence graph.

    InfluenceGraph is a domain-specific NetworkX MultiDiGraph for signed
    regulatory influence graphs. It enforces signed edges and provides analyses
    specific to logical and regulatory graphs, such as feedback circuits, signed
    autoregulations, strongly connected components and graph collapse.

    Edge signs must be stored in the edge attribute `sign`, using either:
        - 1, "+", or "positive" for positive influences,
        - -1, "-", or "negative" for negative influences.

    At most one positive edge and one negative edge may exist between a pair of
    nodes.
    """

    __slots__ = ()

    def __init__(
        self,
        graph: Optional[Union[nx.DiGraph[Any], nx.MultiDiGraph[Any]]] = None,
        **attr: Any,
    ) -> None:
        """
        Initialize a signed influence graph.

        Parameters
        ----------
        graph: nx.DiGraph, optional
            Signed directed graph. Edges must define a `sign` attribute equal to
            -1 or 1. MultiDiGraph inputs are accepted and normalized.
        **attr: Any
            Graph attributes passed to the underlying NetworkX graph.

        Raises
        ------
        ValueError
            If an input edge is missing `sign`, contains an invalid sign, or
            duplicates an existing signed edge.
        """

        super().__init__(**attr)

        if graph is not None:
            self._replace_with_graph(self._plain_graph(graph))
        else:
            self._validate_graph()

    def __str__(self) -> str:
        """
        Return a compact human-readable representation of the influence graph.

        Returns
        -------
        str
            Compact graph representation displaying the number of nodes and edges.
        """

        return (
            f"{type(self).__name__}"
            f"(nodes={self.number_of_nodes()}, "
            f"edges={self.number_of_edges()})"
        )

    __repr__ = __str__

    def copy(self, as_view: bool = False) -> "InfluenceGraph":
        """
        Return a copy of the influence graph.

        Returns
        -------
        InfluenceGraph
            Copy of the influence graph.
        """

        if as_view:
            raise NotImplementedError("InfluenceGraph does not support view copies.")

        return type(self)(self._plain_graph(self))

    def to_directed(self, as_view: bool = False) -> "InfluenceGraph":
        """
        Return a directed copy of the influence graph.

        Returns
        -------
        InfluenceGraph
            Directed copy of the influence graph.

        Notes
        -----
        InfluenceGraph objects are already directed.
        """

        if as_view:
            raise NotImplementedError("InfluenceGraph does not support directed views.")

        return self.copy()

    def to_undirected(self, *args: Any, **kwargs: Any) -> NoReturn:
        """
        Disable conversion to undirected graphs.

        InfluenceGraph relies on directed signed influences and therefore does not
        support conversion to undirected graphs.

        Raises
        ------
        NotImplementedError
            Always raised.
        """

        raise NotImplementedError(
            "InfluenceGraph does not support conversion to undirected graphs."
        )

    # InfluenceGraph intentionally exposes source/target names instead of
    # NetworkX's u_for_edge/v_for_edge names.
    def add_edge(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        source: Any,
        target: Any,
        key: Any = None,
        sign: Optional[InfluenceSign] = None,
        **attr: Any,
    ) -> int:
        """
        Add a signed influence edge between two nodes.

        Nodes are automatically added to the graph if they do not already exist.
        The edge sign must be passed explicitly through `sign` and is normalized to
        -1 or 1.

        At most one positive edge and one negative edge may exist between a pair of
        nodes.

        Examples
        --------
        >>> ig = InfluenceGraph()
        >>> ig.add_edge("A", "B", sign=1)

        >>> sorted(ig.edges(data="sign"))
        [('A', 'B', 1)]

        >>> ig.add_edge("A", "B", sign=-1)
        >>> sorted(ig.edges(data="sign"))
        [('A', 'B', 1), ('A', 'B', -1)]

        Parameters
        ----------
        source, target: Any
            Source and target nodes.
        sign: {-1, 1}
            Sign of the influence.
        **attr: Any
            Additional edge attributes.

        Raises
        ------
        ValueError
            If `sign` is invalid or duplicates an existing signed edge.
        """

        if sign is None:
            sign = key
            key = None

        if sign is None:
            raise TypeError("missing required argument: 'sign'")

        sign = cast(CircuitSign, self._normalize_sign(sign))

        if self.has_edge(source, target):
            signs = self._edge_signs(source, target)

            if sign in signs:
                raise ValueError(
                    f"duplicated edge sign for edge "
                    f"{source!r} -> {target!r} with sign {sign!r}"
                )

        attr["sign"] = sign

        return super().add_edge(
            source,
            target,
            key=key,
            **attr,
        )

    def add_edges_from(
        self,
        ebunch_to_add: Iterable[Any],
        **attr: Any,
    ) -> None:
        """
        Add multiple signed edges.

        Each added edge must contain a valid `sign` attribute, either directly in
        the edge data mapping or through shared keyword attributes. Edge-specific
        attributes take precedence over shared keyword attributes.

        The operation is transactional: if one edge is invalid, the graph is left
        unchanged.

        Examples
        --------
        >>> ig = InfluenceGraph()
        >>> ig.add_edges_from(
        ...     [
        ...         ("A", "B", {"sign": 1}),
        ...         ("B", "C", {"sign": -1}),
        ...     ]
        ... )

        >>> sorted(ig.edges(data="sign"))
        [('A', 'B', 1), ('B', 'C', -1)]

        Shared edge attributes may also be provided:

        >>> ig.add_edges_from(
        ...     [
        ...         ("C", "D"),
        ...         ("D", "E"),
        ...     ],
        ...     sign=1,
        ... )
        >>> sorted(ig.edges(data="sign"))
        [('A', 'B', 1), ('B', 'C', -1), ('C', 'D', 1), ('D', 'E', 1)]

        Invalid edges are rejected without modifying the graph:

        >>> ig.add_edges_from([("E", "F")])
        Traceback (most recent call last):
            ...
        ValueError: missing edge attribute 'sign' for edge 'E' -> 'F'

        Parameters
        ----------
        ebunch_to_add: Iterable
            Iterable of edges to add. Edges may be provided as:
                - `(source, target, {"sign": sign})`;
                - `(source, target, key, {"sign": sign})`.
            Each edge must contain a valid `sign` attribute, unless `sign` is
            provided as a shared keyword attribute.
        **attr: Any
            Shared edge attributes added to all edges. Edge-specific attributes
            take precedence over these shared attributes.

        Raises
        ------
        ValueError
            If one added edge has no `sign` attribute or an invalid sign value.

        Notes
        -----
        Edge attributes specified directly in `ebunch_to_add` take precedence over
        shared keyword attributes.

        When adding edges from an iterator over the same graph being modified, a
        `RuntimeError` may occur because the graph changes during iteration. To
        avoid this issue, materialize the iterator first, for example with
        `list(...)`.
        """

        candidate = nx.MultiDiGraph()
        nx.MultiDiGraph.add_edges_from(candidate, ebunch_to_add, **attr)
        candidate = self._validated_update_graph(candidate)

        self._check_update_graph(candidate)
        self._apply_update_graph(candidate)

    def add_weighted_edges_from(self, ebunch_to_add, weight="weight", **attr) -> None:
        """
        Disable weighted edge insertion.

        Raises
        ------
        NotImplementedError
            Always raised.
        """
        raise NotImplementedError(
            "InfluenceGraph does not support add_weighted_edges_from(); "
            "use add_edges_from(..., sign=...) instead."
        )

    def update(
        self,
        edges: Any = None,
        nodes: Any = None,
    ) -> None:
        """
        Update the influence graph while preserving signed-edge invariants.

        This method adds nodes and/or edges from another graph-like object or from
        explicit node and edge collections. All added edges must contain a valid
        `sign` attribute. The update is transactional: if one edge is missing a
        valid sign, the current graph is left unchanged.

        Valid edge signs are normalized to -1 or 1.

        Examples
        --------
        Update from another signed influence graph:

            G1:
                A → B

            G2:
                B ─| C
                C → A

        >>> g1 = InfluenceGraph()
        >>> g1.add_edge("A", "B", sign=1)

        >>> g2 = InfluenceGraph()
        >>> g2.add_edge("B", "C", sign=-1)
        >>> g2.add_edge("C", "A", sign=1)

        >>> g1.update(g2)
        >>> sorted(g1.edges(data="sign"))
        [('A', 'B', 1), ('B', 'C', -1), ('C', 'A', 1)]

        Update from an adjacency-like dictionary by explicitly converting it into
        signed edges:

            adjacency:
                A → B
                A → C
                B → C

        >>> adjacency = {
        ...     "A": ["B", "C"],
        ...     "B": ["C"],
        ... }
        >>> edges = [
        ...     (source, target, {"sign": 1})
        ...     for source, targets in adjacency.items()
        ...     for target in targets
        ... ]

        >>> ig = InfluenceGraph()
        >>> ig.update(edges=edges, nodes=adjacency)
        >>> sorted(ig.edges(data="sign"))
        [('A', 'B', 1), ('A', 'C', 1), ('B', 'C', 1)]

        Invalid updates are rejected without modifying the graph:

        >>> ig = InfluenceGraph()
        >>> ig.add_edge("A", "B", sign=1)
        >>> ig.update(edges=[("B", "C")])
        Traceback (most recent call last):
            ...
        ValueError: missing edge attribute 'sign' for edge 'B' -> 'C'

        Parameters
        ----------
        edges: graph-like object or collection of edges, optional
            Graph-like object or edge collection used to update the influence graph.
            If `edges` has `nodes` and `edges` attributes, it is interpreted as a
            graph-like object. Otherwise, it is interpreted as an edge collection.

            Edge collections may contain:
                - `(source, target, data)` triples;
                - `(source, target, key, data)` quadruples for multiedges.

            Each edge data mapping must contain a valid `sign` attribute.
        nodes: collection of nodes, optional
            Nodes to add. Ignored when `edges` is a graph-like object.

        Raises
        ------
        ValueError
            If an added edge has no `sign` attribute or an invalid sign value.
        """

        candidate = nx.MultiDiGraph()
        nx.MultiDiGraph.update(candidate, edges=edges, nodes=nodes)
        candidate = self._validated_update_graph(candidate)

        self._check_update_graph(candidate)
        self._apply_update_graph(candidate)

    def rename(self, old: str, new: str) -> None:
        """
        Rename a node and merge signed edges created by the rename.

        The operation is transactional. If `old` is renamed to an existing
        node, both nodes are fused. Signed edges that become exact duplicates
        after the rename are collapsed, while opposite signs are preserved as
        distinct influences.

        Examples
        --------
        >>> graph = InfluenceGraph()
        >>> graph.add_edge("A", "B", sign=1)
        >>> graph.add_edge("A", "C", sign=1)
        >>> graph.add_edge("C", "D", sign=1)
        >>> graph.rename("C", "B")
        >>> sorted(graph.edges(data="sign"))
        [('A', 'B', 1), ('B', 'D', 1)]

        Parameters
        ----------
        old: str
            Node to rename.
        new: str
            New node name.

        Raises
        ------
        TypeError
            If `old` or `new` is not a string.
        KeyError
            If `old` is not present in the graph.
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

        if old == new:
            return None

        if old not in self:
            raise KeyError(f"node {old!r} not found")

        self.relabel({old: new})

    def relabel(self, mapping: Mapping[str, str]) -> None:
        """
        Relabel several nodes in one operation.

        Nodes absent from the graph are ignored. If relabeling fuses nodes,
        signed edges that become exact duplicates are collapsed, while opposite
        signs are preserved as distinct influences.

        Examples
        --------
        >>> graph = InfluenceGraph()
        >>> graph.add_edge("Trp53", "Sox2", sign=1)
        >>> graph.relabel({"Trp53": "TP53", "Sox2": "SOX2"})
        >>> sorted(graph.edges(data="sign"))
        [('TP53', 'SOX2', 1)]

        Parameters
        ----------
        mapping: Mapping[str, str]
            Node rename mapping.

        Raises
        ------
        TypeError
            If `mapping` is not a mapping from strings to strings.
        """

        if not isinstance(mapping, MappingInstance):
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

        plain_graph = self._plain_graph(self)
        relabeled_graph = nx.relabel_nodes(
            plain_graph,
            mapping=active_mapping,
            copy=True,
        )

        relabeled = self._empty_like()
        relabeled.graph.update(relabeled_graph.graph)

        for node, data in relabeled_graph.nodes(data=True):
            relabeled.add_node(node, **data.copy())

        for source, target, data in relabeled_graph.edges(data=True):
            edge_data = data.copy()
            sign = edge_data.pop("sign")

            try:
                relabeled.add_edge(source, target, sign=sign, **edge_data)
            except ValueError as error:
                if "duplicated edge sign" not in str(error):
                    raise

        self._replace_with_graph(relabeled)

    def strongly_connected_components(
        self,
        include_singleton_selfloops: bool = True,
    ) -> List[Set[str]]:
        """
        Return strongly connected components involved in feedback structures.

        A strongly connected component (SCC) is a maximal set of nodes such that every
        node is reachable from every other node through directed paths.

        By default, singleton SCCs are excluded unless the corresponding node contains
        a self-loop, since isolated nodes without autoregulation do not participate in
        feedback circuits.

        Examples
        --------
        Consider the following influence graph:

            A → B → C → D → E
            ↑       ↓   ↑   ↓
            └───────┘   └───┘

            F → G

            H ↺

        The graph contains three feedback SCCs:
            - {A, B, C}
            - {D, E}
            - {H}

        >>> ig.strongly_connected_components()
        [{'A', 'B', 'C'}, {'D', 'E'}, {'H'}]

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should be
            included.

        Returns
        -------
        list of set of str
            List of strongly connected components represented as node sets.
        """
        sccs = []

        for scc in nx.strongly_connected_components(self):
            if len(scc) > 1:
                sccs.append(set(scc))

            elif include_singleton_selfloops:
                node = next(iter(scc))
                if self.has_edge(node, node):
                    sccs.append(set(scc))

        return sccs

    def feedback_nodes(
        self,
        include_singleton_selfloops: bool = True,
    ) -> Set[str]:
        """
        Return nodes belonging to feedback-relevant strongly connected components.

        Feedback nodes correspond to nodes participating in directed feedback
        structures, i.e. non-singleton strongly connected components or singleton
        components containing a self-loop.

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should be
            included.

        Returns
        -------
        set of str
            Set of nodes participating in feedback structures.
        """

        nodes = set()

        for scc in self.strongly_connected_components(
            include_singleton_selfloops=include_singleton_selfloops,
        ):
            nodes.update(scc)

        return nodes

    def regulators(self) -> List[Any]:
        """
        Return regulator-only nodes.

        Regulator nodes are nodes with outgoing influences and no incoming
        influences.

        Examples
        --------
        Consider the following influence graph:

            A → B → C
            D ─| B
            E

        Nodes `A` and `D` regulate other nodes but are not themselves regulated.

        >>> ig.regulators()
        ['A', 'D']

        Isolated nodes are excluded:

        >>> "E" in ig.regulators()
        False

        Returns
        -------
        list
            Regulator-only nodes.
        """

        return [
            node
            for node in self.nodes()
            if self.out_degree(node) > 0 and self.in_degree(node) == 0
        ]

    def targets(self) -> List[Any]:
        """
        Return terminal target nodes.

        Target nodes are nodes with incoming influences and no outgoing
        influences.

        Examples
        --------
        Consider the following influence graph:

            A → B → C
            D ─| B
            E

        Node `C` is regulated but does not regulate any other node.

        >>> ig.targets()
        ['C']

        Isolated nodes are excluded:

        >>> "E" in ig.targets()
        False

        Returns
        -------
        list
            Terminal target nodes.
        """

        return [
            node
            for node in self.nodes()
            if self.in_degree(node) > 0 and self.out_degree(node) == 0
        ]

    def structural_families(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
    ) -> Dict[StructuralSignature, Set[str]]:
        """
        Return families of structurally equivalent nodes.

        Two nodes are considered structurally equivalent if they share identical
        signed predecessor signatures and, optionally, identical signed successor
        signatures.

        Structural signatures preserve edge signs. Therefore, nodes influenced by
        the same regulators but with different signs are assigned to different
        families.

        By default, nodes involved in feedback structures are excluded from family
        detection, since feedback circuits often carry important dynamical
            information that should not be collapsed.

        Examples
        --------
        Consider the following influence graph:

            TF1 → g1
            TF1 → g2
            TF1 → g3
            TF1 → g4

            TF2 ─| g3
            TF2 ─| g4

        where:
            - `→` denotes positive regulation,
            - `─|` denotes negative regulation.

        Nodes `g1` and `g2` belong to the same structural family because they share
        identical signed regulators:

        >>> ig.structural_families()
        {
            (
                frozenset({('TF1', 1)}),
                frozenset()
            ): {'g1', 'g2'},

            (
                frozenset({('TF1', 1), ('TF2', -1)}),
                frozenset()
            ): {'g3', 'g4'}
        }

        Parameters
        ----------
        include_successors: bool (default: True)
            Whether signed successor signatures should also be included in the
            structural equivalence criterion.
        exclude_feedback_nodes: bool (default: True)
            Whether nodes participating in feedback structures should be excluded.
        min_size: int (default: 2)
            Minimum family size required to keep a structural family.

        Returns
        -------
        dict
            Mapping from structural signatures to node families.
        """

        signatures = {}

        feedback_nodes = self.feedback_nodes() if exclude_feedback_nodes else set()

        for node in self.nodes():

            if node in feedback_nodes:
                continue

            predecessors = frozenset(
                (
                    predecessor,
                    self._normalize_sign(data.get("sign", 1)),
                )
                for predecessor, _, data in self.in_edges(
                    node,
                    data=True,
                )
            )

            if include_successors:

                successors = frozenset(
                    (
                        successor,
                        self._normalize_sign(data.get("sign", 1)),
                    )
                    for _, successor, data in self.out_edges(
                        node,
                        data=True,
                    )
                )

            else:
                successors = frozenset()

            signature = (predecessors, successors)

            if signature not in signatures:
                signatures[signature] = set()

            signatures[signature].add(node)

        return {
            signature: family
            for signature, family in signatures.items()
            if len(family) >= min_size
        }

    def family_collapsed_graph(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
        sep: str = "|",
    ) -> "InfluenceGraph":
        """
        Return a graph where structurally equivalent nodes are collapsed.

        Nodes belonging to the same structural family are replaced by a single
        composite node. The composite node name is obtained by joining member names
        with `sep`, and the original members are stored in the node attribute
        `members`.

        Edges incident to collapsed nodes are rewired to the corresponding composite
        node. Since structural families are defined from signed predecessor and
        successor signatures, this collapse preserves signed influence structure.

        Examples
        --------
        Consider the following influence graph:

            TF1 → g1
            TF1 → g2
            TF1 → g3
            TF1 → g4

            TF2 ─| g3
            TF2 ─| g4

        >>> collapsed = ig.family_collapsed_graph()
        >>> sorted(collapsed.nodes())
        ['TF1', 'TF2', 'g1|g2', 'g3|g4']

        The composite nodes store their original members:

        >>> collapsed.nodes["g1|g2"]["members"]
        {'g1', 'g2'}

        Parameters
        ----------
        include_successors: bool (default: True)
            Whether signed successor signatures should also be included in the
            structural equivalence criterion.
        exclude_feedback_nodes: bool (default: True)
            Whether nodes participating in feedback structures should be excluded
            from collapsing.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.

        Returns
        -------
        InfluenceGraph
            Family-collapsed influence graph.
        """

        families = self.structural_families(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
        )

        node_mapping = {}

        for family in families.values():
            family_name = sep.join(sorted(family))

            for node in family:
                node_mapping[node] = family_name

        collapsed = type(self)()

        for node, data in self.nodes(data=True):
            collapsed_node = node_mapping.get(node, node)

            if collapsed_node not in collapsed:
                collapsed.add_node(collapsed_node)

            if node in node_mapping:
                collapsed.nodes[collapsed_node].setdefault("members", set())
                collapsed.nodes[collapsed_node]["members"].add(node)

            else:
                collapsed.nodes[collapsed_node].update(data)
                collapsed.nodes[collapsed_node].setdefault("members", {node})

        for source, target, data in self.edges(data=True):
            collapsed_source = node_mapping.get(source, source)
            collapsed_target = node_mapping.get(target, target)

            if collapsed_source == collapsed_target:
                continue

            sign = self._normalize_sign(data["sign"])

            if collapsed.has_edge(
                collapsed_source, collapsed_target
            ) and sign in collapsed._edge_signs(
                collapsed_source,
                collapsed_target,
            ):
                continue

            collapsed.add_edge(
                collapsed_source,
                collapsed_target,
                **data,
            )

        return collapsed

    def family_compressed_graph(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
        sep: str = "|",
    ) -> "InfluenceGraph":
        """
        Deprecated alias for `family_collapsed_graph()`.
        """

        if type(self) is InfluenceGraph:
            _warn_deprecated(
                "family_compressed_graph(...)",
                replacement="family_collapsed_graph(...)",
                stacklevel=2,
            )

        return self.family_collapsed_graph(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            sep=sep,
        )

    def edge_sign(
        self,
        source: str,
        target: str,
    ) -> int:
        """
        Return the aggregated sign between two nodes.

        If all edges between `source` and `target` are positive, the returned sign
        is 1. If all are negative, the returned sign is -1. If both positive and
        negative influences coexist, the returned sign is 0.

        Examples
        --------
        Consider the following influence graph:

            A → B
            A ─| B

            C → D

        >>> ig.edge_sign("A", "B")
        0

        >>> ig.edge_sign("C", "D")
        1

        Parameters
        ----------
        source: str
            Source node.
        target: str
            Target node.

        Returns
        -------
        int
            Aggregated edge sign:
                - 1 for positive influence,
                - -1 for negative influence,
                - 0 for ambiguous bi-signed influence.

        Raises
        ------
        KeyError
            If no edge exists between `source` and `target`.
        """

        if not self.has_edge(source, target):
            raise KeyError(f"no edge found between {source!r} and {target!r}")

        signs = set(self._edge_signs(source, target))

        if signs == {1}:
            return 1

        if signs == {-1}:
            return -1

        return 0

    def autoregulations(
        self,
        sign: Optional[InfluenceSign] = None,
    ) -> List[Tuple[str, int]]:
        """
        Return signed autoregulations.

        Autoregulations correspond to self-loops in the influence graph. Edge signs are
        normalized to -1 or 1.

        Examples
        --------
        Consider the following influence graph:

            A ↺
            B ─| B
            C → D

        >>> ig.autoregulations()
        [('A', 1), ('B', -1)]

        >>> ig.autoregulations(sign=1)
        [('A', 1)]

        Parameters
        ----------
        sign: {-1, 1}, optional
            If provided, return only autoregulations with the requested sign.

        Returns
        -------
        list of tuple
            List of `(node, sign)` autoregulations.
        """

        loops = []

        for node in nx.nodes_with_selfloops(self):
            for edge_data in self.get_edge_data(node, node).values():
                edge_sign = self._normalize_sign(edge_data.get("sign", 1))

                if sign is None or edge_sign == sign:
                    loops.append((node, edge_sign))

        return loops

    def path_sign(
        self,
        path: Iterable[str],
    ) -> int:
        """
        Return the aggregated sign of a directed path.

        The path sign corresponds to the product of aggregated edge signs along the
        path:
            - positive paths have sign +1,
            - negative paths have sign -1.

        If one of the traversed edges is bi-signed, the returned path sign is 0.

        Examples
        --------
        Consider the following influence graph:

            A → B ─| C

            D → E
            D ─| E

        >>> ig.path_sign(["A", "B", "C"])
        -1

        >>> ig.path_sign(["D", "E"])
        0

        Parameters
        ----------
        path: Iterable[str]
            Ordered sequence of nodes describing a directed path.

        Returns
        -------
        int
            Aggregated path sign:
                - 1 for positive effect,
                - -1 for negative effect,
                - 0 for ambiguous bi-signed effect.

        Raises
        ------
        KeyError
            If one of the traversed edges does not exist.
        ValueError
            If the path contains fewer than two nodes.
        """

        path = list(path)

        if len(path) < 2:
            raise ValueError("path must contain at least two nodes")

        sign = 1

        for source, target in zip(path, path[1:]):

            edge_sign = self.edge_sign(
                source,
                target,
            )

            if edge_sign == 0:
                return 0

            sign *= edge_sign

        return sign

    @overload
    def signed_path_string(
        self,
        *nodes: str,
    ) -> str: ...

    @overload
    def signed_path_string(
        self,
        *nodes: Union[str, List[str], Tuple[str, ...]],
    ) -> str: ...

    def signed_path_string(
        self,
        *nodes: Union[str, List[str], Tuple[str, ...]],
    ) -> str:
        """
        Return a human-readable signed representation of a directed path.

        Positive influences are represented with `->`, negative influences with
        `-|`, and bi-signed influences with `--`.

        Examples
        --------
        Consider the following influence graph:

            A → B
            B ─| C
            C → D

        >>> ig.signed_path_string("A", "B", "C", "D")
        'A -> B -| C -> D'

        Parameters
        ----------
        *nodes: str
            Nodes describing a directed path in traversal order.

        Returns
        -------
        str
            Human-readable signed path representation.

        Raises
        ------
        ValueError
            If the path contains fewer than two nodes.
        KeyError
            If one traversed edge does not exist.
        """

        path = nodes

        if len(nodes) == 1 and isinstance(nodes[0], (list, tuple)):
            path = cast(Tuple[str, ...], tuple(nodes[0]))
        else:
            path = cast(Tuple[str, ...], nodes)

        if len(path) < 2:
            raise ValueError("path must contain at least two nodes")

        string = str(path[0])

        for source, target in zip(path, path[1:]):
            sign = self.edge_sign(source, target)

            if sign == 1:
                string += f" -> {target}"

            elif sign == -1:
                string += f" -| {target}"

            else:
                string += f" -- {target}"

        return string

    def marker_paths(
        self,
        markers: Iterable[str],
        direction: Direction = "both",
        sccs: Optional[Iterable[Iterable[str]]] = None,
        sign: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Return shortest paths between markers and feedback SCCs.

        Paths are computed between each marker and each feedback-relevant strongly
        connected component (SCC). For a given SCC/marker pair, all shortest paths
        with minimal length are returned.

        No graph collapse is performed, so alternative minimal routes are
        preserved whenever several SCC exit or entry points exist.

        Examples
        --------
        Consider the following influence graph:

            marker1 → A → B → C → D → E → marker2
                      ↑       ↓   ↑   ↓
                      └───────┘   └───┘

            D → marker2

        The graph contains two feedback SCCs:
            - {A, B, C}
            - {D, E}

        The SCC `{D, E}` reaches `marker2` through two distinct shortest paths:

        >>> ig.marker_paths(["marker2"], direction="downstream")
        [
            {
                "marker": "marker2",
                "scc": frozenset({"A", "B", "C"}),
                "direction": "downstream",
                "paths": [
                    {
                        "path": ["C", "D", "marker2"],
                        "sign": 1,
                    },
                ],
            },
            {
                "marker": "marker2",
                "scc": frozenset({"D", "E"}),
                "direction": "downstream",
                "paths": [
                    {
                        "path": ["D", "marker2"],
                        "sign": 1,
                    },
                    {
                        "path": ["E", "marker2"],
                        "sign": 1,
                    },
                ],
            },
        ]

        Parameters
        ----------
        markers: Iterable[str]
            Marker nodes of interest.
        direction: {"upstream", "downstream", "both"} (default: "both")
            Direction of marker/SCC relationships.
        sccs: Iterable[Iterable[str]], optional
            SCCs to use. If None, feedback SCCs are computed.
        sign: {-1, 0, 1}, optional
            If provided, return only paths with the requested aggregated sign.

        Returns
        -------
        list of dict
            SCC/marker path records with keys:
                - `marker`,
                - `scc`,
                - `direction`,
                - `paths`.

            Each `paths` entry contains:
                - `path`: shortest directed path,
                - `sign`: aggregated path sign.
        """

        if direction not in {"upstream", "downstream", "both"}:
            raise ValueError(
                "unsupported direction: expected 'upstream', 'downstream' "
                f"or 'both', but received {direction!r}"
            )

        if sccs is None:
            sccs = self.strongly_connected_components()

        sccs = list(sccs)

        records = []

        for marker in markers:
            if marker not in self:
                continue

            directions = (
                ["downstream", "upstream"] if direction == "both" else [direction]
            )
            distances_by_direction = {}

            for current_direction in directions:
                if current_direction == "downstream":
                    search_graph = self.reverse(copy=False)
                else:
                    search_graph = self

                distances_by_direction[current_direction] = (
                    nx.single_source_shortest_path_length(search_graph, marker)
                )

            for scc in sccs:
                scc = frozenset(scc)

                for current_direction in directions:
                    candidate_paths = []
                    distances = distances_by_direction[current_direction]
                    reachable_nodes = [node for node in scc if node in distances]

                    if not reachable_nodes:
                        continue

                    min_distance = min(distances[node] for node in reachable_nodes)

                    for node in reachable_nodes:
                        if distances[node] != min_distance:
                            continue

                        if current_direction == "downstream":
                            paths = nx.all_shortest_paths(
                                self,
                                source=node,
                                target=marker,
                            )

                        else:
                            paths = nx.all_shortest_paths(
                                self,
                                source=marker,
                                target=node,
                            )

                        for path in paths:

                            path_sign = self.path_sign(path)

                            if sign is None or path_sign == sign:
                                candidate_paths.append(
                                    {
                                        "path": path,
                                        "sign": path_sign,
                                    }
                                )

                    if candidate_paths:

                        min_length = min(
                            len(record["path"]) for record in candidate_paths
                        )

                        shortest_paths = [
                            record
                            for record in candidate_paths
                            if len(record["path"]) == min_length
                        ]

                        records.append(
                            {
                                "marker": marker,
                                "scc": scc,
                                "direction": current_direction,
                                "paths": shortest_paths,
                            }
                        )

        return records

    def circuits(
        self,
        sign: Optional[InfluenceSign] = None,
    ) -> List[Tuple[List[str], int]]:
        """
        Return signed feedback circuits.

        Circuits are detected as simple directed cycles. For multigraph edges,
        all possible edge-sign combinations are considered.

        The sign of a circuit corresponds to the product of edge signs along the
        cycle:
            - positive circuits have sign +1,
            - negative circuits have sign -1.

        Examples
        --------
        Consider the following influence graph:

            A → B ─| C
            ↑       ↓
            └───────┘

            D ─| E
            ↑   ↓
            └───┘

        The graph contains:
            - one negative circuit: A → B ─| C → A,
            - one positive circuit: D ─| E ─| D.

        >>> ig.circuits()
        [
            (['A', 'B', 'C'], -1),
            (['D', 'E'], 1),
        ]

        >>> ig.circuits(sign=1)
        [
            (['D', 'E'], 1),
        ]

        >>> ig.circuits(sign=-1)
        [
            (['A', 'B', 'C'], -1),
        ]

        Parameters
        ----------
        sign: {-1, 1}, optional
            If provided, return only circuits with the requested sign.

        Returns
        -------
        list of tuple
            Each entry is `(cycle, sign)`, where `cycle` is a list of nodes and
            `sign` is the product of edge signs along the circuit.
        """

        simple_graph = nx.DiGraph()
        simple_graph.add_nodes_from(self.nodes())
        simple_graph.add_edges_from(self.edges())

        circuits = []

        for cycle in nx.simple_cycles(simple_graph):
            edge_signs = []

            for source, target in self._cycle_edges(cycle):
                signs = self._edge_signs(source, target)

                if not signs:
                    break

                edge_signs.append(signs)

            else:
                for signs in product(*edge_signs):
                    circuit_sign = 1

                    for edge_sign in signs:
                        circuit_sign *= edge_sign

                    if sign is None or circuit_sign == sign:
                        circuits.append((cycle, circuit_sign))

        return circuits

    def positive_circuits(self) -> List[List[str]]:
        """
        Return positive feedback circuits.

        Returns
        -------
        list of list of str
            Positive feedback circuits.
        """

        return [cycle for cycle, _ in self.circuits(sign=1)]

    def negative_circuits(self) -> List[List[str]]:
        """
        Return negative feedback circuits.

        Returns
        -------
        list of list of str
            Negative feedback circuits.
        """

        return [cycle for cycle, _ in self.circuits(sign=-1)]

    def feedback_induced_graph(
        self,
        include_singleton_selfloops: bool = True,
    ) -> "InfluenceGraph":
        """
        Return the graph induced by feedback nodes.

        Feedback nodes are nodes participating in feedback structures, i.e.
        non-singleton strongly connected components or singleton components
        containing a self-loop. The returned graph is the subgraph induced by these
        nodes, so edges between feedback nodes are preserved.

        Examples
        --------
        Consider the following influence graph:

            A → B → C → D → E
            ↑       ↓   ↑   ↓
            └───────┘   └───┘

            C → X → Y

            F → G

            H ↺

        The graph contains three feedback structures:
            - {A, B, C}
            - {D, E}
            - {H}

        >>> feedback = ig.feedback_induced_graph()

        >>> sorted(feedback.nodes())
        ['A', 'B', 'C', 'D', 'E', 'H']

        >>> sorted(feedback.edges())
        [
            ('A', 'B'),
            ('B', 'C'),
            ('C', 'A'),
            ('C', 'D'),
            ('D', 'E'),
            ('E', 'D'),
            ('H', 'H'),
        ]

        Nodes `X`, `Y`, `F` and `G` are excluded because they do not participate in
        feedback structures.

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should be
            included.

        Returns
        -------
        InfluenceGraph
            Influence graph induced by feedback nodes.
        """

        feedback_nodes = self.feedback_nodes(
            include_singleton_selfloops=include_singleton_selfloops,
        )

        return type(self)(
            cast(_MultiDiGraphBase, self.subgraph(feedback_nodes).copy()),
        )

    def collapsed_graph(
        self,
        include_singleton_selfloops: bool = True,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = False,
        min_size: int = 2,
        sep: str = "|",
    ) -> "InfluenceGraph":
        """
        Return a collapsed influence graph.

        Collapse is applied sequentially through:
            1. feedback-induced graph extraction,
            2. structural family collapse.

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should be
            included in the feedback-induced graph.
        include_successors: bool (default: True)
            Whether signed successor signatures should be included in structural
            family detection.
        exclude_feedback_nodes: bool (default: False)
            Whether nodes participating in feedback structures should be excluded
            from structural family collapse.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.

        Returns
        -------
        InfluenceGraph
            Collapsed influence graph.

        See Also
        --------
        feedback_induced_graph
        family_collapsed_graph
        """

        graph = self.feedback_induced_graph(
            include_singleton_selfloops=include_singleton_selfloops,
        )

        return graph.family_collapsed_graph(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            sep=sep,
        )

    def compressed_graph(
        self,
        include_singleton_selfloops: bool = True,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = False,
        min_size: int = 2,
        sep: str = "|",
    ) -> "InfluenceGraph":
        """
        Deprecated alias for `collapsed_graph()`.
        """

        if type(self) is InfluenceGraph:
            _warn_deprecated(
                "compressed_graph(...)",
                replacement="collapsed_graph(...)",
                stacklevel=2,
            )

        return self.collapsed_graph(
            include_singleton_selfloops=include_singleton_selfloops,
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            sep=sep,
        )

    def to_graphviz(
        self,
        program: str = "dot",
        edge_style: Optional[Callable[[Mapping[str, Any]], Mapping[str, Any]]] = None,
        **kwargs: Any,
    ):
        """
        Convert the influence graph to a native graphviz Digraph.

        Positive influences are displayed as green activating edges, while
        negative influences are displayed as red inhibitory edges. This method
        uses the `graphviz` Python package directly and does not depend on
        pydot.

        Parameters
        ----------
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting graph.
        edge_style: Callable, optional
            Optional callable used to update edge attributes. The callable
            receives edge attribute dictionaries and must return a mapping of
            graphviz edge attributes.
        **kwargs: Any
            Graph attributes assigned to the resulting graphviz object.

        Returns
        -------
        graphviz.Digraph
            Native graphviz influence graph.

        Raises
        ------
        ImportError
            If the `graphviz` Python package is not installed.
        """

        graph = self.copy()

        for _, _, edge_data in graph.edges(data=True):
            sign = self._normalize_sign(edge_data.get("sign", 1))

            edge_data.update(
                color="green4" if sign == 1 else "red2",
                arrowhead="normal" if sign == 1 else "tee",
                penwidth=2,
            )

            if edge_style is not None:
                edge_data.update(edge_style(edge_data))

        return _networkx_to_graphviz(graph, program=program, **kwargs)

    def to_pydot(
        self,
        program: str = "dot",
        edge_style: Optional[Callable[[Mapping[str, Any]], Mapping[str, Any]]] = None,
        **kwargs: Any,
    ) -> "Dot":
        """
        Convert the influence graph to a pydot graph.

        Positive influences are displayed as green activating edges, while negative
        influences are displayed as red inhibitory edges.

        Parameters
        ----------
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting pydot graph.
        edge_style: Callable, optional
            Optional callable used to update edge attributes. The callable receives
            edge attribute dictionaries and must return a mapping of pydot edge
            attributes.
        **kwargs: Any
            Keyword arguments passed to the resulting pydot graph using
            `dot.set(key, value)`.

        Returns
        -------
        Dot
            Pydot influence graph.
        """

        graph = self.copy()

        for _, _, edge_data in graph.edges(data=True):

            sign = self._normalize_sign(
                edge_data.get("sign", 1),
            )

            edge_data.update(
                color="green4" if sign == 1 else "red2",
                arrowhead="normal" if sign == 1 else "tee",
                penwidth=2,
            )

            if edge_style is not None:
                edge_data.update(edge_style(edge_data))

        dot = nx.drawing.nx_pydot.to_pydot(graph)

        dot.set_prog(program)

        for key, value in kwargs.items():
            dot.set(key, value)

        return dot

    def show(
        self,
        program: str = "dot",
        edge_style: Optional[Callable[[Mapping[str, Any]], Mapping[str, Any]]] = None,
        width: Optional[SvgLength] = None,
        height: Optional[SvgLength] = None,
        **kwargs: Any,
    ) -> None:
        """
        Display the influence graph in a Jupyter/IPython environment.

        The graph is rendered through Graphviz using the `to_pydot()` method and
        displayed as an SVG image.

        Parameters
        ----------
        program: str (default: "dot")
            Graphviz layout program used for rendering.
        edge_style: Callable, optional
            Optional callable used to update edge attributes before rendering.
        width: str or int or float, optional
            Display width assigned to the rendered SVG.
        height: str or int or float, optional
            Display height assigned to the rendered SVG.
        **kwargs: Any
            Keyword arguments passed to the underlying pydot graph through
            `dot.set(key, value)`.

        Raises
        ------
        RuntimeError
            If IPython is not available.
        """

        try:
            from IPython.display import SVG, display

        except ImportError:
            raise RuntimeError("show() requires an IPython/Jupyter environment.")

        dot = cast(
            Any,
            self.to_pydot(
                program=program,
                edge_style=edge_style,
                **kwargs,
            ),
        )
        svg = scale_svg(dot.create_svg().decode(), width=width, height=height)

        display(SVG(svg))

    def _empty_like(self) -> "InfluenceGraph":
        """
        Return an empty graph preserving the concrete influence graph type.
        """

        return type(self)()

    def _validated_update_graph(self, graph: Any) -> "InfluenceGraph":
        """
        Return a validated graph containing candidate update data only.
        """

        return type(self)(graph)

    def _check_update_graph(self, graph: "InfluenceGraph") -> None:
        """
        Ensure candidate update edges preserve current signed-edge invariants.
        """

        for source, target, data in graph.edges(data=True):
            sign = self._normalize_sign(data["sign"])

            if self.has_edge(source, target) and sign in self._edge_signs(
                source,
                target,
            ):
                raise ValueError(
                    f"duplicated edge sign for edge "
                    f"{source!r} -> {target!r} with sign {sign!r}"
                )

    def _apply_update_graph(self, graph: "InfluenceGraph") -> None:
        """
        Apply a prevalidated update graph.
        """

        self.graph.update(graph.graph)
        nx.MultiDiGraph.add_nodes_from(
            self,
            ((node, data.copy()) for node, data in graph.nodes(data=True)),
        )

        for source, target, data in graph.edges(data=True):
            edge_data = data.copy()
            sign = edge_data.pop("sign")
            self.add_edge(source, target, sign=sign, **edge_data)

    @staticmethod
    def _plain_graph(graph: Any) -> nx.MultiDiGraph[Any]:
        """
        Return a plain NetworkX MultiDiGraph copy.
        """

        source_graph = nx.MultiDiGraph(graph)
        plain_graph = nx.MultiDiGraph()
        plain_graph.graph.update(source_graph.graph)
        plain_graph.add_nodes_from(
            (node, data.copy()) for node, data in source_graph.nodes(data=True)
        )
        for source, target, key, data in source_graph.edges(keys=True, data=True):
            nx.MultiDiGraph.add_edge(
                plain_graph,
                source,
                target,
                key=key,
                **data.copy(),
            )

        return plain_graph

    def _replace_with_graph(self, graph: Any) -> None:
        """
        Replace graph contents while preserving InfluenceGraph invariants.
        """

        plain_graph = self._plain_graph(graph)

        self.clear()
        self.graph.update(plain_graph.graph)
        nx.MultiDiGraph.add_nodes_from(
            self,
            ((node, data.copy()) for node, data in plain_graph.nodes(data=True)),
        )

        for source, target, data in plain_graph.edges(data=True):
            data = data.copy()

            if "sign" not in data:
                raise ValueError(
                    f"missing edge attribute 'sign' for edge "
                    f"{source!r} -> {target!r}"
                )

            sign = data.pop("sign")
            self.add_edge(source, target, sign=sign, **data)

    @staticmethod
    def _cycle_edges(cycle: List[str]) -> List[Tuple[str, str]]:
        """
        Return ordered directed edges from a cycle node list.
        """

        return list(zip(cycle, cycle[1:] + cycle[:1]))

    def _edge_signs(self, source: str, target: str) -> List[int]:
        """
        Return normalized signs for all edges between two nodes.
        """

        edge_data = self.get_edge_data(source, target, default={})

        return [
            self._normalize_sign(data.get("sign", 1)) for data in edge_data.values()
        ]

    @staticmethod
    def _normalize_sign(sign: Any) -> int:
        """
        Convert common sign encodings to -1 or 1.
        """

        if isinstance(sign, bool):
            raise ValueError(
                "unsupported edge sign: expected -1, 1, '+', '-', "
                f"'positive' or 'negative', but received {sign!r}"
            )

        if sign in [1, "+", "positive"]:
            return 1

        if sign in [-1, "-", "negative"]:
            return -1

        raise ValueError(
            "unsupported edge sign: expected -1, 1, '+', '-', "
            f"'positive' or 'negative', but received {sign!r}"
        )

    def _validate_graph(self) -> None:
        """
        Validate influence graph edge attributes.
        """

        seen = set()

        for source, target, data in self.edges(data=True):
            if "sign" not in data:
                raise ValueError(
                    f"missing edge attribute 'sign' for edge "
                    f"{source!r} -> {target!r}"
                )

            data["sign"] = self._normalize_sign(data["sign"])

            key = (source, target, data["sign"])

            if key in seen:
                raise ValueError(
                    f"duplicated edge sign for edge "
                    f"{source!r} -> {target!r} with sign {data['sign']!r}"
                )

            seen.add(key)


class AggregatedInfluenceGraph(InfluenceGraph):
    """
    Influence graph aggregated across several Boolean networks or influence graphs.

    An `AggregatedInfluenceGraph` extends `InfluenceGraph` by requiring each
    edge to store an occurrence count. The count represents the number of
    source graphs in which the signed influence was observed. The total number
    of aggregated graphs is exposed through the `total` property and updated
    through its validated setter.

    Edge frequencies are not stored as edge attributes. They are computed as:

        count / total

    Examples
    --------
    >>> graph = AggregatedInfluenceGraph(total=4)
    >>> graph.add_edge("A", "B", sign=1, count=3)
    >>> graph.add_edge("B", "C", sign=-1, count=1)
    >>> graph.edge_frequency("A", "B")
    0.75

    >>> graph.autoregulations()
    []

    Parameters
    ----------
    graph: DiGraph or MultiDiGraph, optional
        Initial graph. Every edge must define both `sign` and `count`.
    total: int
        Total number of graphs used to construct the aggregated influence graph.

    Raises
    ------
    TypeError
        If `total` is not an integer, or if an edge count is not an integer.
    ValueError
        If `total` is not positive.
        If an edge count is negative or greater than `total`.
    KeyError
        If an edge does not define a required `sign` or `count` attribute.
    """

    __slots__ = ("__total",)

    def __init__(
        self,
        graph: Optional[_MultiDiGraphBase] = None,
        total: int = 1,
        **attr: Any,
    ) -> None:
        self.__total = self._validate_total(total)
        super().__init__(graph=graph, **attr)
        self.validate_counts()

    @property
    def total(self) -> int:
        """
        Return the total number of aggregated source graphs.

        Assignment through this property is validated, so `graph.total = value`
        cannot set a total smaller than an existing edge count.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.total
        4

        >>> graph.total = 6
        >>> graph.total
        6

        Returns
        -------
        int
            Total number of aggregated source graphs.
        """

        return self.__total

    @total.setter
    def total(self, value: int) -> None:
        """
        Set the total number of aggregated source graphs.

        The new total must be a positive integer and must remain greater than or
        equal to every edge count already stored in the graph.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        >>> graph.total = 6
        >>> graph.total
        6

        >>> graph.edge_frequency("A", "B")
        0.5

        Counts cannot exceed the new total:

        >>> graph.total = 2
        Traceback (most recent call last):
            ...
        ValueError: invalid edge count for edge 'A' -> 'B':
            expected value between 0 and 2 but received 3

        Raises
        ------
        TypeError
            If `value` is not an integer.
        ValueError
            If `value` is not positive or if one existing edge count is greater
            than `value`.
        """

        value = self._validate_total(value)

        for source, target, data in self.edges(data=True):
            self._validate_count(source, target, data.get("count"), total=value)

        self.__total = value

    def copy(self, as_view: bool = False) -> "AggregatedInfluenceGraph":
        """
        Return a copy of the aggregated influence graph.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        >>> copied = graph.copy()
        >>> copied.total
        4

        >>> copied.edge_count("A", "B")
        3

        Parameters
        ----------
        as_view: bool (default: False)
            Whether to return a view instead of a copy. Views are not supported.

        Returns
        -------
        AggregatedInfluenceGraph
            Copy of the aggregated influence graph.

        Raises
        ------
        NotImplementedError
            If `as_view` is True.
        """

        if as_view:
            raise NotImplementedError(
                "AggregatedInfluenceGraph does not support view copies."
            )

        return type(self)(graph=self._plain_graph(self), total=self.total)

    def add_edge(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        source: Any,
        target: Any,
        key: Any = None,
        sign: Optional[InfluenceSign] = None,
        count: Optional[int] = None,
        **attr: Any,
    ) -> int:
        """
        Add an aggregated signed influence.

        Nodes are automatically added to the graph if they do not already exist.
        The edge sign follows the same normalization and duplicate rules as
        `InfluenceGraph`: at most one positive edge and one negative edge may
        exist between a source and target pair.

        The edge count must be an integer between 0 and `total`.

        For compatibility with `InfluenceGraph`, the sign may be passed as the
        third positional argument when no explicit `sign` keyword is provided.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        0

        >>> graph.add_edge("A", "B", -1, count=1)
        1

        >>> sorted(
        ...     (source, target, data["sign"], data["count"])
        ...     for source, target, data in graph.edges(data=True)
        ... )
        [('A', 'B', -1, 1), ('A', 'B', 1, 3)]

        Parameters
        ----------
        source, target: Any
            Source and target nodes.
        sign: {-1, 1}
            Influence sign. If omitted, `key` is interpreted as the sign, as in
            `InfluenceGraph.add_edge("A", "B", 1)`.
        count: int
            Number of source graphs in which the signed influence is observed.
            Integer-valued floats such as `1.0` are accepted and normalized.
        **attr: Any
            Additional edge attributes.

        Returns
        -------
        int
            NetworkX multiedge key assigned to the added edge.

        Raises
        ------
        TypeError
            If `count` is missing or has an unsupported type.
        ValueError
            If `sign` is invalid, if the signed edge already exists, if
            `count` is fractional, or if `count` is outside the valid range.
        """

        if sign is None:
            sign = key
            key = None

        if count is None:
            raise TypeError("missing required argument: 'count'")

        count = self._validate_count(source, target, count)

        return super().add_edge(
            source,
            target,
            key=key,
            sign=sign,
            count=count,
            **attr,
        )

    def add_edges_from(
        self,
        ebunch_to_add: Iterable[Any],
        **attr: Any,
    ) -> None:
        """
        Add multiple aggregated signed edges transactionally.

        All added edges must define valid `sign` and `count` attributes, either
        directly in the edge data mapping or through shared keyword attributes.
        If one edge is invalid, the graph is left unchanged.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=3)
        >>> graph.add_edges_from(
        ...     [
        ...         ("A", "B", {"sign": 1, "count": 2}),
        ...         ("B", "C", {"sign": -1, "count": 1}),
        ...     ]
        ... )

        >>> sorted(
        ...     (source, target, data["sign"], data["count"])
        ...     for source, target, data in graph.edges(data=True)
        ... )
        [('A', 'B', 1, 2), ('B', 'C', -1, 1)]

        Invalid edges are rejected without modifying the graph:

        >>> graph.add_edges_from([("C", "D", {"sign": 1, "count": 4})])
        Traceback (most recent call last):
            ...
        ValueError: invalid edge count for edge 'C' -> 'D':
            expected value between 0 and 3 but received 4

        Parameters
        ----------
        ebunch_to_add: Iterable
            Iterable of edges to add. Each edge must contain `sign` and `count`,
            unless these attributes are supplied through shared keyword
            arguments.
        **attr: Any
            Shared edge attributes added to all edges.

        Raises
        ------
        KeyError
            If one added edge has no `count` attribute.
        TypeError
            If one added edge count is not an integer.
        ValueError
            If one edge has no valid sign, duplicates an existing signed edge, or
            has a count outside the valid range.
        """

        candidate = nx.MultiDiGraph()
        nx.MultiDiGraph.add_edges_from(candidate, ebunch_to_add, **attr)
        candidate = self._validated_update_graph(candidate)

        self._check_update_graph(candidate)
        self._apply_update_graph(candidate)

    def update(
        self,
        edges: Any = None,
        nodes: Any = None,
    ) -> None:
        """
        Update the graph while preserving signed-edge and count invariants.

        This method follows `InfluenceGraph.update()` but also requires all
        added edges to define a valid `count`. The update is transactional: if an
        added edge is invalid, the current graph is left unchanged.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=2)
        >>> graph.add_edge("A", "B", sign=1, count=1)

        >>> graph.update(
        ...     edges=[("B", "C", {"sign": -1, "count": 2})],
        ... )
        >>> graph.edge_count("B", "C")
        2

        Parameters
        ----------
        edges: graph-like object or collection of edges, optional
            Graph-like object or edge collection used to update the aggregated
            influence graph. Added edges must contain valid `sign` and `count`
            attributes.
        nodes: collection of nodes, optional
            Nodes to add. Ignored when `edges` is a graph-like object.

        Raises
        ------
        KeyError
            If one added edge has no `count` attribute.
        TypeError
            If one added edge count is not an integer.
        ValueError
            If one added edge has no valid sign, duplicates an existing signed
            edge, or has a count outside the valid range.
        """

        candidate = nx.MultiDiGraph()
        nx.MultiDiGraph.update(candidate, edges=edges, nodes=nodes)
        candidate = self._validated_update_graph(candidate)

        self._check_update_graph(candidate)
        self._apply_update_graph(candidate)

    def edge_count(
        self,
        source: str,
        target: str,
        sign: Optional[InfluenceSign] = None,
    ) -> int:
        """
        Return the occurrence count of an aggregated signed influence.

        If both positive and negative edges exist between the same source and
        target, `sign` must be supplied to select one signed influence.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        >>> graph.edge_count("A", "B")
        3

        >>> graph.add_edge("A", "B", sign=-1, count=1)
        >>> graph.edge_count("A", "B", sign=-1)
        1

        Parameters
        ----------
        source: str
            Source node.
        target: str
            Target node.
        sign: {-1, 1}, optional
            Signed influence to select when several signed edges exist between
            `source` and `target`.

        Returns
        -------
        int
            Occurrence count of the selected signed influence.

        Raises
        ------
        KeyError
            If no matching edge exists.
        ValueError
            If both positive and negative edges exist and `sign` is not provided.
        """

        return int(self._aggregated_edge_data(source, target, sign)["count"])

    def edge_frequency(
        self,
        source: str,
        target: str,
        sign: Optional[InfluenceSign] = None,
    ) -> float:
        """
        Return the occurrence frequency of an aggregated signed influence.

        The frequency is computed as `edge_count(source, target, sign) / total`.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        >>> graph.edge_frequency("A", "B")
        0.75

        Parameters
        ----------
        source: str
            Source node.
        target: str
            Target node.
        sign: {-1, 1}, optional
            Signed influence to select when several signed edges exist between
            `source` and `target`.

        Returns
        -------
        float
            Occurrence frequency of the selected signed influence.

        Raises
        ------
        KeyError
            If no matching edge exists.
        ValueError
            If both positive and negative edges exist and `sign` is not provided.
        """

        return self.edge_count(source, target, sign=sign) / self.total

    # Aggregated autoregulations include frequencies by design.
    def autoregulations(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        sign: Optional[InfluenceSign] = None,
    ) -> List[Tuple[str, int, float]]:
        """
        Return autoregulations with their sign and occurrence frequency.

        Autoregulations correspond to self-loops. Edge signs are normalized to
        -1 or 1, and frequencies are computed from edge counts and `total`.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "A", sign=1, count=3)
        >>> graph.add_edge("B", "B", sign=-1, count=1)
        >>> graph.autoregulations()
        [('A', 1, 0.75), ('B', -1, 0.25)]

        >>> graph.autoregulations(sign=1)
        [('A', 1, 0.75)]

        Parameters
        ----------
        sign: {-1, 1}, optional
            If provided, return only autoregulations with the requested sign.

        Returns
        -------
        list of tuple
            Tuples `(node, sign, frequency)`.

        Raises
        ------
        ValueError
            If `sign` is provided and is not a valid edge sign.
        """

        loops = []
        normalized_sign = None if sign is None else self._normalize_sign(sign)

        for node in nx.nodes_with_selfloops(self):
            edge_data = self.get_edge_data(node, node, default={})

            for data in edge_data.values():
                edge_sign = self._normalize_sign(data["sign"])

                if normalized_sign is None or edge_sign == normalized_sign:
                    loops.append((node, edge_sign, data["count"] / self.total))

        return loops

    def frequency_bin(
        self,
        source: str,
        target: str,
        bins: Iterable[float],
        sign: Optional[InfluenceSign] = None,
    ) -> Tuple[float, float]:
        """
        Return the frequency interval containing an edge frequency.

        Intervals are interpreted as closed intervals over consecutive bin
        values: `(bins[i], bins[i + 1])`. If a frequency is exactly on a shared
        boundary, it belongs to both adjacent intervals and the first matching
        interval is returned.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        >>> graph.frequency_bin("A", "B", bins=(0.0, 0.5, 0.75, 1.0))
        (0.5, 0.75)

        Parameters
        ----------
        source: str
            Source node.
        target: str
            Target node.
        bins: Iterable of float
            Ordered bin boundaries used to classify frequencies.
        sign: {-1, 1}, optional
            Signed influence to select when several signed edges exist between
            `source` and `target`.

        Returns
        -------
        tuple of float
            Lower and upper bounds of the matching frequency interval.

        Raises
        ------
        KeyError
            If no matching edge exists.
        ValueError
            If the edge is ambiguous without `sign`, or if its frequency is not
            covered by `bins`.
        """

        bins = tuple(bins)
        frequency = self.edge_frequency(source, target, sign=sign)

        for lower, upper in zip(bins[:-1], bins[1:]):
            if lower <= frequency <= upper:
                return (lower, upper)

        raise ValueError(
            f"edge frequency {frequency!r} is not covered by bins {bins!r}"
        )

    def structural_families(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
        protect_feedback_nodes: Optional[bool] = None,
    ) -> Dict[Tuple[Any, ...], Set[str]]:
        """
        Group structurally equivalent nodes using signed frequency-aware edges.

        Two nodes are grouped only if they have the same signed predecessor and,
        optionally, successor structure with the same frequency bins.

        Frequency bins are part of the structural signature. Therefore two nodes
        with the same signed regulators are separated if the corresponding edge
        counts fall in different frequency intervals.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edges_from(
        ...     [
        ...         ("TF", "g1", {"sign": 1, "count": 3}),
        ...         ("TF", "g2", {"sign": 1, "count": 3}),
        ...         ("TF", "g3", {"sign": 1, "count": 1}),
        ...     ]
        ... )
        >>> families = graph.structural_families(
        ...     include_successors=False,
        ...     bins=(0.0, 0.5, 1.0),
        ... )
        >>> set(map(frozenset, families.values()))
        {frozenset({'g1', 'g2'})}

        Parameters
        ----------
        include_successors: bool (default: True)
            Whether signed successor signatures should also be included in the
            structural equivalence criterion.
        exclude_feedback_nodes: bool (default: True)
            Whether nodes participating in feedback structures should be excluded.
        min_size: int (default: 2)
            Minimum family size required to keep a structural family.
        bins: Iterable of float
            Ordered bin boundaries used to classify edge frequencies.
        protect_feedback_nodes: bool, optional
            Backward-compatible alias for `exclude_feedback_nodes`. If provided,
            it takes precedence over `exclude_feedback_nodes`.

        Returns
        -------
        dict
            Mapping from frequency-aware structural signatures to node families.
        """

        if protect_feedback_nodes is not None:
            exclude_feedback_nodes = protect_feedback_nodes

        bins = tuple(bins)
        feedback_nodes = self.feedback_nodes() if exclude_feedback_nodes else set()
        families = {}

        for node in self.nodes():
            if node in feedback_nodes:
                continue

            predecessors = tuple(
                sorted(
                    (
                        predecessor,
                        self._normalize_sign(data["sign"]),
                        self._frequency_bin_for_count(data["count"], bins),
                    )
                    for predecessor, _, data in self.in_edges(node, data=True)
                )
            )

            if include_successors:
                successors = tuple(
                    sorted(
                        (
                            successor,
                            self._normalize_sign(data["sign"]),
                            self._frequency_bin_for_count(data["count"], bins),
                        )
                        for _, successor, data in self.out_edges(node, data=True)
                    )
                )

            else:
                successors = tuple()

            signature = (predecessors, successors)
            families.setdefault(signature, set()).add(node)

        return {
            signature: family
            for signature, family in families.items()
            if len(family) >= min_size
        }

    def family_collapsed_graph(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
        sep: str = "|",
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    ) -> InfluenceGraph:
        """
        Return a frequency-aware graph where structural families are collapsed.

        Nodes belonging to the same frequency-aware structural family are
        replaced by a composite node. The composite node name is obtained by
        joining member names with `sep`, and original members are stored in the
        node attribute `members`.

        The returned graph is a plain `InfluenceGraph`, not an
        `AggregatedInfluenceGraph`: after node fusion, edge counts no longer have
        the exact meaning “number of source graphs containing this signed
        interaction”. Collapsed edges therefore store `frequency`,
        `min_frequency` and `max_frequency` metadata.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edges_from(
        ...     [
        ...         ("TF", "g1", {"sign": 1, "count": 3}),
        ...         ("TF", "g2", {"sign": 1, "count": 3}),
        ...         ("g1", "out", {"sign": -1, "count": 2}),
        ...         ("g2", "out", {"sign": -1, "count": 2}),
        ...     ]
        ... )
        >>> collapsed = graph.family_collapsed_graph(
        ...     exclude_feedback_nodes=False,
        ... )
        >>> sorted(collapsed.nodes())
        ['TF', 'g1|g2', 'out']

        >>> collapsed.nodes["g1|g2"]["members"]
        {'g1', 'g2'}

        >>> collapsed["TF"]["g1|g2"][0]["frequency"]
        0.75

        Parameters
        ----------
        include_successors: bool (default: True)
            Whether signed successor signatures should also be included in
            structural family detection.
        exclude_feedback_nodes: bool (default: True)
            Whether nodes participating in feedback structures should be excluded
            from collapsing.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.
        bins: Iterable of float
            Ordered bin boundaries used to classify edge frequencies.

        Returns
        -------
        InfluenceGraph
            Plain influence graph with collapsed nodes and frequency metadata
            on edges.

        Raises
        ------
        ValueError
            If an edge frequency is not covered by `bins`.
        """

        families = self.structural_families(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            bins=bins,
        )

        node_mapping = {}

        for family in families.values():
            family_name = sep.join(sorted(family))

            for node in family:
                node_mapping[node] = family_name

        collapsed = InfluenceGraph()

        for node, data in self.nodes(data=True):
            collapsed_node = node_mapping.get(node, node)

            if collapsed_node not in collapsed:
                collapsed.add_node(collapsed_node)

            if node in node_mapping:
                collapsed.nodes[collapsed_node].setdefault("members", set())
                collapsed.nodes[collapsed_node]["members"].add(node)

            else:
                collapsed.nodes[collapsed_node].update(data)
                collapsed.nodes[collapsed_node].setdefault("members", {node})

        edge_frequencies = {}

        for source, target, data in self.edges(data=True):
            collapsed_source = node_mapping.get(source, source)
            collapsed_target = node_mapping.get(target, target)

            if collapsed_source == collapsed_target:
                continue

            sign = self._normalize_sign(data["sign"])
            key = (collapsed_source, collapsed_target, sign)
            edge_frequencies.setdefault(key, []).append(data["count"] / self.total)

        for (source, target, sign), frequencies in edge_frequencies.items():
            frequency = sum(frequencies) / len(frequencies)

            collapsed.add_edge(
                source,
                target,
                sign=sign,
                frequency=frequency,
                min_frequency=min(frequencies),
                max_frequency=max(frequencies),
            )

        return collapsed

    def family_compressed_graph(
        self,
        *_args: Any,
        **_kwargs: Any,
    ) -> NoReturn:
        """
        Disable the deprecated `InfluenceGraph.family_compressed_graph()` alias.

        `AggregatedInfluenceGraph` uses collapse terminology directly. Use
        `family_collapsed_graph()` instead.

        Raises
        ------
        NotImplementedError
            Always raised.
        """

        raise NotImplementedError(
            "AggregatedInfluenceGraph does not support "
            "family_compressed_graph(); use family_collapsed_graph() instead."
        )

    def feedback_induced_graph(
        self,
        include_singleton_selfloops: bool = True,
    ) -> "AggregatedInfluenceGraph":
        """
        Return the aggregated graph induced by feedback nodes.

        Feedback nodes are computed as in `InfluenceGraph.feedback_induced_graph`.
        The returned graph preserves `count` attributes and the aggregation total.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=3)
        >>> graph.add_edges_from(
        ...     [
        ...         ("A", "B", {"sign": 1, "count": 2}),
        ...         ("B", "A", {"sign": 1, "count": 2}),
        ...         ("B", "C", {"sign": -1, "count": 1}),
        ...     ]
        ... )
        >>> feedback = graph.feedback_induced_graph()
        >>> sorted(feedback.nodes())
        ['A', 'B']

        >>> feedback.total
        3

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should
            be included.

        Returns
        -------
        AggregatedInfluenceGraph
            Aggregated influence graph induced by feedback nodes.
        """

        feedback_nodes = self.feedback_nodes(
            include_singleton_selfloops=include_singleton_selfloops,
        )

        return type(self)(
            self._plain_graph(self.subgraph(feedback_nodes)),
            total=self.total,
        )

    def collapsed_graph(
        self,
        include_singleton_selfloops: bool = True,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = False,
        min_size: int = 2,
        sep: str = "|",
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    ) -> InfluenceGraph:
        """
        Return a collapsed influence graph for visualization or summaries.

        Collapse is applied sequentially through:
            1. feedback-induced graph extraction,
            2. frequency-aware structural family collapse.

        The returned graph is a plain `InfluenceGraph` with frequency metadata
        on edges, not an exact `AggregatedInfluenceGraph`.

        Parameters
        ----------
        include_singleton_selfloops: bool (default: True)
            Whether singleton SCCs corresponding to self-regulated nodes should
            be included in the feedback-induced graph.
        include_successors: bool (default: True)
            Whether signed successor signatures should be included in structural
            family detection.
        exclude_feedback_nodes: bool (default: False)
            Whether nodes participating in feedback structures should be excluded
            from structural family collapse.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.
        bins: Iterable of float
            Ordered bin boundaries used to classify edge frequencies.

        Returns
        -------
        InfluenceGraph
            Collapsed influence graph with `frequency`, `min_frequency` and
            `max_frequency` edge attributes.

        See Also
        --------
        feedback_induced_graph
        family_collapsed_graph
        """

        graph = self.feedback_induced_graph(
            include_singleton_selfloops=include_singleton_selfloops,
        )

        return graph.family_collapsed_graph(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            sep=sep,
            bins=bins,
        )

    def compressed_graph(
        self,
        *_args: Any,
        **_kwargs: Any,
    ) -> NoReturn:
        """
        Disable the deprecated `InfluenceGraph.compressed_graph()` alias.

        `AggregatedInfluenceGraph` uses collapse terminology directly. Use
        `collapsed_graph()` instead.

        Raises
        ------
        NotImplementedError
            Always raised.
        """

        raise NotImplementedError(
            "AggregatedInfluenceGraph does not support compressed_graph(); "
            "use collapsed_graph() instead."
        )

    # Aggregated rendering exposes collapse/frequency controls.
    def to_graphviz(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        collapse: Optional[CollapseMode] = None,
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
        protect_feedback_nodes: bool = True,
        include_singleton_selfloops: bool = True,
        min_frequency: float = 0.0,
        show_edge_labels: bool = True,
        edge_style: Union[
            Callable[[float], Mapping[str, Any]],
            bool,
            None,
        ] = ratio_edge_style,
        program: str = "dot",
        **kwargs: Any,
    ):
        """
        Convert the aggregated influence graph to a native graphviz Digraph.

        Rendering is always performed from an `InfluenceGraph`. Without
        collapse, the graph preserves exact edge counts and labels display
        counts. With `collapse="family"`, `collapse="feedback"` or
        `collapse="both"`, the rendered graph uses frequency labels.

        Parameters
        ----------
        collapse: {None, "family", "feedback", "both"} (default: None)
            Optional collapse applied before rendering. `"family"` uses
            `family_collapsed_graph()`. `"feedback"` keeps only the
            feedback-induced graph. `"both"` uses `collapsed_graph()`.
        bins: Iterable of float
            Ordered bin boundaries used by collapsed structural families.
        protect_feedback_nodes: bool (default: True)
            Whether feedback nodes are protected from family collapse.
        include_singleton_selfloops: bool (default: True)
            Whether singleton self-loop SCCs are retained when using
            `collapse="feedback"`.
        min_frequency: float (default: 0.0)
            Minimum edge frequency required for display.
        show_edge_labels: bool (default: True)
            Whether to display counts or frequencies as edge labels.
        edge_style: Callable[[float], Mapping[str, Any]] or bool or None
            Function used to style edges according to frequency. If True, use
            `ratio_edge_style`. If False or None, no frequency style is applied.
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting graph.
        **kwargs: Any
            Graph attributes assigned to the resulting graphviz object.

        Returns
        -------
        graphviz.Digraph
            Styled graphviz object.

        Raises
        ------
        ImportError
            If the `graphviz` Python package is not installed.
        ValueError
            If `collapse` is invalid or `min_frequency` is outside [0, 1].
        """

        graph = self._visualization_graph(
            collapse=collapse,
            bins=bins,
            protect_feedback_nodes=protect_feedback_nodes,
            include_singleton_selfloops=include_singleton_selfloops,
            min_frequency=min_frequency,
            show_edge_labels=show_edge_labels,
            edge_style=edge_style,
        )

        return _networkx_to_graphviz(graph, program=program, **kwargs)

    # Aggregated rendering exposes collapse/frequency controls.
    def to_pydot(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        collapse: Optional[CollapseMode] = None,
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
        protect_feedback_nodes: bool = True,
        include_singleton_selfloops: bool = True,
        min_frequency: float = 0.0,
        show_edge_labels: bool = True,
        edge_style: Union[
            Callable[[float], Mapping[str, Any]],
            bool,
            None,
        ] = ratio_edge_style,
        program: str = "dot",
        **kwargs: Any,
    ) -> "Dot":
        """
        Convert the aggregated influence graph to a pydot graph.

        Rendering is always performed from an `InfluenceGraph`. Without
        collapse, the graph preserves exact edge counts and labels display
        counts. With `collapse="family"`, `collapse="feedback"` or
        `collapse="both"`, the rendered graph uses frequency labels.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        0
        >>> dot = graph.to_pydot(rankdir="LR")
        >>> dot.get_rankdir()
        'LR'

        Parameters
        ----------
        collapse: {None, "family", "feedback", "both"} (default: None)
            Optional collapse applied before rendering. `"family"` uses
            `family_collapsed_graph()`. `"feedback"` keeps only the
            feedback-induced graph. `"both"` uses `collapsed_graph()`.
        bins: Iterable of float
            Ordered bin boundaries used by collapsed structural families.
        protect_feedback_nodes: bool (default: True)
            Whether feedback nodes are protected from family collapse.
        include_singleton_selfloops: bool (default: True)
            Whether singleton self-loop SCCs are retained when using
            `collapse="feedback"`.
        min_frequency: float (default: 0.0)
            Minimum edge frequency required for display.
        show_edge_labels: bool (default: True)
            Whether to display counts or frequencies as edge labels.
        edge_style: Callable[[float], Mapping[str, Any]] or bool or None
            Function used to style edges according to frequency. If True, use
            `ratio_edge_style`. If False or None, no frequency style is applied.
        program: str (default: "dot")
            Graphviz layout program assigned to the resulting pydot graph.
        **kwargs: Any
            Keyword arguments passed to the resulting pydot graph using
            `dot.set(key, value)`.

        Returns
        -------
        Dot
            Styled pydot graph.

        Raises
        ------
        ValueError
            If `collapse` is invalid or `min_frequency` is outside [0, 1].
        """

        graph = self._visualization_graph(
            collapse=collapse,
            bins=bins,
            protect_feedback_nodes=protect_feedback_nodes,
            include_singleton_selfloops=include_singleton_selfloops,
            min_frequency=min_frequency,
            show_edge_labels=show_edge_labels,
            edge_style=edge_style,
        )

        dot = nx.drawing.nx_pydot.to_pydot(graph)

        dot.set_prog(program)

        for key, value in kwargs.items():
            dot.set(key, value)

        return dot

    # Aggregated rendering exposes collapse/frequency controls.
    def show(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        collapse: Optional[CollapseMode] = None,
        bins: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
        protect_feedback_nodes: bool = True,
        include_singleton_selfloops: bool = True,
        min_frequency: float = 0.0,
        show_edge_labels: bool = True,
        edge_style: Union[
            Callable[[float], Mapping[str, Any]],
            bool,
            None,
        ] = ratio_edge_style,
        program: str = "dot",
        width: Optional[SvgLength] = None,
        height: Optional[SvgLength] = None,
        **kwargs: Any,
    ) -> None:
        """
        Display the aggregated influence graph in a Jupyter/IPython environment.

        The graph is rendered through Graphviz using the `to_pydot()` method and
        displayed as an SVG image. Optional collapse modes allow rendering the
        exact graph, the family-collapsed graph, the feedback-induced graph, or
        the full collapsed graph.

        Examples
        --------
        Display the exact aggregated graph:

        >>> graph = AggregatedInfluenceGraph(total=4)
        >>> graph.add_edge("A", "B", sign=1, count=3)
        0
        >>> graph.show(width="700px")  # doctest: +SKIP

        Display a family-collapsed graph:

        >>> graph.show(collapse="family", width=900)  # doctest: +SKIP

        Display only feedback influences:

        >>> graph.show(collapse="feedback")  # doctest: +SKIP

        Parameters
        ----------
        collapse: {None, "family", "feedback", "both"} (default: None)
            Optional graph transformation applied before rendering.
            If None, render the exact aggregated graph and display edge counts.
            If `"family"`, render `family_collapsed_graph()` and display edge
            frequencies.
            If `"feedback"`, render the feedback-induced graph and display edge
            frequencies.
            If `"both"`, render `collapsed_graph()` and display edge
            frequencies.
        bins: Iterable of float
            Ordered bin boundaries used by collapsed structural families.
        protect_feedback_nodes: bool (default: True)
            Whether feedback nodes are protected from family collapse when
            `collapse="family"`.
        include_singleton_selfloops: bool (default: True)
            Whether singleton self-loop SCCs are retained when using
            `collapse="feedback"` or `collapse="both"`.
        min_frequency: float (default: 0.0)
            Minimum edge frequency required for display.
        show_edge_labels: bool (default: True)
            Whether to display counts or frequencies as edge labels.
        edge_style: Callable[[float], Mapping[str, Any]] or bool or None
            Function used to style edges according to frequency. If True, use
            `ratio_edge_style`. If False or None, no frequency style is applied.
        program: str (default: "dot")
            Graphviz layout program used for rendering.
        width: str or int or float, optional
            Display width assigned to the rendered SVG root. Integers and floats
            are passed as raw SVG length values; strings can include CSS units,
            for example `"900px"` or `"80%"`.
        height: str or int or float, optional
            Display height assigned to the rendered SVG root. Integers and
            floats are passed as raw SVG length values; strings can include CSS
            units.
        **kwargs: Any
            Keyword arguments passed to the underlying pydot graph through
            `dot.set(key, value)`.

        Returns
        -------
        None
            The SVG is displayed in the current IPython/Jupyter output cell.

        Raises
        ------
        RuntimeError
            If IPython is not available.
        ValueError
            If `collapse` is invalid or `min_frequency` is outside [0, 1].

        Notes
        -----
        `width` and `height` only affect notebook display by rewriting the root
        SVG attributes after Graphviz rendering. They do not change the graph
        layout computed by Graphviz.
        """

        try:
            from IPython.display import SVG, display

        except ImportError:
            raise RuntimeError("show() requires an IPython/Jupyter environment.")

        dot = cast(
            Any,
            self.to_pydot(
                collapse=collapse,
                bins=bins,
                protect_feedback_nodes=protect_feedback_nodes,
                include_singleton_selfloops=include_singleton_selfloops,
                min_frequency=min_frequency,
                show_edge_labels=show_edge_labels,
                edge_style=edge_style,
                program=program,
                **kwargs,
            ),
        )
        svg = scale_svg(dot.create_svg().decode(), width=width, height=height)

        display(SVG(svg))

    def validate_counts(self) -> None:
        """
        Validate aggregated edge counts.

        Every edge must define a `count` attribute. Counts must be integers
        between 0 and `total`.

        Examples
        --------
        >>> graph = AggregatedInfluenceGraph(total=3)
        >>> graph.add_edge("A", "B", sign=1, count=2)
        >>> graph.validate_counts()

        Invalid counts are rejected:

        >>> invalid = nx.MultiDiGraph()
        >>> invalid.add_edge("A", "B", sign=1, count=4)
        >>> AggregatedInfluenceGraph(invalid, total=3)
        Traceback (most recent call last):
            ...
        ValueError: invalid edge count for edge 'A' -> 'B':
            expected value between 0 and 3 but received 4

        Raises
        ------
        KeyError
            If an edge does not define `count`.
        TypeError
            If an edge count is not an integer.
        ValueError
            If an edge count is negative or greater than `total`.
        """

        for source, target, data in self.edges(data=True):
            if "count" not in data:
                raise KeyError(
                    f"missing required edge attribute 'count' "
                    f"for edge {source!r} -> {target!r}"
                )

            count = data["count"]

            self._validate_count(source, target, count)

    @classmethod
    def from_influence_graphs(
        cls,
        *graphs: InfluenceGraph,
    ) -> "AggregatedInfluenceGraph":
        """
        Build an aggregated influence graph from influence graphs.

        Each input graph contributes at most one occurrence to each signed edge.
        The resulting `total` is the number of input graphs, and every aggregated
        edge stores the number of graphs in which the same signed influence was
        observed.

        Positive and negative influences between the same source and target are
        aggregated as distinct signed edges, matching `InfluenceGraph`
        multiedge semantics.

        Examples
        --------
        >>> ig1 = InfluenceGraph()
        >>> ig1.add_edge("A", "B", sign=1)
        >>> ig1.add_edge("A", "B", sign=-1)

        >>> ig2 = InfluenceGraph()
        >>> ig2.add_edge("A", "B", sign=1)
        >>> ig2.add_edge("B", "C", sign=-1)

        >>> graph = AggregatedInfluenceGraph.from_influence_graphs(ig1, ig2)
        >>> graph.total
        2

        >>> graph.edge_count("A", "B", sign=1)
        2

        >>> graph.edge_count("A", "B", sign=-1)
        1

        Parameters
        ----------
        *graphs: InfluenceGraph
            Influence graphs to aggregate. At least one graph is required.

        Returns
        -------
        AggregatedInfluenceGraph
            Aggregated graph whose edge counts indicate in how many input graphs
            each signed influence was observed.

        Raises
        ------
        ValueError
            If no graph is provided, or if one input graph contains an invalid
            edge sign.
        """

        total = len(graphs)

        if total == 0:
            raise ValueError("expected at least one influence graph")

        aggregated = cls(total=total)

        for graph in graphs:
            aggregated.add_nodes_from(
                (node, data.copy()) for node, data in graph.nodes(data=True)
            )

            for source, target, data in graph.edges(data=True):
                sign = cast(CircuitSign, cls._normalize_sign(data["sign"]))

                try:
                    edge_data = aggregated._aggregated_edge_data(
                        source,
                        target,
                        sign=sign,
                    )
                except KeyError:
                    aggregated.add_edge(source, target, sign=sign, count=1)
                else:
                    edge_data["count"] += 1

        aggregated.validate_counts()

        return aggregated

    @classmethod
    def from_boolean_networks(
        cls,
        *networks: Union["BooleanNetwork", "BooleanNetworkEnsemble"],
    ) -> "AggregatedInfluenceGraph":
        """
        Build an aggregated influence graph from Boolean networks.

        Each Boolean network is converted with `to_influence_graph()` before
        aggregation. A single `BooleanNetworkEnsemble` argument is accepted and
        expanded automatically.

        Examples
        --------
        >>> from bonesistools.boolpy.boolean_network import BooleanNetwork

        >>> bn1 = BooleanNetwork({"A": "B", "B": 0, "C": 1})
        >>> bn2 = BooleanNetwork({"A": "B", "B": "C", "C": 0})
        >>> bn3 = BooleanNetwork({"A": "!B", "B": "C", "C": 1})

        >>> graph = AggregatedInfluenceGraph.from_boolean_networks(
        ...     bn1,
        ...     bn2,
        ...     bn3,
        ... )
        >>> graph.edge_count("B", "A", sign=1)
        2

        >>> graph.edge_count("B", "A", sign=-1)
        1

        A Boolean network ensemble can also be passed directly:

        >>> from bonesistools.boolpy.boolean_network import BooleanNetworkEnsemble
        >>> ensemble = BooleanNetworkEnsemble(bns=[bn1, bn2, bn3])
        >>> graph = AggregatedInfluenceGraph.from_boolean_networks(ensemble)
        >>> graph.total
        3

        Parameters
        ----------
        *networks: BooleanNetwork
            Boolean networks to convert and aggregate, or a single
            `BooleanNetworkEnsemble`. At least one network is required.

        Returns
        -------
        AggregatedInfluenceGraph
            Aggregated influence graph built from the influence graphs inferred
            from the input Boolean networks.

        Raises
        ------
        AttributeError
            If one input object does not provide `to_influence_graph()`.
        ValueError
            If no network is provided, or if conversion produces an influence
            graph with invalid edge signs.
        """

        from ..boolean_network import BooleanNetworkEnsemble

        if len(networks) == 1 and isinstance(networks[0], BooleanNetworkEnsemble):
            boolean_networks = cast(Tuple["BooleanNetwork", ...], tuple(networks[0]))
        else:
            boolean_networks = cast(Tuple["BooleanNetwork", ...], networks)

        graphs = tuple(network.to_influence_graph() for network in boolean_networks)

        return cls.from_influence_graphs(*graphs)

    def _visualization_graph(
        self,
        collapse: Optional[CollapseMode],
        bins: Iterable[float],
        protect_feedback_nodes: bool,
        include_singleton_selfloops: bool,
        min_frequency: float,
        show_edge_labels: bool,
        edge_style: Union[
            Callable[[float], Mapping[str, Any]],
            bool,
            None,
        ],
    ) -> InfluenceGraph:
        """
        Return a styled InfluenceGraph for pydot or graphviz rendering.
        """

        min_frequency = _as_probability(min_frequency, "min_frequency")

        if collapse is None:
            graph = InfluenceGraph(self._plain_graph(self))

        elif collapse == "family":
            graph = self.family_collapsed_graph(
                exclude_feedback_nodes=protect_feedback_nodes,
                bins=bins,
            )

        elif collapse == "feedback":
            graph = InfluenceGraph(
                self._plain_graph(
                    self.feedback_induced_graph(
                        include_singleton_selfloops=include_singleton_selfloops,
                    )
                )
            )

        elif collapse == "both":
            graph = self.collapsed_graph(
                include_singleton_selfloops=include_singleton_selfloops,
                exclude_feedback_nodes=protect_feedback_nodes,
                bins=bins,
            )

        else:
            raise ValueError(
                "unsupported collapse: expected None, 'family', 'feedback' "
                f"or 'both', but received {collapse!r}"
            )

        edge_style_callable = None

        if edge_style is True:
            edge_style_callable = ratio_edge_style

        elif callable(edge_style):
            edge_style_callable = edge_style

        edges_to_remove = []

        for source, target, key, data in graph.edges(keys=True, data=True):
            frequency = data.get("frequency")

            if frequency is None:
                frequency = data["count"] / self.total
                data["frequency"] = frequency

            if frequency < min_frequency:
                edges_to_remove.append((source, target, key))
                continue

            sign = self._normalize_sign(data["sign"])

            data.update(
                color="green4" if sign == 1 else "red2",
                arrowhead="normal" if sign == 1 else "tee",
                penwidth="2",
            )

            if show_edge_labels:
                if collapse is None:
                    data["label"] = str(data["count"])
                else:
                    data["label"] = f"{frequency:.2g}"

            if edge_style_callable is not None:
                data.update(edge_style_callable(frequency))

        graph.remove_edges_from(edges_to_remove)

        return graph

    def _empty_like(self) -> "AggregatedInfluenceGraph":
        """
        Return an empty aggregated graph preserving the current total.
        """

        return type(self)(total=self.total)

    def _validated_update_graph(self, graph: Any) -> "AggregatedInfluenceGraph":
        """
        Return a validated aggregated graph containing candidate update data only.
        """

        return type(self)(graph, total=self.total)

    def _replace_with_graph(self, graph: Any) -> None:
        """
        Replace graph contents while preserving aggregated invariants.
        """

        plain_graph = self._plain_graph(graph)

        self.clear()
        self.graph.update(plain_graph.graph)
        nx.MultiDiGraph.add_nodes_from(
            self,
            ((node, data.copy()) for node, data in plain_graph.nodes(data=True)),
        )

        for source, target, data in plain_graph.edges(data=True):
            data = data.copy()

            if "sign" not in data:
                raise ValueError(
                    f"missing edge attribute 'sign' for edge "
                    f"{source!r} -> {target!r}"
                )

            if "count" not in data:
                raise KeyError(
                    f"missing required edge attribute 'count' "
                    f"for edge {source!r} -> {target!r}"
                )

            sign = data.pop("sign")
            self.add_edge(source, target, sign=sign, **data)

    def _aggregated_edge_data(
        self,
        source: str,
        target: str,
        sign: Optional[InfluenceSign] = None,
    ) -> MutableMapping[str, Any]:
        """
        Return one aggregated edge data mapping.
        """

        edge_data = self.get_edge_data(source, target)

        if edge_data is None:
            raise KeyError(f"no edge found between {source!r} and {target!r}")

        if sign is not None:
            normalized_sign = self._normalize_sign(sign)

            for data in edge_data.values():
                if self._normalize_sign(data["sign"]) == normalized_sign:
                    return data

            raise KeyError(
                f"no edge found between {source!r} and {target!r} "
                f"with sign {normalized_sign!r}"
            )

        if len(edge_data) != 1:
            raise ValueError(
                f"ambiguous aggregated edge between {source!r} and {target!r}: "
                "pass sign=... to select one signed influence"
            )

        return next(iter(edge_data.values()))

    def _frequency_bin_for_count(
        self,
        count: int,
        bins: Iterable[float],
    ) -> Tuple[float, float]:
        """
        Return the frequency interval containing an edge count.
        """

        bins = tuple(bins)
        frequency = count / self.total

        for lower, upper in zip(bins[:-1], bins[1:]):
            if lower <= frequency <= upper:
                return (lower, upper)

        raise ValueError(
            f"edge frequency {frequency!r} is not covered by bins {tuple(bins)!r}"
        )

    def _validate_count(
        self,
        source: Any,
        target: Any,
        count: Any,
        total: Optional[int] = None,
    ) -> int:
        """
        Validate one aggregated edge count.
        """

        if total is None:
            total = self.total

        count = _as_non_negative_integer(count, "count")

        if count > total:
            raise ValueError(
                f"invalid edge count for edge {source!r} -> {target!r}: "
                f"expected value between 0 and {total} "
                f"but received {count!r}"
            )

        return count

    @staticmethod
    def _validate_total(total: Any) -> int:
        """
        Validate an aggregation total.
        """

        return _as_positive_integer(total, "total")
