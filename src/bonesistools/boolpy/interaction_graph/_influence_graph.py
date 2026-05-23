#!/usr/bin/env python

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING, Any, Iterable, Literal, Mapping, Optional, Callable

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import networkx as nx

if TYPE_CHECKING:
    from pydot import Dot

CircuitSign = Literal[-1, 1]
Direction = Literal["upstream", "downstream", "both"]

StructuralSignature = tuple[
    frozenset[tuple[str, int]],
    frozenset[tuple[str, int]],
]


class InfluenceGraph(nx.MultiDiGraph):
    """
    Signed influence graph.

    InfluenceGraph is a domain-specific NetworkX MultiDiGraph for signed
    regulatory influence graphs. It enforces signed edges and provides analyses
    specific to logical and regulatory graphs, such as feedback circuits, signed
    autoregulations, strongly connected components and graph compression.

    Edge signs must be stored in the edge attribute `sign`, using either:
        - 1, "+", or "positive" for positive influences,
        - -1, "-", or "negative" for negative influences.

    At most one positive edge and one negative edge may exist between a pair of
    nodes.
    """

    def __init__(
        self,
        graph: Optional[nx.DiGraph] = None,
        **attr: Any,
    ) -> None:

        super().__init__(**attr)

        if graph is not None:
            self._replace_with_graph(self._plain_graph(graph))
        else:
            self._validate_graph()

    __slots__ = ()

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

    def copy(self, as_view: bool = False) -> Self:
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

    def to_directed(self, as_view: bool = False) -> Self:
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

    def to_undirected(self, *args: Any, **kwargs: Any) -> None:
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

    def add_edge(
        self,
        source: Any,
        target: Any,
        sign: CircuitSign,
        **attr: Any,
    ) -> None:
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

        sign = self._normalize_sign(sign)

        if self.has_edge(source, target):
            signs = self._edge_signs(source, target)

            if sign in signs:
                raise ValueError(
                    f"duplicated edge sign for edge "
                    f"{source!r} -> {target!r} with sign {sign!r}"
                )

        attr["sign"] = sign

        super().add_edge(
            source,
            target,
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

        candidate = self._plain_graph(self)
        candidate.add_edges_from(ebunch_to_add, **attr)
        candidate = type(self)(candidate)

        self._replace_with_graph(candidate)

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

        candidate = self._plain_graph(self)
        nx.MultiDiGraph.update(candidate, edges=edges, nodes=nodes)
        candidate = type(self)(candidate)

        self._replace_with_graph(candidate)

    def strongly_connected_components(
        self,
        include_singleton_selfloops: bool = True,
    ) -> list[set[str]]:
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
    ) -> set[str]:
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

    def regulators(self) -> list[Any]:
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

    def targets(self) -> list[Any]:
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
    ) -> dict[StructuralSignature, set[str]]:
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
        information that should not be compressed.

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

    def family_compressed_graph(
        self,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = True,
        min_size: int = 2,
        sep: str = "|",
    ) -> Self:
        """
        Return a graph where structurally equivalent nodes are collapsed.

        Nodes belonging to the same structural family are replaced by a single
        composite node. The composite node name is obtained by joining member names
        with `sep`, and the original members are stored in the node attribute
        `members`.

        Edges incident to collapsed nodes are rewired to the corresponding composite
        node. Since structural families are defined from signed predecessor and
        successor signatures, this compression preserves signed influence structure.

        Examples
        --------
        Consider the following influence graph:

            TF1 → g1
            TF1 → g2
            TF1 → g3
            TF1 → g4

            TF2 ─| g3
            TF2 ─| g4

        >>> compressed = ig.family_compressed_graph()
        >>> sorted(compressed.nodes())
        ['TF1', 'TF2', 'g1|g2', 'g3|g4']

        The composite nodes store their original members:

        >>> compressed.nodes["g1|g2"]["members"]
        {'g1', 'g2'}

        Parameters
        ----------
        include_successors: bool (default: True)
            Whether signed successor signatures should also be included in the
            structural equivalence criterion.
        exclude_feedback_nodes: bool (default: True)
            Whether nodes participating in feedback structures should be excluded
            from compression.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.

        Returns
        -------
        InfluenceGraph
            Family-compressed influence graph.
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

        compressed = type(self)()

        for node, data in self.nodes(data=True):
            compressed_node = node_mapping.get(node, node)

            if compressed_node not in compressed:
                compressed.add_node(compressed_node)

            if node in node_mapping:
                compressed.nodes[compressed_node].setdefault("members", set())
                compressed.nodes[compressed_node]["members"].add(node)

            else:
                compressed.nodes[compressed_node].update(data)
                compressed.nodes[compressed_node].setdefault("members", {node})

        for source, target, data in self.edges(data=True):
            compressed_source = node_mapping.get(source, source)
            compressed_target = node_mapping.get(target, target)

            if compressed_source == compressed_target:
                continue

            sign = self._normalize_sign(data["sign"])

            if compressed.has_edge(
                compressed_source, compressed_target
            ) and sign in compressed._edge_signs(
                compressed_source,
                compressed_target,
            ):
                continue

            compressed.add_edge(
                compressed_source,
                compressed_target,
                **data,
            )

        return compressed

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
        sign: Optional[CircuitSign] = None,
    ) -> list[tuple[str, int]]:
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

    def marker_paths(
        self,
        markers: Iterable[str],
        direction: Direction = "both",
        sccs: Optional[Iterable[Iterable[str]]] = None,
        sign: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """
        Return shortest paths between markers and feedback SCCs.

        Paths are computed between each marker and each feedback-relevant strongly
        connected component (SCC). For a given SCC/marker pair, all shortest paths
        with minimal length are returned.

        No graph compression is performed, so alternative minimal routes are
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

        records = []

        for marker in markers:
            if marker not in self:
                continue

            for scc in sccs:
                scc = frozenset(scc)

                for current_direction in (
                    ["downstream", "upstream"] if direction == "both" else [direction]
                ):

                    candidate_paths = []

                    for node in scc:

                        if current_direction == "downstream":

                            if not nx.has_path(
                                self,
                                node,
                                marker,
                            ):
                                continue

                            paths = nx.all_shortest_paths(
                                self,
                                source=node,
                                target=marker,
                            )

                        else:

                            if not nx.has_path(
                                self,
                                marker,
                                node,
                            ):
                                continue

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
        sign: Optional[CircuitSign] = None,
    ) -> list[tuple[list[str], int]]:
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

    def positive_circuits(self) -> list[list[str]]:
        """
        Return positive feedback circuits.

        Returns
        -------
        list of list of str
            Positive feedback circuits.
        """

        return [cycle for cycle, _ in self.circuits(sign=1)]

    def negative_circuits(self) -> list[list[str]]:
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
    ) -> Self:
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
            self.subgraph(feedback_nodes).copy(),
        )

    def compressed_graph(
        self,
        include_singleton_selfloops: bool = True,
        include_successors: bool = True,
        exclude_feedback_nodes: bool = False,
        min_size: int = 2,
        sep: str = "|",
    ) -> Self:
        """
        Return a compressed influence graph.

        Compression is applied sequentially through:
            1. feedback-induced graph extraction,
            2. structural family compression.

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
            from structural family compression.
        min_size: int (default: 2)
            Minimum family size required to collapse nodes.
        sep: str (default: "|")
            Separator used to build composite node names.

        Returns
        -------
        InfluenceGraph
            Compressed influence graph.

        See Also
        --------
        feedback_induced_graph
        family_compressed_graph
        """

        graph = self.feedback_induced_graph(
            include_singleton_selfloops=include_singleton_selfloops,
        )

        return graph.family_compressed_graph(
            include_successors=include_successors,
            exclude_feedback_nodes=exclude_feedback_nodes,
            min_size=min_size,
            sep=sep,
        )

    def to_pydot(
        self,
        program: str = "dot",
        edge_style: Optional[Callable[[Mapping[str, Any]], Mapping[str, Any]]] = None,
        **kwargs: Mapping[str, Any],
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
        **kwargs: Mapping[str, Any]
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
        **kwargs: Mapping[str, Any],
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
        **kwargs: Mapping[str, Any]
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

        svg = (
            self.to_pydot(
                program=program,
                edge_style=edge_style,
                **kwargs,
            )
            .create_svg()
            .decode()
        )

        display(SVG(svg))

    @staticmethod
    def _plain_graph(graph: Any) -> nx.MultiDiGraph:
        """
        Return a plain NetworkX MultiDiGraph copy.
        """

        source_graph = nx.MultiDiGraph(graph)
        plain_graph = nx.MultiDiGraph()
        plain_graph.graph.update(source_graph.graph)
        plain_graph.add_nodes_from(
            (node, data.copy()) for node, data in source_graph.nodes(data=True)
        )
        plain_graph.add_edges_from(
            (
                (source, target, key, data.copy())
                for source, target, key, data in source_graph.edges(
                    keys=True, data=True
                )
            ),
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
    def _cycle_edges(cycle: list[str]) -> list[tuple[str, str]]:
        """
        Return ordered directed edges from a cycle node list.
        """

        return list(zip(cycle, cycle[1:] + cycle[:1]))

    def _edge_signs(self, source: str, target: str) -> list[int]:
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
