#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from typing import (
    TYPE_CHECKING,
    Iterable,
    Iterator,
    List,
    Optional,
    Tuple,
    Union,
    cast,
)

from ..._compat import Literal
from ..._validation import _as_literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import Configuration, HypercubeLike
from ._asynchronous import (
    _AsynchronousTransition,
    _bdd_asynchronous_forward_chaining,
    _bdd_asynchronous_predecessors,
    _bdd_asynchronous_relation_partitions,
    _bdd_asynchronous_successors,
    _bdd_asynchronous_transition_partitions,
    _bdd_asynchronous_transition_sources,
)
from ._bdd import (
    _bdd_equivalence,
    _bdd_exclusive_or,
    _bdd_forward_quantification,
)
from ._dynamics import (
    _BottomBasinReduction,
    _compile_bitset_rules,
    _configuration_set_from_configuration_bits,
    _ConfigurationBits,
    _ExtendedComponentReduction,
    _ForwardBasinReduction,
    _ForwardReduction,
    _import_dd_backend,
    _restrict_transition_guided_reduction,
    _successor_configuration_bits,
    _terminal_strongly_connected_components,
    _transition_guided_reduction_size,
    _TransitionGuidedReduction,
    _validate_hypercube,
)
from ._general import _bdd_general_transition_partitions
from ._synchronous import _bdd_synchronous_transition_partitions

if TYPE_CHECKING:
    from ._network import BooleanNetwork
    from ._typing import (
        _BDDConfigurationSetNode,
        _BDDManager,
        _BDDNode,
        _BDDTransitionRelationNode,
    )

_SYNCHRONOUS_LARGE_NETWORK_SIZE = 90
_SYNCHRONOUS_UNATENESS_DAG_SIZE = 256
_SYMBOLIC_SCC_THRESHOLD = 128


class SymbolicTransitionSystem:
    """
    Represent finite-state Boolean dynamics for composable symbolic analyses.

    All symbolic configuration sets created by one transition system share its
    compiled dynamics. Sets from distinct systems are intentionally
    incompatible, even when both systems originate from the same Boolean
    network.

    The network is captured when the system is compiled. Later modifications
    to the original `BooleanNetwork` do not alter this transition system.

    Typical workflow
    ----------------
    ::

        BooleanNetwork
              |
              v
          symbolic()
              |
              v
        SymbolicTransitionSystem
              |
              | configurations(...)
              v
        SymbolicConfigurationSet
              |
              | post() / pre() / reachable() / ...
              |
              | configurations()
              v
        ConfigurationSet

    Examples
    --------
    The following asynchronous cascade reaches `00`, then `10`, then its
    unique fixed point `11`:

    >>> from bonesistools.logic.bn import BooleanNetwork
    >>> bn = BooleanNetwork({"A": 1, "B": "A"})
    >>> system = bn.symbolic(update="asynchronous")
    >>> initial = system.configurations({"A": 0, "B": 0})
    >>> reachable = initial.reachable()
    >>> attractors = reachable.terminal_sccs()
    >>> attractors[0].enumerate()
    ({'A': 1, 'B': 1},)
    >>> explicit = reachable.configurations()
    >>> explicit.enumerate()
    ({'A': 0, 'B': 0}, {'A': 1, 'B': 0}, {'A': 1, 'B': 1})

    Parameters
    ----------
    network: BooleanNetwork
        Closed Boolean network to compile.
    update: {"asynchronous", "synchronous", "general"} (default: "asynchronous")
        Finite-state update semantics.

    Raises
    ------
    ValueError
        If the network is not closed or `update` is invalid.
    """

    __slots__ = (
        "_and_exists",
        "_bdd",
        "_components",
        "_current_variable_set",
        "_current_variables",
        "_direct_forward_quantification",
        "_direct_transition_clusters",
        "_direct_transition_relation",
        "_disjunctive_transitions",
        "_forward_quantification",
        "_network",
        "_next_variables",
        "_rename_current_to_next",
        "_rename_next_to_current",
        "_transition_clusters",
        "_update",
    )

    def __init__(
        self,
        network: "BooleanNetwork",
        *,
        update: Literal[
            "asynchronous",
            "synchronous",
            "general",
        ] = "asynchronous",
    ) -> None:

        resolved_update = _as_literal(
            update,
            choices=("asynchronous", "synchronous", "general"),
            name="update",
        )
        network.validate()
        self._initialize(network.copy(), update=resolved_update)

    def __repr__(self) -> str:
        """Return a compact description of the compiled transition system."""

        return (
            f"SymbolicTransitionSystem(components={len(self._components)}, "
            f"update={self._update!r})"
        )

    def configurations(
        self,
        configurations: Union[
            HypercubeLike,
            ConfigurationSet,
            Iterable[HypercubeLike],
        ],
    ) -> "SymbolicConfigurationSet":
        """
        Create a symbolic configuration set from one or more configurations.

        Partial configurations represent all compatible complete
        configurations. A mapping is interpreted as one possibly partial
        configuration. A `ConfigurationSet` is encoded directly from its
        internal disjoint hypercubes. Any other iterable is interpreted as a
        union of possibly partial configurations.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> system.configurations({"A": 0}).enumerate()
        ({'A': 0, 'B': 0}, {'A': 0, 'B': 1})
        >>> states = system.configurations(
        ...     [{"A": 0, "B": 0}, {"A": 1, "B": 1}]
        ... )
        >>> states.enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 1})

        Parameters
        ----------
        configurations: HypercubeLike, ConfigurationSet or iterable
            Configuration, configuration set or union of Boolean subspaces.

        Returns
        -------
        SymbolicConfigurationSet
            Exact symbolic union of the supplied configurations.

        Raises
        ------
        TypeError
            If `configurations` is a symbolic set or cannot be interpreted as
            explicit configurations.
        ValueError
            If a configuration contains unknown components or a
            `ConfigurationSet` uses different components.
        """

        if isinstance(configurations, SymbolicConfigurationSet):
            raise TypeError(
                "symbolic configuration sets cannot be implicitly converted "
                "between transition systems"
            )
        if isinstance(configurations, ConfigurationSet):
            if set(configurations.components) != set(self._components):
                raise ValueError(
                    "configuration set components do not match transition "
                    "system components"
                )
            node = self._encode_configurations(configurations)
        elif isinstance(configurations, MappingABC):
            hypercube = _validate_hypercube(
                self._components,
                cast(HypercubeLike, configurations),
                name="configurations",
            )
            node = self._encode_hypercube(hypercube)
        else:
            node = self._bdd.false
            for configuration in configurations:
                hypercube = _validate_hypercube(
                    self._components,
                    configuration,
                    name="configurations",
                )
                node |= self._encode_hypercube(hypercube)

        return self._wrap(node)

    def empty(self) -> "SymbolicConfigurationSet":
        """
        Return the empty symbolic configuration set.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A"}).symbolic()
        >>> system.empty().is_empty()
        True

        Returns
        -------
        SymbolicConfigurationSet
            Empty set owned by this transition system.
        """

        return self._wrap(self._bdd.false)

    def universe(self) -> "SymbolicConfigurationSet":
        """
        Return all Boolean configurations over `components`.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> system.universe().enumerate()
        ({'A': 0, 'B': 0}, {'A': 0, 'B': 1}, {'A': 1, 'B': 0}, {'A': 1, 'B': 1})

        Returns
        -------
        SymbolicConfigurationSet
            Complete Boolean configuration space.
        """

        return self._wrap(self._bdd.true)

    def boolean_network(self) -> "BooleanNetwork":
        """
        Return the Boolean network captured at compilation.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> original = BooleanNetwork({"A": "~A"})
        >>> restored = original.symbolic().boolean_network()
        >>> restored == original
        True
        >>> restored is original
        False

        Returns
        -------
        BooleanNetwork
            New Boolean network instance that can be modified independently.
        """

        return self._network.copy()

    @property
    def update(self) -> Literal["asynchronous", "synchronous", "general"]:
        """Return the compiled update semantics."""

        return self._update

    @property
    def components(self) -> Tuple[str, ...]:
        """Return components in symbolic variable order."""

        return self._components

    def post(
        self,
        configurations: "SymbolicConfigurationSet",
    ) -> "SymbolicConfigurationSet":
        """
        Return successors after exactly one transition.

        Asynchronous and general transitions are non-reflexive, so a stable
        configuration has no successor. Synchronous dynamics return `f(x)`;
        consequently, a synchronous fixed point is its own successor.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A"}).symbolic()
        >>> fixed = system.configurations({"A": 1})
        >>> fixed.post().is_empty()
        True
        >>> synchronous = BooleanNetwork({"A": "A"}).symbolic(
        ...     update="synchronous"
        ... )
        >>> synchronous.configurations({"A": 1}).post().enumerate()
        ({'A': 1},)

        Parameters
        ----------
        configurations: SymbolicConfigurationSet
            Source configurations owned by this transition system.

        Returns
        -------
        SymbolicConfigurationSet
            Exact union of one-step successors.

        Raises
        ------
        TypeError
            If `configurations` is not a `SymbolicConfigurationSet`.
        ValueError
            If `configurations` belongs to another transition system.
        """

        node = self._configuration_node(configurations)
        return self._wrap(self._direct_successors(node))

    def pre(
        self,
        configurations: "SymbolicConfigurationSet",
    ) -> "SymbolicConfigurationSet":
        """
        Return predecessors after exactly one transition.

        Asynchronous and general transitions are non-reflexive. Synchronous
        dynamics retain the transition from a fixed point to itself.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic()
        >>> target = system.configurations({"A": 1})
        >>> system.pre(target).enumerate()
        ({'A': 0},)

        Parameters
        ----------
        configurations: SymbolicConfigurationSet
            Target configurations owned by this transition system.

        Returns
        -------
        SymbolicConfigurationSet
            Exact union of one-step predecessors.

        Raises
        ------
        TypeError
            If `configurations` is not a `SymbolicConfigurationSet`.
        ValueError
            If `configurations` belongs to another transition system.
        """

        node = self._configuration_node(configurations)
        return self._wrap(self._direct_predecessors(node))

    def reachable(
        self,
        initial: "SymbolicConfigurationSet",
        *,
        within: Optional["SymbolicConfigurationSet"] = None,
    ) -> "SymbolicConfigurationSet":
        """
        Return all configurations reachable from the initial set.

        This forward reachability closure includes the initial configurations
        and therefore uses paths of zero or more transitions. If `within` is
        supplied, initial configurations outside that region and transitions
        leaving it are excluded.

        Under asynchronous semantics, transition partitions are chained so
        configurations reached by one partition are immediately available to
        subsequent partitions [1].

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": 1, "B": "A"}).symbolic()
        >>> initial = system.configurations({"A": 0, "B": 0})
        >>> system.reachable(initial).enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 0}, {'A': 1, 'B': 1})

        Keeping `B = 0` excludes the final configuration:

        >>> region = system.configurations({"B": 0})
        >>> system.reachable(initial, within=region).enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 0})

        Parameters
        ----------
        initial: SymbolicConfigurationSet
            Initial configurations owned by this transition system.
        within: SymbolicConfigurationSet, optional
            Region in which paths must remain. If `None`, use the complete
            configuration space.

        Returns
        -------
        SymbolicConfigurationSet
            Configurations reachable through zero or more transitions.

        Raises
        ------
        TypeError
            If an argument is not a `SymbolicConfigurationSet`.
        ValueError
            If an argument belongs to another transition system.

        References
        ----------
        [1] Ciardo et al. (2006). The saturation algorithm for symbolic
        state-space exploration. International Journal on Software Tools for
        Technology Transfer, 8(1), 4-25.
        """

        initial_node = self._configuration_node(initial)
        within_node = self._optional_configuration_node(within)
        if within_node is not None:
            initial_node &= within_node
        return self._wrap(
            self._forward_reachable_states(initial_node, within=within_node)
        )

    def coreachable(
        self,
        target: "SymbolicConfigurationSet",
        *,
        within: Optional["SymbolicConfigurationSet"] = None,
    ) -> "SymbolicConfigurationSet":
        """
        Return all configurations from which the target set is reachable.

        This backward reachability closure includes the target configurations
        and therefore uses paths of zero or more transitions. If `within` is
        supplied, target configurations outside that region and transitions
        leaving it are excluded.

        Examples
        --------
        Both configurations of the Boolean switch can reach `A = 1`:

        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic()
        >>> target = system.configurations({"A": 1})
        >>> system.coreachable(target).enumerate()
        ({'A': 0}, {'A': 1})

        Parameters
        ----------
        target: SymbolicConfigurationSet
            Target configurations owned by this transition system.
        within: SymbolicConfigurationSet, optional
            Region in which paths must remain. If `None`, use the complete
            configuration space.

        Returns
        -------
        SymbolicConfigurationSet
            Configurations from which a target is reachable.

        Raises
        ------
        TypeError
            If an argument is not a `SymbolicConfigurationSet`.
        ValueError
            If an argument belongs to another transition system.
        """

        target_node = self._configuration_node(target)
        within_node = self._optional_configuration_node(within)
        if within_node is None:
            within_node = self._bdd.true
        else:
            target_node &= within_node
        return self._wrap(self._backward_reachable_states(target_node, within_node))

    def terminal_sccs(
        self,
        configurations: Optional["SymbolicConfigurationSet"] = None,
    ) -> Tuple["SymbolicConfigurationSet", ...]:
        """
        Return strongly connected components with no outgoing transition in a region.

        When no region is supplied, terminality is evaluated in the complete
        transition graph. Otherwise, it is evaluated in the subgraph induced
        by the supplied configurations.

        For a transition-closed region, these terminal components correspond
        to the attractors contained in that region.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic(
        ...     update="asynchronous"
        ... )
        >>> system.terminal_sccs()[0].enumerate()
        ({'A': 0}, {'A': 1})

        Notes
        -----
        Small regions are handled explicitly. Larger regions use the same
        semantics-specific symbolic selection as
        `BooleanNetwork.attractors()`.

        For transition-closed asynchronous regions, candidates are reduced
        with interleaved transition-guided reduction [1], then terminal
        components are extracted with the Xie-Beerel procedure [2]. Arbitrary
        non-closed regions retain generic induced-subgraph SCC semantics.

        Parameters
        ----------
        configurations: SymbolicConfigurationSet, optional
            Region whose induced transition graph is analyzed. If `None`, use
            the complete configuration space.

        Returns
        -------
        tuple of SymbolicConfigurationSet
            Exact terminal strongly connected components.

        Raises
        ------
        TypeError
            If `configurations` is not a `SymbolicConfigurationSet`.
        ValueError
            If `configurations` belongs to another transition system.

        References
        ----------
        [1] Benes et al. (2021). Computing bottom SCCs symbolically using
        transition guided reduction. International Conference on Computer
        Aided Verification, 505-528.

        [2] Xie and Beerel (2000). Implicit enumeration of strongly connected
        components and an application to formal verification. IEEE
        Transactions on Computer-Aided Design of Integrated Circuits and
        Systems, 19(10), 1225-1230.
        """

        if configurations is None:
            states = self._bdd.true
        else:
            states = self._configuration_node(configurations)

        components = self._terminal_sccs_for_region(
            states,
            transition_closed=True if states == self._bdd.true else None,
        )

        return tuple(self._wrap(component) for component in components)

    def _transition(
        self,
        source: HypercubeLike,
        target: HypercubeLike,
        *,
        quantifier: Literal["exists", "robust", "universal"],
    ) -> bool:
        """Test the one-step transition relation between two subspaces."""

        source_hypercube = _validate_hypercube(
            self._components,
            source,
            name="source",
        )
        target_hypercube = _validate_hypercube(
            self._components,
            target,
            name="target",
        )
        initial = self._encode_hypercube(source_hypercube)
        target_node = self._encode_hypercube(target_hypercube)

        if quantifier == "universal":
            self._ensure_next_variables()
            target_next = self._bdd.let(
                self._rename_current_to_next,
                target_node,
            )
            required_transitions = initial & target_next
            return (
                required_transitions & ~self._transition_relation() == self._bdd.false
            )

        predecessors = self._transition_predecessors(target_node)
        if quantifier == "exists":
            return initial & predecessors != self._bdd.false

        return initial & ~predecessors == self._bdd.false

    def _reachability(
        self,
        source: HypercubeLike,
        target: HypercubeLike,
        *,
        quantifier: Literal["exists", "robust", "universal"],
    ) -> bool:
        """Test whether initial configurations can reach a target subspace."""

        source_hypercube = _validate_hypercube(
            self._components,
            source,
            name="source",
        )
        target_hypercube = _validate_hypercube(
            self._components,
            target,
            name="target",
        )
        initial = self._encode_hypercube(source_hypercube)
        target_node = self._encode_hypercube(target_hypercube)
        reachable_from_initial = self._forward_reachable_states(initial)
        reachable_target = reachable_from_initial & target_node

        if reachable_target == self._bdd.false:
            return False

        if quantifier == "universal":
            return self._all_initial_configurations_reach_all_target_configurations(
                initial,
                target_hypercube,
                reachable_from_initial,
            )

        source_is_concrete = (
            source_hypercube.components == frozenset(self._components)
            and source_hypercube.is_fully_specified
        )
        if quantifier == "exists" or source_is_concrete:
            return True

        states_reaching_target = reachable_target

        while True:
            if self._matches_reachability_quantifier(
                initial,
                states_reaching_target,
                quantifier=quantifier,
            ):
                return True

            predecessors = self._predecessors(states_reaching_target)
            predecessors &= reachable_from_initial
            updated = states_reaching_target | predecessors
            if updated == states_reaching_target:
                return False

            states_reaching_target = updated

    def _reachable_attractors(
        self,
        initial_configurations: ConfigurationSet,
    ) -> Tuple[ConfigurationSet, ...]:
        """Return reachable terminal SCCs without decoding transient states."""

        initial = self._encode_configurations(initial_configurations)
        reachable = self._forward_reachable_states(initial)
        return tuple(
            self._decode_configuration_set(component)
            for component in self._terminal_sccs_for_region(
                reachable,
                transition_closed=True,
            )
        )

    def _terminal_sccs_for_region(
        self,
        states: _BDDConfigurationSetNode,
        *,
        transition_closed: Optional[bool] = None,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """Select one shared terminal-SCC algorithm for a symbolic region."""

        if states == self._bdd.false:
            return ()

        n_states = self._bdd.count(
            states,
            nvars=len(self._components),
        )
        if n_states <= _SYMBOLIC_SCC_THRESHOLD:
            return self._explicit_terminal_scc_nodes(states)

        if self._update == "general":
            return self._terminal_sccs(states)

        if transition_closed is None:
            transition_closed = self._is_transition_closed(states)
        if not transition_closed:
            return self._terminal_sccs(states)

        if self._update == "synchronous":
            return self._synchronous_terminal_sccs(states)

        candidates = self._transition_guided_reduction(states)
        return self._xie_beerel_terminal_sccs(candidates)

    def _is_transition_closed(self, states: _BDDConfigurationSetNode) -> bool:
        """Return whether every successor of a region remains in that region."""

        return self._successors(states) & ~states == self._bdd.false

    def _reachable_configurations(
        self,
        initial_configuration: Hypercube,
    ) -> Iterator[Configuration]:
        """Iterate over configurations in the symbolic reachable-state set."""

        initial = self._encode_hypercube(initial_configuration)
        reachable = self._forward_reachable_states(initial)

        def iterate() -> Iterator[Configuration]:
            for assignment in self._bdd.pick_iter(
                reachable,
                care_vars=self._current_variable_set,
            ):
                yield {
                    component: int(assignment[variable])
                    for component, variable in zip(
                        self._components,
                        self._current_variables,
                    )
                }

        return iterate()

    def _initialize(
        self,
        network: "BooleanNetwork",
        *,
        update: Literal["asynchronous", "synchronous", "general"],
    ) -> None:
        """Compile a validated network into the selected BDD backend."""

        bdd_module = _import_dd_backend()
        self._network = network
        self._bdd: _BDDManager = bdd_module.BDD()
        self._and_exists = getattr(bdd_module, "and_exists", None)
        self._update: Literal["asynchronous", "synchronous", "general"] = update
        self._components = tuple(network.keys())
        self._current_variables = tuple(
            f"x{index}" for index in range(len(self._components))
        )
        self._current_variable_set = set(self._current_variables)
        self._next_variables = tuple(
            f"y{index}" for index in range(len(self._components))
        )
        if update == "asynchronous":
            ordered_variables = self._current_variables
        else:
            ordered_variables = tuple(
                variable
                for pair in zip(self._current_variables, self._next_variables)
                for variable in pair
            )
        self._bdd.declare(*ordered_variables)
        self._disjunctive_transitions: Tuple[_AsynchronousTransition, ...] = ()
        self._transition_clusters: Tuple[_BDDTransitionRelationNode, ...] = ()
        self._forward_quantification: Tuple[Tuple[str, ...], ...] = ()
        if update in {"synchronous", "general"}:
            transition_partitions = (
                _bdd_synchronous_transition_partitions
                if update == "synchronous"
                else _bdd_general_transition_partitions
            )
            (
                self._transition_clusters,
                self._forward_quantification,
            ) = transition_partitions(
                network,
                bdd=self._bdd,
                components=self._components,
                current_vars=self._current_variables,
                next_vars=self._next_variables,
            )
        else:
            self._disjunctive_transitions = _bdd_asynchronous_transition_partitions(
                network,
                bdd=self._bdd,
                components=self._components,
                current_vars=self._current_variables,
            )
        self._rename_current_to_next = dict(
            zip(self._current_variables, self._next_variables)
        )
        self._rename_next_to_current = dict(
            zip(self._next_variables, self._current_variables)
        )
        self._direct_transition_clusters = None
        self._direct_forward_quantification = None
        self._direct_transition_relation = None

    def _wrap(
        self,
        node: _BDDConfigurationSetNode,
    ) -> "SymbolicConfigurationSet":
        """Wrap a BDD node without copying or validating it."""

        return SymbolicConfigurationSet._from_node(self, node)

    def _configuration_node(
        self,
        configurations: "SymbolicConfigurationSet",
    ) -> _BDDConfigurationSetNode:
        """Return the BDD node of a compatible symbolic configuration set."""

        if not isinstance(configurations, SymbolicConfigurationSet):
            raise TypeError(
                "expected a SymbolicConfigurationSet belonging to this "
                "transition system"
            )
        if configurations.system is not self:
            raise ValueError(
                "symbolic configuration sets belong to different transition systems"
            )
        return configurations._node

    def _optional_configuration_node(
        self,
        configurations: Optional["SymbolicConfigurationSet"],
    ) -> Optional[_BDDConfigurationSetNode]:
        """Return a compatible node or None for an omitted symbolic set."""

        if configurations is None:
            return None
        return self._configuration_node(configurations)

    def _direct_successors(
        self,
        configurations: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """
        Return exact one-step successors for the public `post()` operation.

        General reachability uses identity transitions internally to simplify
        its partitioned relation. This direct variant adds the constraint that
        at least one component changes, as required by the public semantics.
        """

        if self._update != "general":
            return self._successors(configurations)

        clusters, forward_quantification = self._general_direct_partitions()
        successors = configurations
        for cluster, quantified_variables in zip(
            clusters,
            forward_quantification,
        ):
            if quantified_variables:
                successors = self._relational_product(
                    successors,
                    cluster,
                    quantified_variables,
                )
            else:
                successors &= cluster
        return self._bdd.let(
            self._rename_next_to_current,
            successors,
        )

    def _direct_predecessors(
        self,
        configurations: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Return exact one-step predecessors without a monolithic relation."""

        if self._update != "general":
            return self._predecessors(configurations)

        configurations_next = self._bdd.let(
            self._rename_current_to_next,
            configurations,
        )
        changed = self._general_direct_partitions()[0][-1]
        predecessors = configurations_next & changed
        for cluster, next_variable in zip(
            self._transition_clusters,
            self._next_variables,
        ):
            predecessors = self._relational_product(
                predecessors,
                cluster,
                (next_variable,),
            )
        return predecessors

    def _general_direct_partitions(
        self,
    ) -> Tuple[
        Tuple[_BDDTransitionRelationNode, ...],
        Tuple[Tuple[str, ...], ...],
    ]:
        """Return lazily compiled partitions for non-reflexive general steps."""

        if (
            self._direct_transition_clusters is None
            or self._direct_forward_quantification is None
        ):
            clusters = self._transition_clusters + (self._general_changed_partition(),)
            self._direct_transition_clusters = clusters
            self._direct_forward_quantification = _bdd_forward_quantification(
                self._bdd,
                clusters,
                self._current_variables,
            )

        return (
            self._direct_transition_clusters,
            self._direct_forward_quantification,
        )

    def _general_changed_partition(self) -> _BDDTransitionRelationNode:
        """Return the lazy BDD constraint requiring one changed component."""

        if self._direct_transition_clusters is not None:
            return self._direct_transition_clusters[-1]

        changed = self._bdd.false
        for current_variable, next_variable in zip(
            self._current_variables,
            self._next_variables,
        ):
            changed |= _bdd_exclusive_or(
                self._bdd.var(current_variable),
                self._bdd.var(next_variable),
            )
        return changed

    def _encode_hypercube(
        self,
        hypercube: Hypercube,
    ) -> _BDDConfigurationSetNode:
        """Encode one validated Boolean hypercube."""

        return self._encode_hypercube_variables(
            hypercube,
            self._current_variables,
        )

    def _encode_hypercube_variables(
        self,
        hypercube: Hypercube,
        variables: Tuple[str, ...],
    ) -> _BDDConfigurationSetNode:
        """Encode one hypercube over a selected BDD variable family."""

        encoded = self._bdd.true
        for component, variable_name in zip(self._components, variables):
            if component not in hypercube:
                continue

            value = hypercube[component]
            if value.is_fixed:
                variable = self._bdd.var(variable_name)
                encoded &= variable if value.value else ~variable

        return encoded

    def _encode_configurations(
        self,
        configurations: ConfigurationSet,
    ) -> _BDDConfigurationSetNode:
        """Encode a ConfigurationSet as a BDD over current variables."""

        source_indices = {
            component: index
            for index, component in enumerate(configurations.components)
        }
        states = self._bdd.false
        for fixed_mask, value_mask in configurations._iter_encoded_hypercubes():
            encoded = self._bdd.true
            for component, variable_name in zip(
                self._components,
                self._current_variables,
            ):
                bit = 1 << source_indices[component]
                if not fixed_mask & bit:
                    continue

                variable = self._bdd.var(variable_name)
                encoded &= variable if value_mask & bit else ~variable
            states |= encoded

        return states

    def _forward_reachable_states(
        self,
        initial: _BDDConfigurationSetNode,
        within: Optional[_BDDConfigurationSetNode] = None,
    ) -> _BDDConfigurationSetNode:
        """
        Return the symbolic forward-reachability closure of an initial set.

        Successors are accumulated until no transition discovers a new
        configuration. If `within` is provided, successors outside that region
        are excluded. Asynchronous semantics use transition-partition chaining
        [1], refined with one workset per partition.

        References
        ----------
        [1] Ciardo et al. (2006). The saturation algorithm for symbolic
        state-space exploration. International Journal on Software Tools for
        Technology Transfer, 8(1), 4-25.
        """

        if self._update == "asynchronous":
            return _bdd_asynchronous_forward_chaining(
                self._bdd,
                initial,
                self._disjunctive_transitions,
                within=within,
            )

        reachable = initial
        frontier = initial

        while frontier != self._bdd.false:
            successors = self._successors(frontier)
            if within is not None:
                successors &= within
            new_frontier = successors & ~reachable
            if new_frontier == self._bdd.false:
                break

            reachable |= new_frontier
            frontier = new_frontier

        return reachable

    def _backward_reachable_states(
        self,
        initial: _BDDConfigurationSetNode,
        within: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """
        Return the symbolic backward-reachability closure within a region.

        Predecessors contained in `within` are accumulated until no transition
        discovers a new configuration.
        """

        reachable = initial
        frontier = initial

        while frontier != self._bdd.false:
            predecessors = self._predecessors(frontier) & within
            new_frontier = predecessors & ~reachable
            if new_frontier == self._bdd.false:
                break

            reachable |= new_frontier
            frontier = new_frontier

        return reachable

    def _synchronous_terminal_sccs(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """Select a conservative synchronous terminal-SCC algorithm."""

        if len(self._components) >= _SYNCHRONOUS_LARGE_NETWORK_SIZE:
            return self._terminal_sccs(states)

        transition_dag_size = sum(
            transition.dag_size for transition in self._transition_clusters
        )
        if (
            transition_dag_size > _SYNCHRONOUS_UNATENESS_DAG_SIZE
            and self._has_many_non_unate_synchronous_influences()
        ):
            return self._terminal_sccs(states)
        return self._synchronous_terminal_cycles(states)

    def _has_many_non_unate_synchronous_influences(self) -> bool:
        """Test whether over one third of rule influences are non-unate."""

        n_influences = 0
        n_non_unate = 0
        for transition, next_variable in zip(
            self._transition_clusters,
            self._next_variables,
        ):
            rule = self._bdd.let({next_variable: True}, transition)
            for current_variable in self._bdd.support(rule):
                low = self._bdd.let({current_variable: False}, rule)
                high = self._bdd.let({current_variable: True}, rule)
                positive = high & ~low != self._bdd.false
                negative = low & ~high != self._bdd.false
                n_influences += 1
                n_non_unate += positive and negative

        return 3 * n_non_unate > n_influences

    def _synchronous_terminal_cycles(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """Extract cycles from a deterministic synchronous transition system."""

        remaining = self._synchronous_recurrent_states(states)
        cycles = []

        while remaining != self._bdd.false:
            seed = self._pick_state(remaining)
            cycle = seed
            successor = self._successors(seed)

            while successor & cycle == self._bdd.false:
                cycle |= successor
                successor = self._successors(successor)

            cycles.append(cycle)
            remaining &= ~cycle

        return tuple(cycles)

    def _synchronous_recurrent_states(
        self,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Remove transient states until only synchronous cycles remain."""

        recurrent = states

        while True:
            updated = recurrent & self._successors(recurrent)
            if updated == recurrent:
                return recurrent
            recurrent = updated

    def _transition_guided_reduction(
        self,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """
        Remove states that cannot belong to asynchronous terminal SCCs.

        This implements interleaved transition-guided reduction [1]. Each
        asynchronous component update is treated as one transition label.

        References
        ----------
        [1] Benes et al. (2021). Computing bottom SCCs symbolically using
        transition guided reduction. CAV 2021, 505-528.
        """

        remaining = states
        reductions: List[_TransitionGuidedReduction] = []
        for transition_index, transition in enumerate(self._disjunctive_transitions):
            forward = self._transition_sources(transition, remaining)
            if forward != self._bdd.false:
                reductions.append(
                    _ForwardReduction(
                        transition_index=transition_index,
                        forward=forward,
                    )
                )

        to_discard = self._bdd.false
        while reductions:
            if to_discard != self._bdd.false:
                remaining &= ~to_discard
                for reduction in reductions:
                    _restrict_transition_guided_reduction(reduction, remaining)
                to_discard = self._bdd.false

            if not isinstance(
                reductions[-1],
                (_ForwardBasinReduction, _BottomBasinReduction),
            ):
                reductions.sort(
                    key=_transition_guided_reduction_size,
                    reverse=True,
                )

            reduction = reductions[-1]
            if isinstance(reduction, _ForwardReduction):
                successors = self._saturation_successors(
                    reduction.forward,
                    within=remaining,
                )
                if successors != self._bdd.false:
                    reduction.forward |= successors
                    continue

                reductions.pop()
                transition = self._disjunctive_transitions[reduction.transition_index]
                reductions.append(
                    _ExtendedComponentReduction(
                        transition_index=reduction.transition_index,
                        forward=reduction.forward,
                        extended_component=self._transition_sources(
                            transition,
                            remaining,
                        ),
                    )
                )
                if remaining & ~reduction.forward != self._bdd.false:
                    reductions.append(
                        _ForwardBasinReduction(
                            forward=reduction.forward,
                            basin=reduction.forward,
                        )
                    )
                continue

            if isinstance(reduction, _ExtendedComponentReduction):
                predecessors = self._saturation_predecessors(
                    reduction.extended_component,
                    within=reduction.forward,
                )
                if predecessors != self._bdd.false:
                    reduction.extended_component |= predecessors
                    continue

                reductions.pop()
                bottom = reduction.forward & ~reduction.extended_component
                transition = self._disjunctive_transitions[reduction.transition_index]
                escaping = self._transition_sources_outside(
                    transition,
                    remaining,
                )
                if bottom != self._bdd.false or escaping != self._bdd.false:
                    reductions.append(
                        _BottomBasinReduction(
                            bottom=bottom,
                            basin=bottom | escaping,
                        )
                    )
                continue

            if isinstance(reduction, _ForwardBasinReduction):
                predecessors = self._saturation_predecessors(
                    reduction.basin,
                    within=remaining,
                )
                if predecessors != self._bdd.false:
                    reduction.basin |= predecessors
                    continue

                reductions.pop()
                to_discard = reduction.basin & ~reduction.forward
                continue

            predecessors = self._saturation_predecessors(
                reduction.basin,
                within=remaining,
            )
            if predecessors != self._bdd.false:
                reduction.basin |= predecessors
                continue

            reductions.pop()
            to_discard = reduction.basin & ~reduction.bottom

        if to_discard != self._bdd.false:
            remaining &= ~to_discard

        return remaining

    def _xie_beerel_terminal_sccs(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """
        Extract global terminal SCCs from a reduced candidate state set.

        The search follows the symbolic bottom-SCC procedure of Xie and Beerel
        [1]. Successors are evaluated against the original transition relation
        because transition-guided reduction does not return a closed set.

        References
        ----------
        [1] Xie and Beerel (2000). Implicit enumeration of strongly connected
        components and an application to formal verification. IEEE
        Transactions on Computer-Aided Design of Integrated Circuits and
        Systems, 19(10), 1225-1230.
        """

        remaining = states
        pivot_hint = self._bdd.false
        terminal_components = []

        while remaining != self._bdd.false:
            hinted_states = pivot_hint & remaining
            pivot_hint = self._bdd.false
            pivot = self._pick_state(
                hinted_states if hinted_states != self._bdd.false else remaining
            )
            basin = self._saturated_backward_reachable_states(
                pivot,
                within=remaining,
            )
            component = pivot

            while True:
                successors = self._saturation_successors(component)
                if successors == self._bdd.false:
                    terminal_components.append(component)
                    break

                component |= successors
                escaped = successors & ~basin
                if escaped != self._bdd.false:
                    pivot_hint = escaped
                    break

            remaining &= ~basin

        return tuple(terminal_components)

    def _saturated_backward_reachable_states(
        self,
        initial: _BDDConfigurationSetNode,
        *,
        within: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Compute backward reachability one transition label at a time."""

        reachable = initial
        while True:
            predecessors = self._saturation_predecessors(
                reachable,
                within=within,
            )
            if predecessors == self._bdd.false:
                return reachable
            reachable |= predecessors

    def _saturation_successors(
        self,
        states: _BDDConfigurationSetNode,
        within: Optional[_BDDConfigurationSetNode] = None,
    ) -> _BDDConfigurationSetNode:
        """Return new successors for the first productive transition label."""

        for transition in reversed(self._disjunctive_transitions):
            successors = self._partition_successors(states, transition)
            if within is not None:
                successors &= within
            successors &= ~states
            if successors != self._bdd.false:
                return successors

        return self._bdd.false

    def _saturation_predecessors(
        self,
        states: _BDDConfigurationSetNode,
        *,
        within: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Return new predecessors for the first productive transition label."""

        for transition in reversed(self._disjunctive_transitions):
            predecessors = self._partition_predecessors(states, transition)
            predecessors &= within & ~states
            if predecessors != self._bdd.false:
                return predecessors

        return self._bdd.false

    def _transition_sources(
        self,
        transition: _AsynchronousTransition,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Return states enabling one labelled asynchronous transition."""

        return _bdd_asynchronous_transition_sources(
            states,
            transition,
        )

    def _transition_sources_outside(
        self,
        transition: _AsynchronousTransition,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Return sources whose labelled transition leaves a state set."""

        return states & _bdd_asynchronous_predecessors(
            self._bdd,
            ~states,
            transition,
        )

    def _partition_successors(
        self,
        states: _BDDConfigurationSetNode,
        transition: _AsynchronousTransition,
    ) -> _BDDConfigurationSetNode:
        """Return successors through one asynchronous transition partition."""

        return _bdd_asynchronous_successors(
            self._bdd,
            states,
            transition,
        )

    def _partition_predecessors(
        self,
        states: _BDDConfigurationSetNode,
        transition: _AsynchronousTransition,
    ) -> _BDDConfigurationSetNode:
        """Return predecessors through one asynchronous transition partition."""

        return _bdd_asynchronous_predecessors(
            self._bdd,
            states,
            transition,
        )

    def _terminal_sccs(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """Extract terminal SCCs while keeping state regions symbolic."""

        remaining = states
        terminal_components = []
        while remaining != self._bdd.false:
            component = self._terminal_scc(remaining)
            terminal_components.append(component)
            basin = self._backward_reachable_states(component, remaining)
            remaining &= ~basin

        return tuple(terminal_components)

    def _terminal_scc(
        self,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Follow the symbolic SCC graph until reaching one terminal SCC."""

        region = states
        seed = self._pick_state(region)

        while True:
            forward = self._forward_reachable_states(seed, within=region)
            backward = self._backward_reachable_states(seed, forward)
            component = forward & backward
            successors = self._successors(component) & region & ~component
            if successors == self._bdd.false:
                return component

            region = forward & ~backward
            seed = self._pick_state(successors)

    def _explicit_terminal_scc_nodes(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_BDDConfigurationSetNode, ...]:
        """Compute terminal SCCs of a small symbolic region explicitly."""

        configuration_bits = frozenset(self._decode_configuration_bits(states))
        compiled_rules = _compile_bitset_rules(self._network, self._components)
        successors = {
            state: tuple(
                successor
                for successor in _successor_configuration_bits(
                    compiled_rules,
                    state,
                    update=self._update,
                    n_components=len(self._components),
                )
                if successor in configuration_bits
            )
            for state in configuration_bits
        }
        return tuple(
            self._encode_configurations(
                _configuration_set_from_configuration_bits(
                    component,
                    self._components,
                )
            )
            for component in _terminal_strongly_connected_components(successors)
        )

    def _successors(
        self,
        configurations: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """
        Return closure-oriented successors under the configured semantics.

        A transition system is compiled for one semantics, stored in
        `self._update`, so the semantics is not passed to each operation.
        Asynchronous dynamics take the union of each component's local update;
        synchronous and general dynamics apply their compiled transition
        partitions. General dynamics may retain identity transitions because
        they do not change reachability closures or terminal SCCs. The public
        `post()` operation uses `_direct_successors()` to exclude them.
        """

        if self._update == "asynchronous":
            successors = self._bdd.false
            for transition in self._disjunctive_transitions:
                successors |= _bdd_asynchronous_successors(
                    self._bdd,
                    configurations,
                    transition,
                )
            return successors
        else:
            successors = configurations
            for cluster, quantified_variables in zip(
                self._transition_clusters,
                self._forward_quantification,
            ):
                if quantified_variables:
                    successors = self._relational_product(
                        successors,
                        cluster,
                        quantified_variables,
                    )
                else:
                    successors &= cluster
        return self._bdd.let(
            self._rename_next_to_current,
            successors,
        )

    def _predecessors(
        self,
        configurations: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """
        Return symbolic predecessors under the configured update semantics.

        A transition system is compiled for one semantics, stored in
        `self._update`, so the semantics is not passed to each operation.
        Asynchronous dynamics take the union of each component's local update;
        synchronous and general dynamics apply their compiled transition
        partitions.
        """

        if self._update == "asynchronous":
            predecessors = self._bdd.false
            for transition in self._disjunctive_transitions:
                predecessors |= _bdd_asynchronous_predecessors(
                    self._bdd,
                    configurations,
                    transition,
                )
            return predecessors

        configurations_next = self._bdd.let(
            self._rename_current_to_next,
            configurations,
        )
        predecessors = configurations_next
        for cluster, next_variable in zip(
            self._transition_clusters,
            self._next_variables,
        ):
            predecessors = self._relational_product(
                predecessors,
                cluster,
                (next_variable,),
            )
        return predecessors

    def _transition_predecessors(
        self,
        configurations: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Return sources with a direct transition into a configuration set."""

        if self._update != "general":
            return self._predecessors(configurations)

        configurations_next = self._bdd.let(
            self._rename_current_to_next,
            configurations,
        )
        return self._bdd.exist(
            self._next_variables,
            self._transition_relation() & configurations_next,
        )

    def _transition_relation(self) -> _BDDTransitionRelationNode:
        """Return the exact, non-reflexive relation where required."""

        if self._direct_transition_relation is not None:
            return self._direct_transition_relation

        if self._update == "asynchronous":
            self._ensure_next_variables()
            relation = self._bdd.false
            for transition in _bdd_asynchronous_relation_partitions(
                self._network,
                bdd=self._bdd,
                components=self._components,
                current_vars=self._current_variables,
                next_vars=self._next_variables,
                functional_transitions=self._disjunctive_transitions,
            ):
                relation |= transition
        else:
            relation = self._bdd.true
            for cluster in self._transition_clusters:
                relation &= cluster

            if self._update == "general":
                changed = self._bdd.false
                for current_variable, next_variable in zip(
                    self._current_variables,
                    self._next_variables,
                ):
                    changed |= _bdd_exclusive_or(
                        self._bdd.var(current_variable),
                        self._bdd.var(next_variable),
                    )
                relation &= changed

        self._direct_transition_relation = relation
        return relation

    def _ensure_next_variables(self) -> None:
        """Declare next-state variables when a full relation is requested."""

        undeclared = tuple(
            variable
            for variable in self._next_variables
            if variable not in self._bdd.vars
        )
        if undeclared:
            self._bdd.declare(*undeclared)

    def _relational_product(
        self,
        left: _BDDNode,
        right: _BDDNode,
        quantified_variables: Tuple[str, ...],
    ) -> _BDDNode:
        """Conjoin two BDDs while existentially quantifying variables."""

        if self._and_exists is not None:
            return self._and_exists(left, right, quantified_variables)

        return self._bdd.exist(quantified_variables, left & right)

    def _all_initial_configurations_reach_all_target_configurations(
        self,
        initial: _BDDConfigurationSetNode,
        target_hypercube: Hypercube,
        reachable_from_initial: _BDDConfigurationSetNode,
    ) -> bool:
        """Test universal reachability without enumerating target states."""

        target_variables = tuple(f"z{index}" for index in range(len(self._components)))
        undeclared_variables = tuple(
            variable for variable in target_variables if variable not in self._bdd.vars
        )
        if undeclared_variables:
            self._bdd.declare(*undeclared_variables)

        target_configurations = self._encode_hypercube_variables(
            target_hypercube,
            target_variables,
        )
        current_matches_target = target_configurations
        for current_variable, target_variable in zip(
            self._current_variables,
            target_variables,
        ):
            current_matches_target &= _bdd_equivalence(
                self._bdd.var(current_variable),
                self._bdd.var(target_variable),
            )

        required_pairs = initial & target_configurations
        states_reaching_targets = current_matches_target & reachable_from_initial

        while True:
            if required_pairs & ~states_reaching_targets == self._bdd.false:
                return True

            predecessors = self._predecessors(states_reaching_targets)
            predecessors &= reachable_from_initial
            predecessors &= target_configurations
            updated = states_reaching_targets | predecessors
            if updated == states_reaching_targets:
                return False

            states_reaching_targets = updated

    def _matches_reachability_quantifier(
        self,
        initial: _BDDConfigurationSetNode,
        states_reaching_target: _BDDConfigurationSetNode,
        *,
        quantifier: Literal["exists", "robust"],
    ) -> bool:
        """Test whether the current backward closure answers the query."""

        if quantifier == "exists":
            return initial & states_reaching_target != self._bdd.false

        return initial & ~states_reaching_target == self._bdd.false

    def _pick_state(
        self,
        states: _BDDConfigurationSetNode,
    ) -> _BDDConfigurationSetNode:
        """Select one concrete state from a non-empty symbolic state set."""

        assignment = next(
            iter(
                self._bdd.pick_iter(
                    states,
                    care_vars=self._current_variable_set,
                )
            )
        )
        selected = self._bdd.true
        for variable_name in self._current_variables:
            variable = self._bdd.var(variable_name)
            selected &= variable if assignment[variable_name] else ~variable

        return selected

    def _decode_configuration_set(
        self,
        states: _BDDConfigurationSetNode,
    ) -> ConfigurationSet:
        """Decode a symbolic state set into exact disjoint hypercubes."""

        hypercubes = []
        for assignment in self._bdd.pick_iter(states):
            fixed_mask = 0
            value_mask = 0
            for index, variable_name in enumerate(self._current_variables):
                if variable_name not in assignment:
                    continue

                bit = 1 << index
                fixed_mask |= bit
                if assignment[variable_name]:
                    value_mask |= bit
            hypercubes.append((fixed_mask, value_mask))

        return ConfigurationSet._from_encoded_hypercubes(
            self._components,
            hypercubes,
        )

    def _decode_configuration_bits(
        self,
        states: _BDDConfigurationSetNode,
    ) -> Tuple[_ConfigurationBits, ...]:
        """Decode a symbolic state set into explicit bitsets."""

        decoded = []
        for assignment in self._bdd.pick_iter(
            states,
            care_vars=self._current_variable_set,
        ):
            configuration_bits = 0
            for index, variable_name in enumerate(self._current_variables):
                if assignment[variable_name]:
                    configuration_bits |= 1 << index
            decoded.append(configuration_bits)

        return tuple(decoded)

    @classmethod
    def _from_validated_network(
        cls,
        network: "BooleanNetwork",
        *,
        update: Literal["asynchronous", "synchronous", "general"],
    ) -> "SymbolicTransitionSystem":
        """Compile a validated network without copying it."""

        transition_system = cls.__new__(cls)
        transition_system._initialize(network, update=update)
        return transition_system


class SymbolicConfigurationSet:
    """
    Immutable symbolic set of configurations tied to one transition system.

    Instances are created by `SymbolicTransitionSystem.configurations()`,
    `SymbolicTransitionSystem.empty()` and `SymbolicTransitionSystem.universe()`.
    Set operations remain symbolic and never materialize concrete
    configurations implicitly. Symbolic configuration sets behave like
    immutable mathematical sets.

    Examples
    --------
    Here, `left` represents `A = 0` and `right` represents `B = 1`:

    >>> from bonesistools.logic.bn import BooleanNetwork
    >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
    >>> left = system.configurations({"A": 0})
    >>> right = system.configurations({"B": 1})
    >>> (left & right).enumerate()
    ({'A': 0, 'B': 1},)
    >>> (left | right).enumerate()
    ({'A': 0, 'B': 0}, {'A': 0, 'B': 1}, {'A': 1, 'B': 1})
    >>> (~left).enumerate()
    ({'A': 1, 'B': 0}, {'A': 1, 'B': 1})
    >>> {"A": 0} in left
    True

    For a partial configuration, membership means that every compatible
    complete configuration belongs to the symbolic set.

    Supported operators
    -------------------
    =========================  ==========================================
    Expression                 Meaning
    =========================  ==========================================
    `left | right`             Union
    `left & right`             Intersection
    `left - right`             Difference
    `left ^ right`             Symmetric difference
    `~left`                    Complement
    `left == right`            Equality
    `left <= right`            Inclusion
    `left < right`             Strict inclusion
    `left >= right`            Reverse inclusion
    `left > right`             Strict reverse inclusion
    `configuration in left`    Full configuration or subspace containment
    =========================  ==========================================

    Notes
    -----
    Iteration and `enumerate()` may require exponentially many configurations.
    Whenever possible, prefer symbolic operations and `count()` over explicit
    materialization.
    """

    __slots__ = ("_node", "_system")

    def __init__(self) -> None:
        raise TypeError(
            "SymbolicConfigurationSet objects must be created by a "
            "SymbolicTransitionSystem"
        )

    def __contains__(self, configuration: object) -> bool:
        """
        Test whether a configuration or subspace is fully represented.

        Invalid or unsupported objects are treated as absent.
        """

        try:
            return self.contains(cast(HypercubeLike, configuration))
        except (TypeError, ValueError):
            return False

    def __repr__(self) -> str:
        """Return a compact representation without counting configurations."""

        return (
            "SymbolicConfigurationSet("
            f"components={len(self._system.components)}, "
            f"update={self._system.update!r})"
        )

    def __iter__(self) -> Iterator[Configuration]:
        """Iterate lazily over complete configurations."""

        for assignment in self._system._bdd.pick_iter(
            self._node,
            care_vars=self._system._current_variables,
        ):
            yield {
                component: int(assignment[variable])
                for component, variable in zip(
                    self._system.components,
                    self._system._current_variables,
                )
            }

    def __or__(self, other: object) -> "SymbolicConfigurationSet":
        """Return the symbolic union with a compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        return self._system._wrap(self._node | self._compatible_node(other))

    def __and__(self, other: object) -> "SymbolicConfigurationSet":
        """Return the symbolic intersection with a compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        return self._system._wrap(self._node & self._compatible_node(other))

    def __sub__(self, other: object) -> "SymbolicConfigurationSet":
        """Return the symbolic set difference with a compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        return self._system._wrap(self._node & ~self._compatible_node(other))

    def __xor__(self, other: object) -> "SymbolicConfigurationSet":
        """Return the symmetric difference with a compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        other_node = self._compatible_node(other)
        return self._system._wrap(
            (self._node & ~other_node) | (~self._node & other_node)
        )

    def __invert__(self) -> "SymbolicConfigurationSet":
        """Return the complement in the transition-system universe."""

        return self._system._wrap(~self._node)

    def __eq__(self, other: object) -> bool:
        """Test equality within the same symbolic transition system."""

        return (
            isinstance(other, SymbolicConfigurationSet)
            and self._system is other._system
            and self._node == other._node
        )

    def __le__(self, other: object) -> bool:
        """Test whether this set is included in another compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        return self._node & ~self._compatible_node(other) == self._system._bdd.false

    def __lt__(self, other: object) -> bool:
        """Test strict inclusion in another compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        other_node = self._compatible_node(other)
        return (
            self._node != other_node
            and self._node & ~other_node == self._system._bdd.false
        )

    def __ge__(self, other: object) -> bool:
        """Test whether this set includes another compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        other_node = self._compatible_node(other)
        return other_node & ~self._node == self._system._bdd.false

    def __gt__(self, other: object) -> bool:
        """Test whether this set strictly includes another compatible set."""

        if not isinstance(other, SymbolicConfigurationSet):
            return NotImplemented
        other_node = self._compatible_node(other)
        return (
            self._node != other_node
            and other_node & ~self._node == self._system._bdd.false
        )

    def configurations(self) -> ConfigurationSet:
        """
        Materialize this symbolic set as an explicit ConfigurationSet.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> symbolic = system.configurations({"A": 0})
        >>> explicit = symbolic.configurations()
        >>> explicit.enumerate()
        ({'A': 0, 'B': 0}, {'A': 0, 'B': 1})

        Returns
        -------
        ConfigurationSet
            Exact materialized configuration set.

        Notes
        -----
        The conversion preserves compact hypercubes when possible, but some
        symbolic sets may require many hypercubes to materialize.
        """

        return self._system._decode_configuration_set(self._node)

    @property
    def system(self) -> SymbolicTransitionSystem:
        """Return the transition system owning this symbolic set."""

        return self._system

    @property
    def components(self) -> Tuple[str, ...]:
        """Return ordered Boolean components."""

        return self._system.components

    def is_empty(self) -> bool:
        """
        Test whether no configuration is represented.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A"}).symbolic()
        >>> system.empty().is_empty()
        True

        Returns
        -------
        bool
            Whether this set is empty.
        """

        return self._node == self._system._bdd.false

    def is_universe(self) -> bool:
        """
        Test whether every configuration is represented.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A"}).symbolic()
        >>> system.universe().is_universe()
        True

        Returns
        -------
        bool
            Whether this set is the complete configuration space.
        """

        return self._node == self._system._bdd.true

    def count(self) -> int:
        """
        Return the exact number of represented configurations.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> system.configurations({"A": 0}).count()
        2

        Returns
        -------
        int
            Number of complete configurations in this set.
        """

        return int(
            self._system._bdd.count(
                self._node,
                nvars=len(self._system.components),
            )
        )

    def contains(self, configuration: HypercubeLike) -> bool:
        """
        Test whether a configuration or Boolean subspace is fully represented.

        For a partial configuration, this returns `True` only if every
        compatible complete configuration belongs to the set. Use
        `intersects()` to test whether at least one compatible configuration is
        shared instead.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> states = system.configurations({"A": 0})
        >>> states.contains({"A": 0, "B": 1})
        True
        >>> states.contains({"B": 1})
        False

        Parameters
        ----------
        configuration: HypercubeLike
            Complete or partially defined configuration to test.

        Returns
        -------
        bool
            Whether the represented subspace is fully included in this set.

        Raises
        ------
        TypeError
            If `configuration` cannot be interpreted as a Boolean subspace.
        ValueError
            If `configuration` contains an unknown component or invalid value.
        """

        hypercube = _validate_hypercube(
            self._system.components,
            configuration,
            name="configuration",
        )
        encoded = self._system._encode_hypercube(hypercube)
        return encoded & ~self._node == self._system._bdd.false

    def intersects(
        self,
        configurations: Union[HypercubeLike, "SymbolicConfigurationSet"],
    ) -> bool:
        """
        Test whether at least one configuration is shared with another set or subspace.

        Unlike `contains()`, this requires only one common complete
        configuration.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> states = system.configurations({"A": 0})
        >>> states.intersects({"B": 1})
        True
        >>> states.intersects({"A": 1})
        False

        Parameters
        ----------
        configurations: HypercubeLike or SymbolicConfigurationSet
            Boolean subspace or symbolic set to compare.

        Returns
        -------
        bool
            Whether the intersection is non-empty.

        Raises
        ------
        TypeError
            If `configurations` cannot be interpreted as a Boolean subspace.
        ValueError
            If a subspace is invalid or a symbolic set belongs to another
            transition system.
        """

        if isinstance(configurations, SymbolicConfigurationSet):
            node = self._compatible_node(configurations)
        else:
            hypercube = _validate_hypercube(
                self._system.components,
                configurations,
                name="configurations",
            )
            node = self._system._encode_hypercube(hypercube)
        return self._node & node != self._system._bdd.false

    def enumerate(self) -> Tuple[Configuration, ...]:
        """
        Return all represented complete configurations.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> system.configurations({"A": 0}).enumerate()
        ({'A': 0, 'B': 0}, {'A': 0, 'B': 1})

        Returns
        -------
        tuple of Configuration
            Materialized complete configurations.

        Notes
        -----
        The result size may be exponential in the number of components.
        """

        return tuple(self)

    def pick(self) -> Configuration:
        """
        Return one represented configuration without uniformity guarantees.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "A", "B": "B"}).symbolic()
        >>> states = system.configurations({"A": 0})
        >>> states.pick() in states
        True

        Returns
        -------
        Configuration
            One complete configuration contained in this set.

        Raises
        ------
        ValueError
            If the symbolic set is empty.
        """

        try:
            assignment = next(
                self._system._bdd.pick_iter(
                    self._node,
                    care_vars=self._system._current_variables,
                )
            )
        except StopIteration as error:
            raise ValueError("cannot pick from an empty symbolic set") from error

        return {
            component: int(assignment[variable])
            for component, variable in zip(
                self._system.components,
                self._system._current_variables,
            )
        }

    def post(self) -> "SymbolicConfigurationSet":
        """
        Return successors after exactly one transition.

        Asynchronous and general transitions are non-reflexive. Synchronous
        dynamics retain the transition from a fixed point to itself.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic()
        >>> system.configurations({"A": 0}).post().enumerate()
        ({'A': 1},)

        Returns
        -------
        SymbolicConfigurationSet
            Exact union of one-step successors.
        """

        return self._system.post(self)

    def pre(self) -> "SymbolicConfigurationSet":
        """
        Return predecessors after exactly one transition.

        Asynchronous and general transitions are non-reflexive. Synchronous
        dynamics retain the transition from a fixed point to itself.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic()
        >>> system.configurations({"A": 1}).pre().enumerate()
        ({'A': 0},)

        Returns
        -------
        SymbolicConfigurationSet
            Exact union of one-step predecessors.
        """

        return self._system.pre(self)

    def reachable(
        self,
        *,
        within: Optional["SymbolicConfigurationSet"] = None,
    ) -> "SymbolicConfigurationSet":
        """
        Return all configurations reachable from this set.

        This forward reachability closure includes this set and therefore uses
        paths of zero or more transitions. If `within` is supplied,
        configurations outside that region and transitions leaving it are
        excluded.

        Examples
        --------
        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": 1, "B": "A"}).symbolic()
        >>> initial = system.configurations({"A": 0, "B": 0})
        >>> initial.reachable().enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 0}, {'A': 1, 'B': 1})

        Parameters
        ----------
        within: SymbolicConfigurationSet, optional
            Region in which paths must remain. If `None`, use the complete
            configuration space.

        Returns
        -------
        SymbolicConfigurationSet
            Configurations reachable from this set through zero or more
            transitions.

        Raises
        ------
        TypeError
            If `within` is not a `SymbolicConfigurationSet`.
        ValueError
            If `within` belongs to another transition system.
        """

        return self._system.reachable(self, within=within)

    def coreachable(
        self,
        *,
        within: Optional["SymbolicConfigurationSet"] = None,
    ) -> "SymbolicConfigurationSet":
        """
        Return all configurations from which this set is reachable.

        This backward reachability closure includes this set and therefore
        uses paths of zero or more transitions. If `within` is supplied,
        configurations outside that region and transitions leaving it are
        excluded.

        Examples
        --------
        Both configurations of this switch can reach `A = 1`:

        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": "~A"}).symbolic()
        >>> target = system.configurations({"A": 1})
        >>> target.coreachable().enumerate()
        ({'A': 0}, {'A': 1})

        Parameters
        ----------
        within: SymbolicConfigurationSet, optional
            Region in which paths must remain. If `None`, use the complete
            configuration space.

        Returns
        -------
        SymbolicConfigurationSet
            Configurations from which this set is reachable.

        Raises
        ------
        TypeError
            If `within` is not a `SymbolicConfigurationSet`.
        ValueError
            If `within` belongs to another transition system.
        """

        return self._system.coreachable(self, within=within)

    def terminal_sccs(self) -> Tuple["SymbolicConfigurationSet", ...]:
        """
        Return strongly connected components with no outgoing transition in this set.

        Terminality is evaluated in the subgraph induced by this symbolic set.
        The returned components therefore need not be terminal in the complete
        transition graph.

        If this region is transition-closed, these terminal components
        correspond to the attractors contained in it.

        Examples
        --------
        The global terminal SCC is `11`. In the subgraph induced by `A = 0`,
        `00` is terminal instead:

        >>> from bonesistools.logic.bn import BooleanNetwork
        >>> system = BooleanNetwork({"A": 1, "B": "A"}).symbolic()
        >>> system.universe().terminal_sccs()[0].enumerate()
        ({'A': 1, 'B': 1},)
        >>> region = system.configurations({"A": 0})
        >>> region.terminal_sccs()[0].enumerate()
        ({'A': 0, 'B': 0},)

        Notes
        -----
        Small regions are handled explicitly. Larger transition-closed
        asynchronous regions use interleaved transition-guided reduction [1]
        followed by the Xie-Beerel terminal-SCC procedure [2].

        Returns
        -------
        tuple of SymbolicConfigurationSet
            Exact terminal strongly connected components.

        References
        ----------
        [1] Benes et al. (2021). Computing bottom SCCs symbolically using
        transition guided reduction. International Conference on Computer
        Aided Verification, 505-528.

        [2] Xie and Beerel (2000). Implicit enumeration of strongly connected
        components and an application to formal verification. IEEE
        Transactions on Computer-Aided Design of Integrated Circuits and
        Systems, 19(10), 1225-1230.
        """

        return self._system.terminal_sccs(self)

    def _compatible_node(
        self,
        other: "SymbolicConfigurationSet",
    ) -> _BDDConfigurationSetNode:
        """Return another set node after checking manager identity."""

        if other._system is not self._system:
            raise ValueError(
                "symbolic configuration sets belong to different transition systems"
            )
        return other._node

    @classmethod
    def _from_node(
        cls,
        system: SymbolicTransitionSystem,
        node: _BDDConfigurationSetNode,
    ) -> "SymbolicConfigurationSet":
        """Wrap one BDD node without validation, copying or counting."""

        configurations = cls.__new__(cls)
        configurations._system = system
        configurations._node = node
        return configurations
