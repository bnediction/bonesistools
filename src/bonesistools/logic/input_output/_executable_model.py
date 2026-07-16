#!/usr/bin/env python

"""
Executable logical-model container.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Tuple,
    Union,
    cast,
)

from ..._compat import Literal
from ..boolean_algebra import ConfigurationSet, Hypercube
from ..boolean_algebra._typing import Configuration, HypercubeLike
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph


class ExecutableModel:
    """
    Protected snapshot of an executable logical model.

    An executable model groups a Boolean network with its influence graph,
    named initial conditions, parameters, perturbations and source metadata.
    All inputs are copied at construction time, and accessors return independent
    copies. Consequently, modifying a returned object cannot invalidate the
    model's internal consistency.

    Parameters
    ----------
    boolean_network: BooleanNetwork, optional
        Executable Boolean network. If supplied without `influence_graph`, its
        influence graph is derived automatically.
    influence_graph: InfluenceGraph, optional
        Signed influence graph. If supplied with `boolean_network`, it must
        contain the same components. Its edges may preserve conservative
        influences from the source format in addition to those induced by
        simplified rules.
    initial_conditions: Mapping[str, HypercubeLike], optional
        Named, possibly partial, Boolean initial configurations.
    parameters: Mapping[str, Any], optional
        Executable or simulation parameters associated with the model.
    perturbations: Mapping[str, Any], optional
        Named perturbation definitions.
    metadata: Mapping[str, Any], optional
        Source-format metadata and information not represented by the other
        model objects.

    Raises
    ------
    TypeError
        If a model object or named mapping has an unsupported type.
    ValueError
        If `boolean_network` and `influence_graph` are inconsistent, or an
        initial condition refers to a component outside the Boolean network.
    """

    __slots__ = (
        "__boolean_network",
        "__influence_graph",
        "__initial_conditions",
        "__parameters",
        "__perturbations",
        "__metadata",
    )

    def __init__(
        self,
        *,
        boolean_network: Optional[BooleanNetwork] = None,
        influence_graph: Optional[InfluenceGraph] = None,
        initial_conditions: Optional[Mapping[str, HypercubeLike]] = None,
        parameters: Optional[Mapping[str, Any]] = None,
        perturbations: Optional[Mapping[str, Any]] = None,
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> None:

        if boolean_network is not None and not isinstance(
            boolean_network,
            BooleanNetwork,
        ):
            raise TypeError(
                "unsupported argument type for 'boolean_network': "
                f"expected {BooleanNetwork} or None but received "
                f"{type(boolean_network)}"
            )

        if influence_graph is not None and not isinstance(
            influence_graph,
            InfluenceGraph,
        ):
            raise TypeError(
                "unsupported argument type for 'influence_graph': "
                f"expected {InfluenceGraph} or None but received "
                f"{type(influence_graph)}"
            )

        copied_network = None if boolean_network is None else boolean_network.copy()
        if influence_graph is None and copied_network is not None:
            copied_graph = copied_network.to_influence_graph()
        else:
            copied_graph = (
                None if influence_graph is None else deepcopy(influence_graph)
            )

        if copied_network is not None and copied_graph is not None:
            self._validate_influence_graph(copied_network, copied_graph)

        copied_initial_conditions = self._copy_initial_conditions(
            initial_conditions,
            network=copied_network,
        )

        self.__boolean_network = copied_network
        self.__influence_graph = copied_graph
        self.__initial_conditions = copied_initial_conditions
        self.__parameters = self._copy_mapping(parameters, name="parameters")
        self.__perturbations = self._copy_mapping(
            perturbations,
            name="perturbations",
        )
        self.__metadata = self._copy_mapping(metadata, name="metadata")

    def __repr__(self) -> str:
        """
        Return a compact representation of the executable model.

        Large containers are summarized by their sizes to keep notebook
        displays readable.
        """

        if self.__boolean_network is None:
            network_repr = "boolean_network=None"
        else:
            network_repr = (
                f"boolean_network(components={len(self.__boolean_network.components)})"
            )

        if self.__influence_graph is None:
            graph_repr = "influence_graph=None"
        else:
            graph_repr = (
                "influence_graph(nodes="
                f"{self.__influence_graph.number_of_nodes()}, "
                f"edges={self.__influence_graph.number_of_edges()})"
            )

        return (
            f"{type(self).__name__}("
            f"{network_repr}, "
            f"{graph_repr}, "
            f"initial_conditions(hypercubes={len(self.__initial_conditions)}), "
            f"parameters(n={len(self.__parameters)}), "
            f"perturbations(n={len(self.__perturbations)}), "
            f"metadata(entries={len(self.__metadata)}))"
        )

    @property
    def boolean_network(self) -> BooleanNetwork:
        """
        Return an independent copy of the Boolean network.

        Raises
        ------
        ValueError
            If the source model could not be converted to a Boolean network.
        """

        return self._require_boolean_network().copy()

    @property
    def influence_graph(self) -> InfluenceGraph:
        """
        Return an independent copy of the influence graph.

        Raises
        ------
        ValueError
            If no signed influence graph is available.
        """

        if self.__influence_graph is None:
            reason = self.__metadata.get(
                "influence_graph_reason",
                "the source model did not contain signed interactions",
            )
            raise ValueError(f"Influence graph unavailable: {reason}")

        return deepcopy(self.__influence_graph)

    def initial_conditions(self) -> Dict[str, Hypercube]:
        """Return independent copies of the named initial conditions."""

        return {
            name: condition.copy()
            for name, condition in self.__initial_conditions.items()
        }

    def parameters(self) -> Dict[str, Any]:
        """Return an independent copy of the executable parameters."""

        return deepcopy(self.__parameters)

    def perturbations(self) -> Dict[str, Any]:
        """Return an independent copy of the perturbation definitions."""

        return deepcopy(self.__perturbations)

    def metadata(self) -> Dict[str, Any]:
        """Return an independent copy of the source metadata."""

        return deepcopy(self.__metadata)

    def reachable_configurations(
        self,
        initial: Union[str, HypercubeLike],
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "general",
    ) -> Iterator[Configuration]:
        """
        Iterate over configurations reachable from an initial condition.

        `initial` may be a named condition stored by the model or a directly
        supplied complete Boolean configuration.
        """

        return self._require_boolean_network().reachable_configurations(
            self._resolve_condition(initial),
            update=update,
        )

    def attractors(
        self,
        initial: Optional[Union[str, HypercubeLike]] = None,
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "general",
    ) -> Tuple[ConfigurationSet, ...]:
        """
        Return attractors reachable from an optional initial condition.

        `initial` may be a named condition stored by the model or a directly
        supplied partial Boolean configuration. If `None`, the complete state
        space is considered initial.
        """

        resolved_initial = None if initial is None else self._resolve_condition(initial)
        return self._require_boolean_network().attractors(
            resolved_initial,
            update=update,
        )

    def transition(
        self,
        source: Union[str, HypercubeLike],
        target: Union[str, HypercubeLike],
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "general",
        quantifier: Literal["exists", "robust", "universal"] = "robust",
    ) -> bool:
        """Test a direct transition between two named or supplied conditions."""

        return self._require_boolean_network().transition(
            self._resolve_condition(source),
            self._resolve_condition(target),
            update=update,
            quantifier=quantifier,
        )

    def reachability(
        self,
        source: Union[str, HypercubeLike],
        target: Union[str, HypercubeLike],
        *,
        update: Literal[
            "asynchronous", "synchronous", "general", "most-permissive"
        ] = "general",
        quantifier: Literal["exists", "robust", "universal"] = "robust",
    ) -> bool:
        """Test reachability between two named or supplied conditions."""

        return self._require_boolean_network().reachability(
            self._resolve_condition(source),
            self._resolve_condition(target),
            update=update,
            quantifier=quantifier,
        )

    def save(self, file: Union[str, Path]) -> None:
        """
        Save the executable model as a ZGINML archive.

        Boolean GINML models retain their parsed logical structure and known
        GINsim companion data whenever possible. Models imported from a
        multi-valued source are exported as their Booleanized executable
        network; automatic reconstruction of the original multi-valued model
        is intentionally not attempted.

        Parameters
        ----------
        file: str or Path
            Output `.zginml` archive.

        Returns
        -------
        None
            Writes the model to `file`.

        Raises
        ------
        ValueError
            If `file` does not use the `.zginml` extension, or neither a
            Boolean network nor reusable GINML source metadata is available.
        """

        from ._zginml import _write_zginml

        file = Path(file)

        if file.suffix.lower() != ".zginml":
            raise ValueError(
                "unsupported executable-model output format: expected a '.zginml' file"
            )

        _write_zginml(
            file,
            boolean_network=self.__boolean_network,
            influence_graph=self.__influence_graph,
            initial_conditions=self.__initial_conditions,
            parameters=self.__parameters,
            perturbations=self.__perturbations,
            metadata=self.__metadata,
        )

    def _require_boolean_network(self) -> BooleanNetwork:
        """Return the protected Boolean network or explain its absence."""

        if self.__boolean_network is None:
            reason = self.__metadata.get(
                "boolean_network_reason",
                "the source model could not be converted to a Boolean network",
            )
            raise ValueError(f"Boolean network unavailable: {reason}")

        return self.__boolean_network

    def _resolve_condition(
        self,
        condition: Union[str, HypercubeLike],
    ) -> HypercubeLike:
        """Resolve a named initial condition without exposing internal data."""

        if not isinstance(condition, str):
            return condition

        try:
            return self.__initial_conditions[condition].copy()
        except KeyError as error:
            available = ", ".join(repr(name) for name in self.__initial_conditions)
            suffix = "" if not available else f"; available names are {available}"
            raise ValueError(
                f"unknown initial condition {condition!r}{suffix}"
            ) from error

    @staticmethod
    def _copy_initial_conditions(
        initial_conditions: Optional[Mapping[str, HypercubeLike]],
        *,
        network: Optional[BooleanNetwork],
    ) -> Dict[str, Hypercube]:
        """Validate and copy named initial conditions."""

        if initial_conditions is None:
            return {}

        if not isinstance(initial_conditions, Mapping):
            raise TypeError(
                "unsupported argument type for 'initial_conditions': "
                f"expected {Mapping} or None but received "
                f"{type(initial_conditions)}"
            )

        copied = {}
        for name, condition in initial_conditions.items():
            if not isinstance(name, str):
                raise TypeError(
                    "unsupported initial-condition name type: "
                    f"expected {str} but received {type(name)}"
                )

            hypercube = Hypercube(condition)

            if network is not None:
                unknown = set(hypercube) - network.components
                if unknown:
                    raise ValueError(
                        f"initial condition {name!r} contains unknown "
                        f"components: {', '.join(sorted(unknown))}"
                    )

            copied[name] = hypercube

        return copied

    @staticmethod
    def _copy_mapping(
        mapping: Optional[Mapping[str, Any]],
        *,
        name: str,
    ) -> Dict[str, Any]:
        """Validate and deeply copy a named mapping."""

        if mapping is None:
            return {}

        if not isinstance(mapping, Mapping):
            raise TypeError(
                f"unsupported argument type for {name!r}: "
                f"expected {Mapping} or None but received {type(mapping)}"
            )

        for key in mapping:
            if not isinstance(key, str):
                raise TypeError(
                    f"unsupported {name} key type: "
                    f"expected {str} but received {type(key)}"
                )

        return cast(Dict[str, Any], deepcopy(dict(mapping)))

    @staticmethod
    def _validate_influence_graph(
        network: BooleanNetwork,
        graph: InfluenceGraph,
    ) -> None:
        """Validate that graph and network use the same Boolean components."""

        if set(graph) != network.components:
            raise ValueError(
                "inconsistent executable model: the influence graph must "
                "contain exactly the Boolean-network components"
            )
