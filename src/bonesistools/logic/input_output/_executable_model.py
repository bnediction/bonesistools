#!/usr/bin/env python

"""
Executable logical-model container.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from ..boolean_algebra import Hypercube
from ..boolean_network import BooleanNetwork
from ..influence_graph import InfluenceGraph


@dataclass
class ExecutableModel:
    """
    Container for executable logical-model data.

    `ExecutableModel` groups a Boolean network, an influence graph, initial
    states, perturbations and metadata without forcing all source formats into
    a single Boolean-network representation.
    """

    boolean_network: Optional[BooleanNetwork] = None
    influence_graph: Optional[InfluenceGraph] = None
    initial_states: Dict[str, Hypercube] = field(default_factory=dict)
    perturbations: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        """
        Return a compact representation of the executable model.

        Large containers are summarized by their sizes to keep notebook
        displays readable while preserving full access through attributes.
        """

        boolean_network = self.boolean_network
        if boolean_network is None:
            boolean_network_repr = "boolean_network=None"
        else:
            boolean_network_repr = (
                f"boolean_network(components={len(boolean_network.components)})"
            )

        influence_graph = self.influence_graph
        if influence_graph is None:
            influence_graph_repr = "influence_graph=None"
        else:
            influence_graph_repr = (
                f"influence_graph(nodes={influence_graph.number_of_nodes()}, "
                f"edges={influence_graph.number_of_edges()})"
            )

        return (
            f"{type(self).__name__}("
            f"{boolean_network_repr}, "
            f"{influence_graph_repr}, "
            f"initial_states(hypercubes={len(self.initial_states)}), "
            f"perturbations(n={len(self.perturbations)}), "
            f"metadata(entries={len(self.metadata)}))"
        )

    def get(self, attribute: str) -> Any:
        """
        Return an executable-model attribute or raise if it is unavailable.

        Supported attributes are `boolean_network`, `influence_graph`,
        `initial_states`, `perturbations` and `metadata`.

        Raises
        ------
        ValueError
            If `attribute` is unknown, or if the requested optional model
            component is unavailable.
        """

        if attribute == "boolean_network":
            if self.boolean_network is not None:
                return self.boolean_network

            reason = self.metadata.get(
                "boolean_network_reason",
                "the source model could not be converted to a Boolean network",
            )

            raise ValueError(f"Boolean network unavailable: {reason}")

        if attribute == "influence_graph":
            if self.influence_graph is not None:
                return self.influence_graph

            reason = self.metadata.get(
                "influence_graph_reason",
                "the source model did not contain signed interactions",
            )

            raise ValueError(f"Influence graph unavailable: {reason}")

        if attribute == "initial_states":
            return self.initial_states

        if attribute == "perturbations":
            return self.perturbations

        if attribute == "metadata":
            return self.metadata

        raise ValueError(
            "unsupported executable model attribute: expected "
            "'boolean_network', 'influence_graph', 'initial_states', "
            "'perturbations' or 'metadata' "
            f"but received {attribute!r}"
        )
