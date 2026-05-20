#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Iterable, Iterator, MutableMapping, MutableSet
from functools import total_ordering
from typing import Dict, Optional, Set, Union

from ._boolean import PartialBoolean
from ._typing import PartialBooleanLike, HypercubeLike, is_hypercube_like

class Hypercube(MutableMapping):
    """
    Partial Boolean hypercube.

    A Hypercube represents a partial Boolean configuration mapping components
    to PartialBoolean values. Missing components are implicitly interpreted as
    free (`"*"`), allowing comparison between hypercubes defined on different
    component subsets.

    Hypercube ordering corresponds to inclusion of represented configurations:
    a hypercube is smaller than another one when it is more specific, i.e. when
    every configuration it represents is also represented by the other
    hypercube.

    Examples
    --------
    The following hypercubes are considered equivalent:

    >>> Hypercube({"A": 0})
    >>> Hypercube({"A": 0, "B": "*"})
    
    Parameters
    ----------
    mapping: Mapping[str, PartialBooleanLike] (optional, default: None)
        Hypercube component mapping.

    Raises
    ------
    ValueError
        If an unsupported PartialBoolean value is provided.
    """

    __slots__ = ("_values",)

    def __init__(
        self,
        mapping: Optional[Dict[str, PartialBooleanLike]] = None,
    ) -> None:

        self._values: Dict[str, PartialBoolean] = {}

        if mapping is not None:
            for component, value in mapping.items():
                self[component] = value

    @property
    def components(self) -> frozenset:
        """
        Return explicitly specified hypercube components.
        """

        return frozenset(self._values)

    @property
    def is_fully_specified(self) -> bool:
        """
        Whether explicitly specified components are all fixed.

        Returns
        -------
        bool
            True if the hypercube contains no explicit free value.
        """

        return all(value.is_fixed for value in self._values.values())

    def copy(self) -> "Hypercube":
        """
        Return a shallow copy of the hypercube.
        """

        return Hypercube(self._values)

    def drop(
        self,
        components: Iterable[str],
        inplace: bool = False,
    ) -> Optional["Hypercube"]:
        """
        Remove components from the hypercube.

        Parameters
        ----------
        components: Iterable[str]
            Components to remove.
        inplace: bool (default: False)
            Whether to modify the current hypercube.

        Returns
        -------
        Hypercube or None
            Modified hypercube if `inplace=False`, otherwise None.
        """

        hypercube = self if inplace else self.copy()

        for component in components:
            hypercube.pop(component, None)

        if not inplace:
            return hypercube

        return None

    def contains(self, other: object) -> bool:
        """
        Test whether the hypercube contains another hypercube.

        Parameters
        ----------
        other: object
            Hypercube-like object to test.

        Returns
        -------
        bool
            Whether every configuration represented by `other` is also
            represented by the current hypercube.
        """

        other = self._coerce_hypercube(other)

        return other <= self

    def identical(self, other: object) -> Set[str]:
        """
        Return components sharing identical values.

        Missing components are interpreted as free values.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare.

        Returns
        -------
        Set[str]
            Components with identical values.
        """

        other = self._coerce_hypercube(other)

        identical = set()

        for component in self.components | other.components:
            if self._get_value(component) == other._get_value(component):
                identical.add(component)

        return identical

    def different(self, other: object) -> Set[str]:
        """
        Return components with different values.

        Missing components are interpreted as free values.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare.

        Returns
        -------
        Set[str]
            Components with different values.
        """

        other = self._coerce_hypercube(other)

        different = set()

        for component in self.components | other.components:
            if self._get_value(component) != other._get_value(component):
                different.add(component)

        return different

    def is_smaller_than(self, other: object) -> bool:
        """
        Test whether the hypercube is included in another hypercube.
        """

        return self <= other

    def is_larger_than(self, other: object) -> bool:
        """
        Test whether the hypercube contains another hypercube.
        """

        return self >= other

    def is_strictly_smaller_than(self, other: object) -> bool:
        """
        Test whether the hypercube is strictly included in another hypercube.
        """

        return self < other

    def is_strictly_larger_than(self, other: object) -> bool:
        """
        Test whether the hypercube strictly contains another hypercube.
        """

        return self > other

    def __getitem__(self, component: str) -> PartialBoolean:

        return self._values[component]

    def __setitem__(
        self,
        component: str,
        value: PartialBooleanLike,
    ) -> None:

        self._values[component] = (
            value if isinstance(value, PartialBoolean) else PartialBoolean(value)
        )

    def __delitem__(self, component: str) -> None:

        del self._values[component]

    def __iter__(self) -> Iterator[str]:

        return iter(self._values)

    def __len__(self) -> int:

        return len(self._values)

    def __repr__(self) -> str:

        return f"Hypercube({dict(self)!r})"

    def __eq__(self, other: object) -> bool:

        other = self._coerce_hypercube(other)

        for component in self.components | other.components:
            if self._get_value(component) != other._get_value(component):
                return False

        return True

    def __le__(self, other: object) -> bool:

        other = self._coerce_hypercube(other)

        for component in self.components | other.components:
            value1 = self._get_value(component)
            value2 = other._get_value(component)

            if not value2.contains(value1):
                return False

        return True

    def __lt__(self, other: object) -> bool:

        other = self._coerce_hypercube(other)

        return self <= other and self != other

    def __ge__(self, other: object) -> bool:

        other = self._coerce_hypercube(other)

        return other <= self

    def __gt__(self, other: object) -> bool:

        other = self._coerce_hypercube(other)

        return other < self

    def __hash__(self) -> int:

        normalized = tuple(
            sorted(
                (component, value)
                for component, value in self._values.items()
                if value.is_fixed
            )
        )

        return hash(normalized)

    def _get_value(self, component: str) -> PartialBoolean:
        """
        Return component value.

        Missing components are interpreted as free values.
        """

        return self._values.get(component, PartialBoolean("*"))

    @staticmethod
    def _coerce_hypercube(other: object) -> "Hypercube":
        """
        Convert supported objects into Hypercube instances.
        """

        if isinstance(other, Hypercube):
            return other

        if isinstance(other, dict):
            return Hypercube(other)

        if is_hypercube_like(other):
            return Hypercube(other)

        raise TypeError(
            "unsupported argument type: expected hypercube-like object, "
            f"but received {type(other)}"
        )

class HypercubeCollection(MutableSet):
    """
    Mutable set of hypercubes.

    A HypercubeCollection stores unique hypercubes and provides operations
    based on the hypercube inclusion order. Hypercubes are interpreted as
    partial Boolean configurations, where missing components and free values
    (`"*"` ) correspond to unconstrained dimensions.

    The collection behaves as a mutable set: duplicate hypercubes are removed,
    iteration order is not part of the public semantics, and membership is based
    on hypercube equality.

    Parameters
    ----------
    hypercubes: Iterable[HypercubeLike] (optional, default: None)
        Hypercubes used to initialise the collection. Non-Hypercube objects are
        coerced into Hypercube instances when possible.

    Raises
    ------
    TypeError
        If one of the provided objects cannot be interpreted as a hypercube.
    ValueError
        If one of the provided hypercubes contains an unsupported
        PartialBoolean value.
    """

    def __init__(
        self,
        hypercubes: Optional[Iterable[HypercubeLike]] = None,
    ) -> None:

        self._hypercubes = set()

        if hypercubes is not None:
            for hypercube in hypercubes:
                self.add(hypercube)

    def __contains__(self, hypercube: object) -> bool:
        """
        Test whether a hypercube belongs to the collection.
        """

        try:
            hypercube = self._coerce_hypercube(hypercube)

        except (TypeError, ValueError):
            return False

        return hypercube in self._hypercubes

    def __iter__(self) -> Iterator[Hypercube]:
        """
        Iterate over stored hypercubes.
        """

        return iter(self._hypercubes)

    def __len__(self) -> int:
        """
        Return the number of hypercubes in the collection.
        """

        return len(self._hypercubes)

    def __eq__(self, other: object) -> bool:
        """
        Test whether two hypercube collections contain identical hypercubes.
        """
        if not isinstance(other, HypercubeCollection):
            try:
                other = HypercubeCollection(other)
            except (TypeError, ValueError):
                return NotImplemented

        return self._hypercubes == other._hypercubes

    @property
    def components(self) -> frozenset:
        """
        Return components appearing in at least one hypercube.
        """

        components = set()

        for hypercube in self:
            components.update(hypercube.components)

        return frozenset(components)

    def add(self, hypercube: HypercubeLike) -> None:
        """
        Add a hypercube to the collection.

        Parameters
        ----------
        hypercube: HypercubeLike
            Hypercube-like object to add.

        Raises
        ------
        TypeError
            If `hypercube` cannot be interpreted as a hypercube.
        ValueError
            If `hypercube` contains an unsupported PartialBoolean value.
        """

        self._hypercubes.add(self._coerce_hypercube(hypercube))

    def discard(self, hypercube: HypercubeLike) -> None:
        """
        Remove a hypercube from the collection if present.

        Parameters
        ----------
        hypercube: HypercubeLike
            Hypercube-like object to remove.
        """

        try:
            hypercube = self._coerce_hypercube(hypercube)

        except (TypeError, ValueError):
            return None

        self._hypercubes.discard(hypercube)

    def copy(self) -> "HypercubeCollection":
        """
        Return a shallow copy of the collection.
        """

        return HypercubeCollection(self._hypercubes)

    def smaller_than(self, other: HypercubeLike) -> "HypercubeCollection":
        """
        Return hypercubes included in another hypercube.

        Parameters
        ----------
        other: HypercubeLike
            Reference hypercube.

        Returns
        -------
        HypercubeCollection
            Collection of hypercubes `h` such that `h <= other`.
        """

        other = self._coerce_hypercube(other)

        return HypercubeCollection(
            hypercube.copy() for hypercube in self if hypercube <= other
        )

    def larger_than(self, other: HypercubeLike) -> "HypercubeCollection":
        """
        Return hypercubes containing another hypercube.

        Parameters
        ----------
        other: HypercubeLike
            Reference hypercube.

        Returns
        -------
        HypercubeCollection
            Collection of hypercubes `h` such that `h >= other`.
        """

        other = self._coerce_hypercube(other)

        return HypercubeCollection(
            hypercube.copy() for hypercube in self if hypercube >= other
        )

    def maximal(self) -> "HypercubeCollection":
        """
        Return maximal hypercubes in the collection.

        A hypercube is maximal if it is not strictly included in any other
        hypercube from the collection. Maximal hypercubes correspond to the most
        general elements of the collection.

        Returns
        -------
        HypercubeCollection
            Collection of maximal hypercubes.
        """

        return HypercubeCollection(
            hypercube.copy()
            for hypercube in self
            if not any(hypercube < other for other in self)
        )

    def minimal(self) -> "HypercubeCollection":
        """
        Return minimal hypercubes in the collection.

        A hypercube is minimal if it does not strictly include any other
        hypercube from the collection. Minimal hypercubes correspond to the most
        specific elements of the collection.

        Returns
        -------
        HypercubeCollection
            Collection of minimal hypercubes.
        """

        return HypercubeCollection(
            hypercube.copy()
            for hypercube in self
            if not any(other < hypercube for other in self)
        )

    def fully_specified(self) -> "HypercubeCollection":
        """
        Return hypercubes with all collection components fixed.
    
        Missing components are interpreted as free values. A hypercube is therefore
        considered fully specified only if every component appearing in the
        collection is fixed in that hypercube.
    
        Returns
        -------
        HypercubeCollection
            Collection containing hypercubes fully specified with respect to the
            collection component set.
        """

        components = self.components

        return HypercubeCollection(
            hypercube.copy()
            for hypercube in self
            if all(hypercube._get_value(component).is_fixed for component in components)
        )

    @staticmethod
    def _coerce_hypercube(hypercube: object) -> Hypercube:
        """
        Convert supported objects into Hypercube instances.

        Parameters
        ----------
        hypercube: object
            Object to coerce.

        Returns
        -------
        Hypercube
            Coerced hypercube.

        Raises
        ------
        TypeError
            If `hypercube` is not hypercube-like.
        ValueError
            If `hypercube` contains unsupported PartialBoolean values.
        """

        if isinstance(hypercube, Hypercube):
            return hypercube.copy()

        if isinstance(hypercube, dict):
            return Hypercube(hypercube)

        if is_hypercube_like(hypercube):
            return Hypercube(hypercube)

        raise TypeError(
            "unsupported argument type: expected hypercube-like object, "
            f"but received {type(hypercube)}"
        )
