#!/usr/bin/env python

from __future__ import annotations

from collections.abc import (
    Iterable,
    Iterator,
    Mapping,
)
from typing import (
    Any,
    Dict,
    FrozenSet,
    MutableMapping,
    MutableSet,
    Optional,
    Set,
    Tuple,
    cast,
    overload,
)

from ._boolean import PartialBoolean
from ._typing import (
    HypercubeLike,
    PartialBooleanLike,
    is_hypercube_like,
)


class Hypercube(MutableMapping[str, PartialBoolean]):
    """
    Partial Boolean hypercube.

    A Hypercube represents a partial Boolean configuration mapping components
    to PartialBoolean values. Missing components are implicitly interpreted as
    free (`"*"`), allowing comparison between hypercubes defined on different
    component subsets.

    Hypercube ordering corresponds to inclusion of represented configurations:
    a hypercube is smaller than another one when it is more specific, i.e. when
    every configuration it represents is also represented by the other
    hypercube. This is the product order induced by the set-theoretic
    PartialBoolean order, not a Kleene truth order.

    Hypercube behaves as a mutable mapping from component names to
    PartialBoolean values. Accepted values are converted to PartialBoolean
    instances when assigned.

    Examples
    --------
    Missing components are interpreted as free dimensions:

    >>> Hypercube({"A": 0}) == Hypercube({"A": 0, "B": "*"})
    True

    More specific hypercubes are smaller than more general hypercubes:

    >>> Hypercube({"A": 0, "B": 1}) < Hypercube({"A": 0})
    True

    Parameters
    ----------
    mapping: Mapping[str, PartialBooleanLike] (optional, default: None)
        Hypercube component mapping. Values may be fixed Boolean values
        (`0`, `1`, `False`, `True`) or free values (`"*"` or compatible
        PartialBoolean-like values).

    Raises
    ------
    ValueError
        If an unsupported PartialBoolean value is provided.
    """

    __slots__ = ("_values",)

    def __init__(
        self,
        mapping: Optional[Mapping[str, PartialBooleanLike]] = None,
    ) -> None:

        self._values: Dict[str, PartialBoolean] = {}

        if mapping is not None:
            for component, value in mapping.items():
                self[component] = value

    def __getitem__(self, component: str) -> PartialBoolean:
        """
        Return the value explicitly associated with a component.

        Parameters
        ----------
        component: str
            Component name.

        Returns
        -------
        PartialBoolean
            Stored PartialBoolean value.

        Raises
        ------
        KeyError
            If `component` is not explicitly specified.
        """

        return self._values[component]

    def __setitem__(
        self,
        component: str,
        value: PartialBooleanLike,
    ) -> None:
        """
        Set a component value.

        Parameters
        ----------
        component: str
            Component name.
        value: PartialBooleanLike
            Value to assign. Supported values are converted to PartialBoolean.

        Raises
        ------
        ValueError
            If `value` cannot be converted to a PartialBoolean.
        """

        self._values[component] = (
            value if isinstance(value, PartialBoolean) else PartialBoolean(value)
        )

    def __delitem__(self, component: str) -> None:
        """
        Remove an explicitly specified component.

        Parameters
        ----------
        component: str
            Component name.

        Raises
        ------
        KeyError
            If `component` is not explicitly specified.
        """

        del self._values[component]

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over explicitly specified component names.

        Returns
        -------
        Iterator[str]
            Iterator over component names.
        """

        return iter(self._values)

    def __len__(self) -> int:
        """
        Return the number of explicitly specified components.

        Returns
        -------
        int
            Number of explicitly specified components.
        """

        return len(self._values)

    @overload
    def update(self, mapping: Mapping[str, PartialBooleanLike]) -> None: ...

    @overload
    def update(self, mapping: Iterable[Tuple[str, PartialBooleanLike]]) -> None: ...

    @overload
    def update(self, **kwargs: PartialBooleanLike) -> None: ...

    def update(  # pyright: ignore[reportIncompatibleMethodOverride]
        self,
        *args: Any,
        **kwargs: PartialBooleanLike,
    ) -> None:
        """
        Update explicitly specified component values.

        Values are converted to `PartialBoolean`, following `__setitem__`.
        """

        if len(args) > 1:
            raise TypeError(
                f"update expected at most 1 positional argument, got {len(args)}"
            )

        if args:
            other = args[0]
            if isinstance(other, Mapping):
                iterable = other.items()
            else:
                iterable = other

            for component, value in iterable:
                self[component] = value

        for component, value in kwargs.items():
            self[component] = value

    def rename(self, old: str, new: str) -> None:
        """
        Rename one explicitly specified hypercube component.

        Examples
        --------
        >>> hc = Hypercube({"Trp53": 1})
        >>> hc.rename("Trp53", "TP53")
        >>> hc
        Hypercube(TP53=1)

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
            If `old` is not explicitly specified in the hypercube.
        ValueError
            If `new` already exists.
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
            raise KeyError(f"component {old!r} not found")

        self.relabel({old: new})

    def relabel(self, mapping: Mapping[str, str]) -> None:
        """
        Relabel several explicitly specified hypercube components.

        Components absent from the hypercube are ignored. Relabeling must not
        merge two explicitly specified components, because that would silently
        drop or contradict one hypercube constraint.

        Examples
        --------
        >>> hc = Hypercube({"Trp53": 1, "Myc": 0})
        >>> hc.relabel({"Trp53": "TP53", "Myc": "MYC"})
        >>> hc
        Hypercube(TP53=1, MYC=0)

        Parameters
        ----------
        mapping: Mapping[str, str]
            Component rename mapping.

        Raises
        ------
        TypeError
            If `mapping` is not a mapping from strings to strings.
        ValueError
            If relabeling would merge two explicitly specified components.
        """

        if not isinstance(mapping, Mapping):
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
            active_mapping.get(component, component) for component in self._values
        ]
        if len(set(renamed_components)) != len(renamed_components):
            raise ValueError("relabeling would merge hypercube components")

        self._values = {
            active_mapping.get(component, component): value
            for component, value in self._values.items()
        }

    def __repr__(self) -> str:
        """
        Return a representation of the hypercube.

        Returns
        -------
        str
            Representation containing explicitly specified component values.
        """

        values = ", ".join(
            f"{_format_hypercube_component(component)}={value}"
            for component, value in self._values.items()
        )

        return f"Hypercube({values})"

    def __eq__(self, other: object) -> bool:
        """
        Test hypercube equality.

        Missing components are interpreted as free values. Therefore,
        `Hypercube({"A": 0})` and `Hypercube({"A": 0, "B": "*"})` are
        considered equal.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if both hypercubes represent the same configurations.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        for component in self.components | other.components:
            if self._get_value(component) != other._get_value(component):
                return False

        return True

    def __le__(self, other: object) -> bool:
        """
        Test whether the hypercube is included in another hypercube.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if every configuration represented by the current hypercube is
            also represented by `other`.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        for component in self.components | other.components:
            value1 = self._get_value(component)
            value2 = other._get_value(component)

            if value1 not in value2:
                return False

        return True

    def __lt__(self, other: object) -> bool:
        """
        Test whether the hypercube is strictly included in another hypercube.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube is included in `other` and both
            hypercubes are not equal.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        return self <= other and self != other

    def __ge__(self, other: object) -> bool:
        """
        Test whether the hypercube contains another hypercube.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if every configuration represented by `other` is also
            represented by the current hypercube.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        return other <= self

    def __gt__(self, other: object) -> bool:
        """
        Test whether the hypercube strictly contains another hypercube.

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube contains `other` and both hypercubes
            are not equal.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        return other < self

    def __hash__(self) -> int:
        """
        Return a hash consistent with hypercube equality.

        Explicit free values are ignored because they are equivalent to missing
        components.

        Returns
        -------
        int
            Hash value based on fixed component assignments.
        """

        normalized = tuple(
            sorted(
                (component, value)
                for component, value in self._values.items()
                if value.is_fixed
            )
        )

        return hash(normalized)

    def copy(self) -> "Hypercube":
        """
        Return a shallow copy of the hypercube.

        Examples
        --------
        >>> hc = Hypercube({"A": 0})
        >>> copied = hc.copy()
        >>> copied
        Hypercube(A=0)

        Returns
        -------
        Hypercube
            New Hypercube containing the same component values.
        """

        return Hypercube(self._values)

    def drop(
        self,
        components: Iterable[str],
        inplace: bool = False,
    ) -> Optional["Hypercube"]:
        """
        Remove components from the hypercube.

        Examples
        --------
        >>> hc = Hypercube({"A": 0, "B": 1})
        >>> hc.drop(["B"])
        Hypercube(A=0)
        >>> hc
        Hypercube(A=0, B=1)

        >>> hc.drop(["B"], inplace=True) is None
        True
        >>> hc
        Hypercube(A=0)

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

    @property
    def components(self) -> FrozenSet[str]:
        """
        Return explicitly specified hypercube components.

        Missing components are not returned, even though they are interpreted
        as free values during comparisons.

        Examples
        --------
        >>> sorted(Hypercube({"A": 0, "B": "*"}).components)
        ['A', 'B']

        Returns
        -------
        frozenset
            Explicitly specified component names.
        """

        return frozenset(self._values)

    @property
    def is_fully_specified(self) -> bool:
        """
        Test whether explicitly specified components are all fixed.

        This property only considers components present in the hypercube.
        Missing components are ignored.

        Examples
        --------
        >>> Hypercube({"A": 0, "B": 1}).is_fully_specified
        True
        >>> Hypercube({"A": 0, "B": "*"}).is_fully_specified
        False

        Returns
        -------
        bool
            True if the hypercube contains no explicit free value.
        """

        return all(value.is_fixed for value in self._values.values())

    def contains(self, other: object) -> bool:
        """
        Test whether the hypercube contains another hypercube.

        A hypercube contains another one when every configuration represented
        by `other` is also represented by the current hypercube. Equivalently,
        the current hypercube is greater than or equal to `other` in the
        inclusion order.

        Examples
        --------
        >>> Hypercube({"A": 0}).contains({"A": 0, "B": 1})
        True
        >>> Hypercube({"A": 0, "B": 1}).contains({"A": 0})
        False

        Parameters
        ----------
        other: object
            Hypercube-like object to test.

        Returns
        -------
        bool
            Whether every configuration represented by `other` is also
            represented by the current hypercube.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        other = self._coerce_hypercube(other)

        return other <= self

    def identical(self, other: object) -> Set[str]:
        """
        Return components sharing identical values.

        Missing components are interpreted as free values.

        Examples
        --------
        >>> hc1 = Hypercube({"A": 0, "B": "*"})
        >>> hc2 = Hypercube({"A": 0, "C": "*"})
        >>> sorted(hc1.identical(hc2))
        ['A', 'B', 'C']

        Parameters
        ----------
        other: object
            Hypercube-like object to compare.

        Returns
        -------
        Set[str]
            Components with identical values.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
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

        Examples
        --------
        >>> hc1 = Hypercube({"A": 0, "B": 1})
        >>> hc2 = Hypercube({"A": 0, "B": "*", "C": 1})
        >>> sorted(hc1.different(hc2))
        ['B', 'C']

        Parameters
        ----------
        other: object
            Hypercube-like object to compare.

        Returns
        -------
        Set[str]
            Components with different values.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
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

        This is an explicit alias for the `<=` operator.

        Examples
        --------
        A fixed value for `B` is more specific than a free or missing `B`, so
        it is smaller:

        >>> Hypercube({"A": 0, "B": 1}).is_smaller_than({"A": 0})
        True

        The reverse direction is false because the more general hypercube is
        not included in the more specific one:

        >>> Hypercube({"A": 0}).is_smaller_than({"A": 0, "B": 1})
        False

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube is included in `other`.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        return self <= other

    def is_larger_than(self, other: object) -> bool:
        """
        Test whether the hypercube contains another hypercube.

        This is an explicit alias for the `>=` operator.

        Examples
        --------
        A missing component is interpreted as free, so this hypercube contains
        the more specific one:

        >>> Hypercube({"A": 0}).is_larger_than({"A": 0, "B": 1})
        True

        The reverse direction is false:

        >>> Hypercube({"A": 0, "B": 1}).is_larger_than({"A": 0})
        False

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube contains `other`.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        return self >= other

    def is_strictly_smaller_than(self, other: object) -> bool:
        """
        Test whether the hypercube is strictly included in another hypercube.

        This is an explicit alias for the `<` operator.

        Examples
        --------
        >>> Hypercube({"A": 0, "B": 1}).is_strictly_smaller_than({"A": 0})
        True

        Equal hypercubes are not strictly smaller:

        >>> Hypercube({"A": 0}).is_strictly_smaller_than({"A": 0, "B": "*"})
        False

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube is strictly included in `other`.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        return self < other

    def is_strictly_larger_than(self, other: object) -> bool:
        """
        Test whether the hypercube strictly contains another hypercube.

        This is an explicit alias for the `>` operator.

        Examples
        --------
        >>> Hypercube({"A": 0}).is_strictly_larger_than({"A": 0, "B": 1})
        True

        Equal hypercubes are not strictly larger:

        >>> Hypercube({"A": 0, "B": "*"}).is_strictly_larger_than({"A": 0})
        False

        Parameters
        ----------
        other: object
            Hypercube-like object to compare against.

        Returns
        -------
        bool
            True if the current hypercube strictly contains `other`.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        return self > other

    def _get_value(self, component: str) -> PartialBoolean:
        """
        Return component value.

        Missing components are interpreted as free values.

        Parameters
        ----------
        component: str
            Component name.

        Returns
        -------
        PartialBoolean
            Explicitly stored value if present, otherwise a free
            PartialBoolean value.
        """

        return self._values.get(component, PartialBoolean("*"))

    @staticmethod
    def _coerce_hypercube(other: object) -> "Hypercube":
        """
        Convert supported objects into Hypercube instances.

        Parameters
        ----------
        other: object
            Object to coerce.

        Returns
        -------
        Hypercube
            Coerced hypercube. Existing Hypercube instances are returned as-is.

        Raises
        ------
        TypeError
            If `other` cannot be interpreted as a hypercube.
        ValueError
            If `other` contains unsupported PartialBoolean values.
        """

        if isinstance(other, Hypercube):
            return other

        if isinstance(other, dict):
            return Hypercube(other)

        if is_hypercube_like(other):
            return Hypercube(other)

        raise TypeError(
            "unsupported argument type for 'other': "
            "expected hypercube-like object "
            f"but received {type(other)}"
        )


class HypercubeCollection(MutableSet[Hypercube]):
    """
    Mutable set of hypercubes.

    A HypercubeCollection stores unique hypercubes and provides operations
    based on the hypercube inclusion order. Hypercubes are interpreted as
    partial Boolean configurations, where missing components and free values
    (`"*"` ) correspond to unconstrained dimensions.

    The collection behaves as a mutable set: duplicate hypercubes are removed,
    iteration order is not part of the public semantics, and membership is based
    on hypercube equality.

    Examples
    --------
    Duplicate and equivalent hypercubes are stored once:

    >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": "*"}])
    >>> len(hcs)
    1

    Collection filtering follows the same inclusion order as Hypercube:

    >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": 1}])
    >>> list(hcs.smaller_than({"A": 0, "B": 1}))
    [Hypercube(A=0, B=1)]

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

        self._hypercubes: Set[Hypercube] = set()

        if hypercubes is not None:
            for hypercube in hypercubes:
                self.add(hypercube)

    def __contains__(self, hypercube: object) -> bool:
        """
        Test whether a hypercube belongs to the collection.

        Invalid or unsupported objects are treated as absent rather than
        raising an exception.

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}])
        >>> {"A": 0, "B": "*"} in hcs
        True
        >>> {"A": 1} in hcs
        False

        Parameters
        ----------
        hypercube: object
            Hypercube-like object to test.

        Returns
        -------
        bool
            True if `hypercube` belongs to the collection.
        """

        try:
            hypercube = self._coerce_hypercube(hypercube)

        except (TypeError, ValueError):
            return False

        return hypercube in self._hypercubes

    def __iter__(self) -> Iterator[Hypercube]:
        """
        Iterate over stored hypercubes.

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}])
        >>> list(hcs)
        [Hypercube(A=0)]

        Returns
        -------
        Iterator[Hypercube]
            Iterator over stored hypercubes.
        """

        return iter(self._hypercubes)

    def __len__(self) -> int:
        """
        Return the number of hypercubes in the collection.

        Examples
        --------
        >>> len(HypercubeCollection([{"A": 0}, {"A": 0, "B": "*"}]))
        1

        Returns
        -------
        int
            Number of unique hypercubes in the collection.
        """

        return len(self._hypercubes)

    def __eq__(self, other: object) -> bool:
        """
        Test whether two hypercube collections contain identical hypercubes.

        `other` may be another HypercubeCollection or an iterable of
        hypercube-like objects.

        Examples
        --------
        >>> HypercubeCollection([{"A": 0}]) == HypercubeCollection(
        ...     [{"A": 0, "B": "*"}]
        ... )
        True

        Parameters
        ----------
        other: object
            Collection-like object to compare against.

        Returns
        -------
        bool or NotImplemented
            True if both collections contain identical hypercubes. Returns
            NotImplemented when `other` cannot be interpreted as a collection.
        """
        if not isinstance(other, HypercubeCollection):
            if not isinstance(other, Iterable):
                return NotImplemented

            try:
                other = HypercubeCollection(cast(Iterable[HypercubeLike], other))
            except (TypeError, ValueError):
                return NotImplemented

        return self._hypercubes == other._hypercubes

    def copy(self) -> "HypercubeCollection":
        """
        Return a shallow copy of the collection.

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}])
        >>> copied = hcs.copy()
        >>> list(copied)
        [Hypercube(A=0)]

        Returns
        -------
        HypercubeCollection
            New collection containing the same hypercubes.
        """

        return HypercubeCollection(self._hypercubes)

    def add(self, value: HypercubeLike) -> None:
        """
        Add a hypercube to the collection.

        Examples
        --------
        >>> hcs = HypercubeCollection()
        >>> hcs.add({"A": 0})
        >>> {"A": 0} in hcs
        True

        Parameters
        ----------
        value: HypercubeLike
            Hypercube-like object to add.

        Raises
        ------
        TypeError
            If `value` cannot be interpreted as a hypercube.
        ValueError
            If `value` contains an unsupported PartialBoolean value.
        """

        self._hypercubes.add(self._coerce_hypercube(value))

    def discard(self, value: HypercubeLike) -> None:
        """
        Remove a hypercube from the collection if present.

        Invalid or unsupported objects are ignored.

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}])
        >>> hcs.discard({"A": 0, "B": "*"})
        >>> len(hcs)
        0

        Parameters
        ----------
        value: HypercubeLike
            Hypercube-like object to remove.
        """

        try:
            hypercube = self._coerce_hypercube(value)

        except (TypeError, ValueError):
            return None

        self._hypercubes.discard(hypercube)

    @property
    def components(self) -> FrozenSet[str]:
        """
        Return components appearing in at least one hypercube.

        Examples
        --------
        >>> sorted(HypercubeCollection([{"A": 0}, {"B": 1}]).components)
        ['A', 'B']

        Returns
        -------
        frozenset
            Components appearing in at least one stored hypercube.
        """

        components = set()

        for hypercube in self:
            components.update(hypercube.components)

        return frozenset(components)

    def smaller_than(self, other: HypercubeLike) -> "HypercubeCollection":
        """
        Return hypercubes included in another hypercube.

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": 1}, {"A": 1}])
        >>> list(hcs.smaller_than({"A": 0, "B": 1}))
        [Hypercube(A=0, B=1)]

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

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}, {"A": 1}])
        >>> list(hcs.larger_than({"A": 0, "B": 1}))
        [Hypercube(A=0)]

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

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": 1}])
        >>> list(hcs.maximal())
        [Hypercube(A=0)]

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

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": 1}])
        >>> list(hcs.minimal())
        [Hypercube(A=0, B=1)]

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

        Examples
        --------
        >>> hcs = HypercubeCollection([{"A": 0}, {"A": 0, "B": 1}])
        >>> list(hcs.fully_specified())
        [Hypercube(A=0, B=1)]

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
            "unsupported argument type for 'hypercube': "
            "expected hypercube-like object "
            f"but received {type(hypercube)}"
        )


def _format_hypercube_component(component: str) -> str:
    """
    Return a compact component name for Hypercube representations.
    """

    if "," in component:
        return repr(component)

    return component
