#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from itertools import product
from typing import List, Optional, Sequence, Tuple, Union, overload

from ..._typing import RandomStateSeed
from ..._validation import _as_positive_integer, _as_seed
from ._hypercube import Hypercube
from ._typing import Configuration, ConfigurationLike, HypercubeLike

_EncodedHypercube = Tuple[int, int]


class ConfigurationSet:
    """
    Exact finite set of Boolean configurations.

    A ConfigurationSet represents a set of complete Boolean configurations over
    a fixed ordered list of components. The internal representation is compact
    and intentionally hidden from the public API: users interact with concrete
    configurations through membership, iteration, counting and sampling.

    The stored representation is exact and non-redundant: adding a more general
    subspace removes more specific stored subspaces that it covers. The
    representation is not guaranteed to use the globally minimal number of
    internal subspaces. Equality therefore compares represented configurations,
    not the internal subspace encoding.

    Examples
    --------
    Create an empty set over two Boolean components:

    >>> states = ConfigurationSet(["A", "B"])
    >>> len(states)
    0

    Add a complete configuration:

    >>> states.add({"A": 0, "B": 0})
    >>> {"A": 0, "B": 0} in states
    True

    Add a partial configuration. Missing or free components are interpreted as
    all admissible Boolean values:

    >>> states.add({"A": "*", "B": 0})
    >>> len(states)
    2
    >>> list(states)
    [{'A': 0, 'B': 0}, {'A': 1, 'B': 0}]

    Iterate lazily over complete configurations, or materialize them with
    `enumerate()`:

    >>> states.enumerate()
    ({'A': 0, 'B': 0}, {'A': 1, 'B': 0})

    Sample configurations when the represented set is too large to enumerate:

    >>> states.sample(seed=0)
    {'A': 0, 'B': 0}

    Parameters
    ----------
    components: Iterable[str]
        Ordered Boolean components defining complete configurations.
    configurations: Iterable[HypercubeLike], optional
        Initial configurations or subspaces to add.

    Raises
    ------
    TypeError
        If components are not strings or configurations cannot be interpreted
        as Boolean subspaces.
    ValueError
        If components are duplicated, a value is unsupported or a configuration
        mentions an unknown component.
    """

    __slots__ = ("_component_to_index", "_components", "_hypercubes")

    def __init__(
        self,
        components: Iterable[str],
        configurations: Optional[Iterable[HypercubeLike]] = None,
    ) -> None:

        self._components = _coerce_components(components)
        self._component_to_index = {
            component: index for index, component in enumerate(self._components)
        }
        self._hypercubes: List[_EncodedHypercube] = []

        if configurations is not None:
            for configuration in configurations:
                self.add(configuration)

    def __contains__(self, configuration: object) -> bool:
        """
        Test whether a configuration or subspace is fully represented.

        Invalid or unsupported objects are treated as absent.
        """

        try:
            hypercube = self._coerce_hypercube(configuration)
        except (TypeError, ValueError):
            return False

        return self._contains_hypercube(hypercube)

    def __eq__(self, other: object) -> bool:
        """
        Test semantic equality of represented configurations.

        Two ConfigurationSet objects are equal when they represent exactly the
        same complete Boolean configurations over the same component names,
        regardless of the internal subspaces used to encode them.
        """

        if not isinstance(other, ConfigurationSet):
            return NotImplemented

        if set(self._components) != set(other._components):
            return False

        if self.count() != other.count():
            return False

        return self._covers(other) and other._covers(self)

    def __iter__(self) -> Iterator[Configuration]:
        """
        Iterate over complete configurations.

        Yields
        ------
        Configuration
            Complete Boolean configurations ordered by `components`.
        """

        for hypercube in self._hypercubes:
            yield from _iter_hypercube_configurations(hypercube, self._components)

    def __len__(self) -> int:
        """
        Return the exact number of represented configurations.
        """

        return self.count()

    def __repr__(self) -> str:
        """
        Return a compact representation.
        """

        return (
            f"ConfigurationSet(components={self._components!r}, "
            f"n_configurations={self.count()})"
        )

    def copy(self) -> "ConfigurationSet":
        """
        Return a shallow copy of the configuration set.

        Examples
        --------
        >>> states = ConfigurationSet(["A"], [{"A": "*"}])
        >>> copied = states.copy()
        >>> copied is states
        False
        >>> copied.enumerate()
        ({'A': 0}, {'A': 1})
        """

        copied = ConfigurationSet(self._components)
        copied._hypercubes = list(self._hypercubes)

        return copied

    def add(
        self,
        configuration: Union[ConfigurationLike, HypercubeLike],
    ) -> None:
        """
        Add a configuration or subspace.

        Adding a subspace that is already covered has no effect. Adding a more
        general subspace removes stored subspaces that become redundant.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B"])
        >>> states.add({"A": 0, "B": 0})
        >>> states.add({"A": "*", "B": 0})
        >>> states.enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 0})

        Parameters
        ----------
        configuration: ConfigurationLike or HypercubeLike
            Configuration or Boolean subspace to add.
        """

        hypercube = self._coerce_hypercube(configuration)
        self._add_hypercube(hypercube)

    def compress(self) -> None:
        """
        Reduce the internal representation without changing the set.

        Compression is opportunistic and does not guarantee a globally minimal
        number of internal subspaces.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B"])
        >>> states.add({"A": 0, "B": 0})
        >>> states.add({"A": 1, "B": 0})
        >>> states.compress()
        >>> states.enumerate()
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 0})
        """

        changed = True
        while changed:
            changed = False
            for i, left in enumerate(self._hypercubes):
                for j in range(i + 1, len(self._hypercubes)):
                    right = self._hypercubes[j]
                    merged = _try_merge_hypercubes(
                        left,
                        right,
                        self._components,
                    )
                    if merged is None:
                        continue

                    self._hypercubes = [
                        hypercube
                        for k, hypercube in enumerate(self._hypercubes)
                        if k not in {i, j}
                    ]
                    self._add_hypercube(merged)
                    changed = True
                    break

                if changed:
                    break

    @property
    def components(self) -> Tuple[str, ...]:
        """
        Return ordered configuration components.
        """

        return self._components

    def contains(
        self,
        configuration: Union[ConfigurationLike, HypercubeLike],
    ) -> bool:
        """
        Test whether a configuration or subspace is fully represented.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B"], [{"A": "*", "B": 0}])
        >>> states.contains({"A": 1, "B": 0})
        True
        >>> states.contains({"B": 0})
        True
        >>> states.contains({"A": 1, "B": 1})
        False
        """

        return self._contains_hypercube(self._coerce_hypercube(configuration))

    def count(self) -> int:
        """
        Return the exact number of represented configurations.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B", "C"], [{"A": 0}])
        >>> states.count()
        4
        """

        return sum(
            _count_hypercube_configurations(hypercube, len(self._components))
            for hypercube in self._hypercubes
        )

    def enumerate(self) -> Tuple[Configuration, ...]:
        """
        Return all represented configurations.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B"], [{"A": 0}])
        >>> states.enumerate()
        ({'A': 0, 'B': 0}, {'A': 0, 'B': 1})
        """

        return tuple(self)

    @overload
    def sample(
        self,
        n: None = None,
        *,
        seed: RandomStateSeed = None,
    ) -> Configuration: ...

    @overload
    def sample(
        self,
        n: int,
        *,
        seed: RandomStateSeed = None,
    ) -> Tuple[Configuration, ...]: ...

    def sample(
        self,
        n: Optional[int] = None,
        *,
        seed: RandomStateSeed = None,
    ) -> Union[Configuration, Tuple[Configuration, ...]]:
        """
        Sample represented configurations with replacement.

        Examples
        --------
        >>> states = ConfigurationSet(["A", "B"], [{"A": "*"}])
        >>> states.sample(seed=0)
        {'A': 0, 'B': 0}
        >>> states.sample(2, seed=0)
        ({'A': 0, 'B': 0}, {'A': 1, 'B': 1})

        Parameters
        ----------
        n: int, optional
            Number of configurations to sample. If omitted, return a single
            configuration.
        seed: int, RandomState or numpy.random, optional
            Random seed.

        Returns
        -------
        Configuration or tuple of Configuration
            One sampled configuration when `n` is omitted, otherwise a tuple of
            sampled configurations.

        Raises
        ------
        ValueError
            If the set is empty or `n` is not positive.
        """

        if not self._hypercubes:
            raise ValueError("cannot sample from an empty ConfigurationSet")

        rng = _as_seed(seed)

        if n is None:
            return self._sample_one(rng)

        n = _as_positive_integer(n, "n")
        return tuple(self._sample_one(rng) for _ in range(n))

    @classmethod
    def _from_encoded_hypercubes(
        cls,
        components: Iterable[str],
        hypercubes: Iterable[_EncodedHypercube],
    ) -> "ConfigurationSet":
        """Build a set from an exact disjoint encoded representation."""

        configurations = cls(components)
        configurations._hypercubes = list(hypercubes)

        return configurations

    def _sample_one(self, rng: object) -> Configuration:

        total = self.count()
        index = int(rng.randint(total))  # pyright: ignore[reportAttributeAccessIssue]

        for hypercube in self._hypercubes:
            size = _count_hypercube_configurations(hypercube, len(self._components))
            if index < size:
                return _configuration_at_index(hypercube, self._components, index)

            index -= size

        raise RuntimeError("failed to sample from ConfigurationSet")

    def _add_hypercube(self, hypercube: _EncodedHypercube) -> None:

        if self._contains_hypercube(hypercube):
            return None

        updated: List[_EncodedHypercube] = []
        for stored in self._hypercubes:
            if _is_encoded_subset(stored, hypercube):
                continue
            updated.extend(
                _subtract_hypercube(stored, hypercube, len(self._components))
            )

        updated.append(hypercube)
        self._hypercubes = updated

    def _contains_hypercube(self, hypercube: _EncodedHypercube) -> bool:

        remaining = [hypercube]

        for stored in self._hypercubes:
            next_remaining = []
            for current in remaining:
                next_remaining.extend(
                    _subtract_hypercube(current, stored, len(self._components))
                )
            remaining = next_remaining
            if not remaining:
                return True

        return False

    def _covers(self, other: "ConfigurationSet") -> bool:

        return all(
            self._contains_hypercube(
                self._coerce_hypercube(_decode_hypercube(hypercube, other._components))
            )
            for hypercube in other._hypercubes
        )

    def _as_hypercubes(self) -> Tuple[Hypercube, ...]:

        return tuple(
            _decode_hypercube(hypercube, self._components)
            for hypercube in self._hypercubes
        )

    def _iter_encoded_hypercubes(self) -> Iterator[_EncodedHypercube]:

        return iter(self._hypercubes)

    def _iter_configuration_bits(self) -> Iterator[int]:

        for hypercube in self._hypercubes:
            yield from _iter_hypercube_configuration_bits(
                hypercube,
                len(self._components),
            )

    def _coerce_hypercube(self, configuration: object) -> _EncodedHypercube:

        if isinstance(configuration, Hypercube):
            hypercube = configuration.copy()
        elif isinstance(configuration, Mapping):
            hypercube = Hypercube(configuration)
        else:
            raise TypeError(
                "unsupported argument type for 'configuration': "
                "expected configuration-like or hypercube-like object "
                f"but received {type(configuration)}"
            )

        unknown = hypercube.components - set(self._components)
        if unknown:
            raise ValueError(
                "configuration contains unknown components: "
                f"{', '.join(sorted(unknown))}"
            )

        return _encode_hypercube(
            _normalize_hypercube(hypercube),
            self._component_to_index,
        )


def _coerce_components(components: Iterable[str]) -> Tuple[str, ...]:

    resolved = tuple(components)

    if not all(isinstance(component, str) for component in resolved):
        raise TypeError("components must be an iterable of strings")

    duplicated = sorted(
        component for component in set(resolved) if resolved.count(component) > 1
    )
    if duplicated:
        raise ValueError(
            f"components must be unique; duplicated components: {', '.join(duplicated)}"
        )

    return resolved


def _normalize_hypercube(hypercube: Hypercube) -> Hypercube:

    return Hypercube(
        {component: value for component, value in hypercube.items() if value.is_fixed}
    )


def _encode_hypercube(
    hypercube: Hypercube,
    component_to_index: Mapping[str, int],
) -> _EncodedHypercube:

    fixed_mask = 0
    value_mask = 0
    for component, value in hypercube.items():
        bit = 1 << component_to_index[component]
        fixed_mask |= bit
        if value.value == 1:
            value_mask |= bit

    return fixed_mask, value_mask


def _count_hypercube_configurations(
    hypercube: _EncodedHypercube,
    n_components: int,
) -> int:

    fixed_mask, _ = hypercube
    return 1 << (n_components - _bit_count(fixed_mask))


def _iter_hypercube_configurations(
    hypercube: _EncodedHypercube,
    components: Sequence[str],
) -> Iterator[Configuration]:

    fixed_mask, value_mask = hypercube
    free_indices = [
        index for index in range(len(components)) if not fixed_mask & (1 << index)
    ]

    for free_values in product((0, 1), repeat=len(free_indices)):
        configuration = {}
        free_position = 0
        for index, component in enumerate(components):
            bit = 1 << index
            if fixed_mask & bit:
                configuration[component] = 1 if value_mask & bit else 0
            else:
                configuration[component] = free_values[free_position]
                free_position += 1
        yield configuration


def _iter_hypercube_configuration_bits(
    hypercube: _EncodedHypercube,
    n_components: int,
) -> Iterator[int]:

    fixed_mask, value_mask = hypercube
    free_mask = ((1 << n_components) - 1) ^ fixed_mask
    free_values = free_mask

    while True:
        yield value_mask | free_values
        if free_values == 0:
            break
        free_values = (free_values - 1) & free_mask


def _configuration_at_index(
    hypercube: _EncodedHypercube,
    components: Sequence[str],
    index: int,
) -> Configuration:

    fixed_mask, value_mask = hypercube
    configuration = {}
    for component_index in reversed(range(len(components))):
        bit = 1 << component_index
        component = components[component_index]
        if fixed_mask & bit:
            configuration[component] = 1 if value_mask & bit else 0
            continue

        configuration[component] = index % 2
        index //= 2

    return {component: configuration[component] for component in components}


def _decode_hypercube(
    hypercube: _EncodedHypercube,
    components: Sequence[str],
) -> Hypercube:

    fixed_mask, value_mask = hypercube
    return Hypercube(
        {
            component: 1 if value_mask & (1 << index) else 0
            for index, component in enumerate(components)
            if fixed_mask & (1 << index)
        }
    )


def _subtract_hypercube(
    hypercube: _EncodedHypercube,
    removed: _EncodedHypercube,
    n_components: int,
) -> Tuple[_EncodedHypercube, ...]:

    intersection = _intersect_hypercubes(hypercube, removed)
    if intersection is None:
        return (hypercube,)

    if _is_encoded_subset(hypercube, removed):
        return ()

    fragments = []
    prefix_fixed = 0
    prefix_value = 0
    fixed_mask, value_mask = hypercube
    intersection_fixed, intersection_value = intersection
    all_mask = (1 << n_components) - 1

    for index in range(n_components):
        bit = 1 << index
        hypercube_values = _component_values(fixed_mask, value_mask, bit)
        intersection_values = _component_values(
            intersection_fixed,
            intersection_value,
            bit,
        )
        excluded_values = hypercube_values - intersection_values

        for value in sorted(excluded_values):
            tail_mask = all_mask ^ ((1 << (index + 1)) - 1)
            tail_fixed = fixed_mask & tail_mask
            fragment_fixed = prefix_fixed | bit | tail_fixed
            fragment_value = (
                prefix_value | (bit if value == 1 else 0) | (value_mask & tail_fixed)
            )
            fragments.append((fragment_fixed, fragment_value & fragment_fixed))

        if len(intersection_values) == 1:
            prefix_fixed |= bit
            if next(iter(intersection_values)) == 1:
                prefix_value |= bit

    return tuple(fragments)


def _intersect_hypercubes(
    left: _EncodedHypercube,
    right: _EncodedHypercube,
) -> Optional[_EncodedHypercube]:

    left_fixed, left_value = left
    right_fixed, right_value = right
    if left_fixed & right_fixed & (left_value ^ right_value):
        return None

    fixed_mask = left_fixed | right_fixed
    value_mask = (left_value & left_fixed) | (right_value & right_fixed)
    return fixed_mask, value_mask & fixed_mask


def _try_merge_hypercubes(
    left: _EncodedHypercube,
    right: _EncodedHypercube,
    components: Sequence[str],
) -> Optional[_EncodedHypercube]:

    differences = []
    fixed_mask = 0
    value_mask = 0
    left_fixed, left_value = left
    right_fixed, right_value = right

    for index in range(len(components)):
        bit = 1 << index
        left_values = _component_values(left_fixed, left_value, bit)
        right_values = _component_values(right_fixed, right_value, bit)

        if left_values == right_values:
            if len(left_values) == 1:
                fixed_mask |= bit
                if next(iter(left_values)) == 1:
                    value_mask |= bit
            continue

        if left_values | right_values == frozenset({0, 1}):
            differences.append(index)
            continue

        return None

    if len(differences) != 1:
        return None

    return fixed_mask, value_mask


def _is_encoded_subset(
    left: _EncodedHypercube,
    right: _EncodedHypercube,
) -> bool:

    left_fixed, left_value = left
    right_fixed, right_value = right
    return (
        left_fixed & right_fixed == right_fixed
        and (left_value ^ right_value) & right_fixed == 0
    )


def _component_values(
    fixed_mask: int,
    value_mask: int,
    bit: int,
) -> frozenset[int]:

    if fixed_mask & bit:
        return frozenset({1 if value_mask & bit else 0})

    return frozenset({0, 1})


def _bit_count(mask: int) -> int:

    return bin(mask).count("1")
