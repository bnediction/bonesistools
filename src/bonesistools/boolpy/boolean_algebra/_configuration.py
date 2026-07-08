#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple, Union, overload

from ..._typing import RandomStateSeed
from ..._validation import _as_positive_integer, _as_seed
from ._boolean import PartialBoolean
from ._hypercube import Hypercube
from ._typing import ConfigurationLike, HypercubeLike, PartialBooleanLike


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
    internal subspaces.

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

    __slots__ = ("_components", "_hypercubes")

    def __init__(
        self,
        components: Iterable[str],
        configurations: Optional[Iterable[HypercubeLike]] = None,
    ) -> None:

        self._components = _coerce_components(components)
        self._hypercubes: List[Hypercube] = []

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

        return self.contains(hypercube)

    def __iter__(self) -> Iterator[Dict[str, int]]:
        """
        Iterate over complete configurations.

        Yields
        ------
        dict[str, int]
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
        copied._hypercubes = [hypercube.copy() for hypercube in self._hypercubes]

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

        if self.contains(hypercube):
            return None

        updated: List[Hypercube] = []
        for stored in self._hypercubes:
            if stored <= hypercube:
                continue
            updated.extend(_subtract_hypercube(stored, hypercube, self._components))

        updated.append(hypercube)
        self._hypercubes = updated

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
                    self.add(merged)
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

        hypercube = self._coerce_hypercube(configuration)
        remaining = [hypercube]

        for stored in self._hypercubes:
            next_remaining = []
            for current in remaining:
                next_remaining.extend(
                    _subtract_hypercube(current, stored, self._components)
                )
            remaining = next_remaining
            if not remaining:
                return True

        return False

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
            _count_hypercube_configurations(hypercube, self._components)
            for hypercube in self._hypercubes
        )

    def enumerate(self) -> Tuple[Dict[str, int], ...]:
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
    ) -> Dict[str, int]: ...

    @overload
    def sample(
        self,
        n: int,
        *,
        seed: RandomStateSeed = None,
    ) -> Tuple[Dict[str, int], ...]: ...

    def sample(
        self,
        n: Optional[int] = None,
        *,
        seed: RandomStateSeed = None,
    ) -> Union[Dict[str, int], Tuple[Dict[str, int], ...]]:
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
        dict[str, int] or tuple of dict[str, int]
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

    def _sample_one(self, rng: object) -> Dict[str, int]:

        total = self.count()
        index = int(rng.randint(total))  # pyright: ignore[reportAttributeAccessIssue]

        for hypercube in self._hypercubes:
            size = _count_hypercube_configurations(hypercube, self._components)
            if index < size:
                return _configuration_at_index(hypercube, self._components, index)

            index -= size

        raise RuntimeError("failed to sample from ConfigurationSet")

    def _coerce_hypercube(self, configuration: object) -> Hypercube:

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

        return _normalize_hypercube(hypercube)


def _coerce_components(components: Iterable[str]) -> Tuple[str, ...]:

    resolved = tuple(components)

    if not all(isinstance(component, str) for component in resolved):
        raise TypeError("components must be an iterable of strings")

    duplicated = sorted(
        component for component in set(resolved) if resolved.count(component) > 1
    )
    if duplicated:
        raise ValueError(
            "components must be unique; duplicated components: "
            f"{', '.join(duplicated)}"
        )

    return resolved


def _normalize_hypercube(hypercube: Hypercube) -> Hypercube:

    return Hypercube(
        {component: value for component, value in hypercube.items() if value.is_fixed}
    )


def _count_hypercube_configurations(
    hypercube: Hypercube,
    components: Sequence[str],
) -> int:

    count = 1
    for component in components:
        count *= len(hypercube._get_value(component).as_set())

    return count


def _iter_hypercube_configurations(
    hypercube: Hypercube,
    components: Sequence[str],
) -> Iterator[Dict[str, int]]:

    values = [
        sorted(hypercube._get_value(component).as_set()) for component in components
    ]

    for configuration in product(*values):
        yield dict(zip(components, configuration))


def _configuration_at_index(
    hypercube: Hypercube,
    components: Sequence[str],
    index: int,
) -> Dict[str, int]:

    configuration = {}
    for component in reversed(components):
        values = sorted(hypercube._get_value(component).as_set())
        value = values[index % len(values)]
        configuration[component] = value
        index //= len(values)

    return {component: configuration[component] for component in components}


def _subtract_hypercube(
    hypercube: Hypercube,
    removed: Hypercube,
    components: Sequence[str],
) -> Tuple[Hypercube, ...]:

    intersection = _intersect_hypercubes(hypercube, removed, components)
    if intersection is None:
        return (hypercube.copy(),)

    if hypercube <= removed:
        return ()

    fragments = []
    prefix: Dict[str, PartialBooleanLike] = {}

    for index, component in enumerate(components):
        hypercube_values = hypercube._get_value(component).as_set()
        intersection_values = intersection._get_value(component).as_set()
        excluded_values = hypercube_values - intersection_values

        for value in sorted(excluded_values):
            fragment = dict(prefix)
            fragment[component] = value
            for tail_component in components[index + 1 :]:
                tail_value = hypercube._get_value(tail_component)
                if tail_value.is_fixed:
                    fragment[tail_component] = tail_value.value
            fragments.append(Hypercube(fragment))

        if len(intersection_values) == 1:
            prefix[component] = next(iter(intersection_values))

    return tuple(fragments)


def _intersect_hypercubes(
    left: Hypercube,
    right: Hypercube,
    components: Sequence[str],
) -> Optional[Hypercube]:

    values = {}
    for component in components:
        intersection = (
            left._get_value(component).as_set() & right._get_value(component).as_set()
        )
        if not intersection:
            return None
        if len(intersection) == 1:
            values[component] = next(iter(intersection))

    return Hypercube(values)


def _try_merge_hypercubes(
    left: Hypercube,
    right: Hypercube,
    components: Sequence[str],
) -> Optional[Hypercube]:

    differences = []
    values: Dict[str, PartialBoolean] = {}

    for component in components:
        left_value = left._get_value(component)
        right_value = right._get_value(component)

        if left_value == right_value:
            if left_value.is_fixed:
                values[component] = left_value
            continue

        if left_value.as_set() | right_value.as_set() == frozenset({0, 1}):
            differences.append(component)
            continue

        return None

    if len(differences) != 1:
        return None

    return Hypercube(values)
