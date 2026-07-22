#!/usr/bin/env python

from __future__ import annotations

import hashlib
import json
import os
import pickle
import tempfile
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Mapping,
    NamedTuple,
    Optional,
    Tuple,
    cast,
)

import pandas as pd
from typing_extensions import Literal

from ...logic.influence_graph import InfluenceGraph
from ..cache import _discard_cached_download, _unlink
from ..cache import path as cache_path

_GRAPH_CACHE_SCHEMA = 1
_PICKLE_PROTOCOL = 4

_NormalizedSign = Literal[-1, 1]
_OrderedEdge = Tuple[Any, Any, _NormalizedSign]


class _GraphSnapshot(NamedTuple):
    graph: InfluenceGraph
    edge_order: Tuple[_OrderedEdge, ...]


class _GraphBuild(NamedTuple):
    snapshot: _GraphSnapshot
    downloads: Tuple[Path, ...] = ()


def _cached_graph(
    key: Mapping[str, Any],
    factory: Callable[[], _GraphBuild],
) -> _GraphSnapshot:
    """Return a fresh canonical graph snapshot from the persistent cache."""

    serialized_key = json.dumps(
        dict(key),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    cache_file = _graph_cache_file(serialized_key)

    snapshot = _read_graph(cache_file, serialized_key)
    if snapshot is not None:
        return snapshot

    build = factory()
    _write_graph(cache_file, serialized_key, build.snapshot)
    for download in dict.fromkeys(build.downloads):
        _discard_cached_download(download)
    return build.snapshot


def _snapshot_from_dataframe(
    dataframe: pd.DataFrame,
    *,
    downloads: Iterable[Path] = (),
) -> _GraphBuild:
    """Construct a graph snapshot while retaining global interaction order."""

    graph = InfluenceGraph.from_dataframe(dataframe)
    edge_order = []
    for source, target, sign in dataframe[["source", "target", "sign"]].itertuples(
        index=False,
        name=None,
    ):
        if bool(pd.isna(sign)):
            continue
        edge_order.append(
            (
                source,
                target,
                cast(_NormalizedSign, InfluenceGraph._normalize_sign(sign)),
            )
        )

    return _GraphBuild(
        snapshot=_GraphSnapshot(graph, tuple(edge_order)),
        downloads=tuple(downloads),
    )


def _ordered_subgraph(
    snapshot: _GraphSnapshot,
    edges: Iterable[_OrderedEdge],
) -> InfluenceGraph:
    """Rebuild selected edges in their original global insertion order."""

    selected_edges = tuple(edges)
    if selected_edges == snapshot.edge_order:
        return snapshot.graph

    source_graph = snapshot.graph
    graph = InfluenceGraph()
    graph.graph.update(source_graph.graph)

    for source, target, sign in selected_edges:
        data = _edge_data(source_graph, source, target, sign).copy()
        data.pop("sign")
        graph.add_edge(source, target, sign=sign, **data)
        graph.nodes[source].update(source_graph.nodes[source])
        graph.nodes[target].update(source_graph.nodes[target])

    return graph


def _edge_data(
    graph: InfluenceGraph,
    source: Any,
    target: Any,
    sign: _NormalizedSign,
) -> Dict[str, Any]:
    """Return edge attributes for one canonical signed influence."""

    for data in graph[source][target].values():
        if data["sign"] == sign:
            return data

    raise ValueError(
        "invalid cached influence graph: missing ordered edge "
        f"{source!r} -> {target!r} with sign {sign!r}"
    )


def _cache_source_identity(value: Any) -> str:
    """Return a stable cache identity for a label or local source file."""

    if isinstance(value, Path) or (
        isinstance(value, str) and value not in ("bundled", "latest")
    ):
        file = Path(value).expanduser()
        try:
            stat = file.stat()
        except OSError:
            return str(file)
        return f"{file.resolve()}:{stat.st_size}:{stat.st_mtime_ns}"

    return str(value)


def _graph_cache_file(serialized_key: str) -> Path:

    digest = hashlib.sha256(serialized_key.encode("utf-8")).hexdigest()
    return (
        cache_path("omnipath")
        / "graphs"
        / f"v{_GRAPH_CACHE_SCHEMA}-{digest}.pickle"
    )


def _read_graph(
    cache_file: Path,
    serialized_key: str,
) -> Optional[_GraphSnapshot]:

    try:
        with cache_file.open("rb") as handle:
            payload = pickle.load(handle)

        if not isinstance(payload, dict):
            raise ValueError("invalid graph cache payload")
        if payload.get("schema") != _GRAPH_CACHE_SCHEMA:
            raise ValueError("unsupported graph cache schema")
        if payload.get("key") != serialized_key:
            raise ValueError("graph cache key mismatch")

        graph = payload.get("graph")
        edge_order = payload.get("edge_order")
        if not isinstance(graph, InfluenceGraph) or not isinstance(edge_order, tuple):
            raise ValueError("invalid graph cache contents")

        return _GraphSnapshot(graph, edge_order)
    except FileNotFoundError:
        return None
    except Exception:
        _unlink(cache_file)
        return None


def _write_graph(
    cache_file: Path,
    serialized_key: str,
    snapshot: _GraphSnapshot,
) -> None:

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{cache_file.name}.",
        suffix=".tmp",
        dir=str(cache_file.parent),
    )
    temporary_file = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            pickle.dump(
                {
                    "schema": _GRAPH_CACHE_SCHEMA,
                    "key": serialized_key,
                    "graph": snapshot.graph,
                    "edge_order": snapshot.edge_order,
                },
                handle,
                protocol=_PICKLE_PROTOCOL,
            )
            handle.flush()
            os.fsync(handle.fileno())
        temporary_file.replace(cache_file)
    finally:
        _unlink(temporary_file)
