#!/usr/bin/env python

import networkx as nx
import pytest

import bonesistools as bt
from bonesistools.boolpy.influence_graph import _scoring


def test_influence_graph_scoring_is_exported_from_public_namespace():
    score = bt.bpy.ig.InteractionScore(score=0.0, total_weight=0.0, path_number=0)

    assert score.normalized_score == 0.0
    assert callable(bt.bpy.ig.infer_signed_interactions)
    assert callable(bt.bpy.ig.infer_signed_interactions_from_walks)
    assert callable(bt.bpy.ig.interaction_scores_from_walks)


def test_interaction_scores_from_walks_cancels_opposite_signed_paths():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "D", sign=1)
    graph.add_edge("A", "C", sign=1)
    graph.add_edge("C", "D", sign=-1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "D"],
        max_depth=2,
    )

    score = scores["A"]["D"]

    assert score.path_number == 2
    assert score.score == 0.0
    assert score.total_weight == 0.5
    assert score.normalized_score == 0.0


def test_interaction_scores_from_walks_detects_majority_positive_influence():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "D", sign=1)
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "D", sign=1)
    graph.add_edge("A", "C", sign=1)
    graph.add_edge("C", "D", sign=-1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "D"],
        max_depth=2,
    )

    score = scores["A"]["D"]

    assert score.path_number == 3
    assert score.score == 1.0
    assert score.total_weight == 1.5
    assert score.normalized_score == pytest.approx(2 / 3)


def test_interaction_scores_from_walks_expands_parallel_signed_edges():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("A", "B", sign=-1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "B"],
        max_depth=1,
    )

    score = scores["A"]["B"]

    assert score.path_number == 2
    assert score.score == 0.0
    assert score.total_weight == 2.0
    assert score.normalized_score == 0.0


def test_interaction_scores_from_walks_ignores_self_scores():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "A", sign=1)
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "A", sign=-1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "B"],
        max_depth=2,
    )

    assert "A" not in scores["A"]
    assert "B" not in scores["B"]
    assert scores["A"]["B"].path_number == 2


def test_scores_from_walks_handles_digraphs_and_missing_sources():
    graph = nx.DiGraph()
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "C", sign=1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "C", "missing"],
        max_depth=2,
        weights=[1.0, 0.5],
    )

    assert scores["A"]["C"] == bt.bpy.ig.InteractionScore(
        score=0.5,
        total_weight=0.5,
        path_number=1,
    )
    assert scores["missing"] == {}
    assert _scoring._edge_signs(graph, "C", "A") == []
    assert _scoring._walk_signs(graph, ["A", "C"]) == []

    with pytest.raises(ValueError, match="expected at least 2 values"):
        bt.bpy.ig.interaction_scores_from_walks(
            graph,
            genes=["A", "C"],
            max_depth=2,
            weights=[1.0],
        )


def test_infer_signed_interactions_from_walks_matches_two_step_pipeline():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "D", sign=1)
    graph.add_edge("A", "B", sign=1)
    graph.add_edge("B", "D", sign=1)
    graph.add_edge("A", "C", sign=1)
    graph.add_edge("C", "D", sign=-1)

    scores = bt.bpy.ig.interaction_scores_from_walks(
        graph,
        genes=["A", "D"],
        max_depth=2,
    )
    expected = bt.bpy.ig.infer_signed_interactions(
        scores,
        genes=["A", "D"],
        threshold=0.6,
    )

    assert (
        bt.bpy.ig.infer_signed_interactions_from_walks(
            graph,
            genes=["A", "D"],
            max_depth=2,
            threshold=0.6,
        )
        == expected
        == [("A", "D", {"sign": 1})]
    )


def test_infer_signed_interactions_keeps_clear_dominant_direction():
    scores = {
        "A": {
            "B": bt.bpy.ig.InteractionScore(
                score=3.0,
                total_weight=4.0,
                path_number=4,
            )
        },
        "B": {
            "A": bt.bpy.ig.InteractionScore(
                score=-1.0,
                total_weight=4.0,
                path_number=4,
            )
        },
    }

    assert bt.bpy.ig.infer_signed_interactions(
        scores,
        genes=["A", "B"],
        threshold=0.75,
    ) == [("A", "B", {"sign": 1})]


def test_infer_signed_interactions_accepts_threshold_boundaries():
    positive_scores = {
        "A": {
            "B": bt.bpy.ig.InteractionScore(
                score=1.0,
                total_weight=1.0,
                path_number=1,
            )
        },
        "B": {},
    }
    negative_scores = {
        "A": {
            "B": bt.bpy.ig.InteractionScore(
                score=-0.25,
                total_weight=1.0,
                path_number=1,
            )
        },
        "B": {},
    }

    assert bt.bpy.ig.infer_signed_interactions(
        positive_scores,
        genes=["A", "B"],
        threshold=1,
    ) == [("A", "B", {"sign": 1})]
    assert bt.bpy.ig.infer_signed_interactions(
        negative_scores,
        genes=["A", "B"],
        threshold=0,
    ) == [("A", "B", {"sign": -1})]


def test_infer_signed_interactions_can_keep_bidirectional_edges():
    scores = {
        "A": {
            "B": bt.bpy.ig.InteractionScore(
                score=4.0,
                total_weight=4.0,
                path_number=4,
            )
        },
        "B": {
            "A": bt.bpy.ig.InteractionScore(
                score=-4.0,
                total_weight=4.0,
                path_number=4,
            )
        },
    }

    assert bt.bpy.ig.infer_signed_interactions(
        scores,
        genes=["A", "B"],
        threshold=0.75,
        allow_bidirectional=True,
    ) == [
        ("A", "B", {"sign": 1}),
        ("B", "A", {"sign": -1}),
    ]


def test_infer_signed_interactions_rejects_weak_or_ambiguous_scores():
    scores = {
        "A": {
            "B": bt.bpy.ig.InteractionScore(
                score=3.0,
                total_weight=4.0,
                path_number=4,
            )
        },
        "B": {
            "A": bt.bpy.ig.InteractionScore(
                score=3.0,
                total_weight=4.0,
                path_number=4,
            )
        },
        "C": {
            "D": bt.bpy.ig.InteractionScore(
                score=1.0,
                total_weight=1.0,
                path_number=1,
            )
        },
        "D": {
            "C": bt.bpy.ig.InteractionScore(
                score=0.0,
                total_weight=0.0,
                path_number=0,
            )
        },
    }

    assert (
        bt.bpy.ig.infer_signed_interactions(
            scores,
            genes=["A", "B"],
            threshold=0.75,
        )
        == []
    )
    assert (
        bt.bpy.ig.infer_signed_interactions(
            scores,
            genes=["C", "D"],
            threshold=0.75,
            minimum_path_number=2,
        )
        == []
    )


def test_interaction_scoring_rejects_invalid_arguments_and_signs():
    graph = nx.MultiDiGraph()
    graph.add_edge("A", "B", sign=0)

    with pytest.raises(ValueError):
        bt.bpy.ig.interaction_scores_from_walks(
            graph,
            genes=["A", "B"],
            max_depth=1,
        )

    with pytest.raises(ValueError):
        bt.bpy.ig.interaction_scores_from_walks(
            graph,
            genes=["A", "B"],
            max_depth=0,
        )

    with pytest.raises(ValueError):
        bt.bpy.ig.interaction_scores_from_walks(
            graph,
            genes=["A", "B"],
            max_depth=2,
            weights=[1.0],
        )

    with pytest.raises(ValueError):
        bt.bpy.ig.infer_signed_interactions({}, threshold=2)

    with pytest.raises(ValueError):
        bt.bpy.ig.infer_signed_interactions({}, minimum_path_number=0)
