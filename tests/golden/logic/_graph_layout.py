#!/usr/bin/env python

from pathlib import Path

import bonesistools as bt

EXPECTED_PATH = Path(__file__).parent / "expected" / "aggregated_influence_graph.dot"


def graphviz_input() -> str:
    networks = [
        {
            "APEX": "BCL2 | MCL1",
            "BCL2": "MCL1 & ~CASP3",
            "CASP3": "BAX & ~BCL2",
            "BAX": "TP53 | CASP3",
            "MCL1": "AKT & ~TP53",
            "AKT": "PI3K",
            "G1": "APEX",
            "G2": "APEX",
            "OUT": "G1 | G2",
            "PI3K": "APEX",
            "TP53": "CASP3 | MCL1",
        },
        {
            "APEX": "BCL2 & MCL1",
            "BCL2": "MCL1 | ~CASP3",
            "CASP3": "BAX & ~BCL2",
            "BAX": "TP53 & CASP3",
            "MCL1": "AKT | ~TP53",
            "AKT": "PI3K",
            "G1": "APEX",
            "G2": "APEX",
            "OUT": "G1 | G2",
            "PI3K": "APEX",
            "TP53": "CASP3 & MCL1",
        },
    ]
    graph = bt.logic.bn.BooleanNetworkEnsemble(*networks).to_influence_graph()

    return graph.to_pydot(
        edge_label=None,
        node_style=None,
        edge_style=None,
    ).to_string()


def save_expected(graphviz: str) -> None:
    EXPECTED_PATH.write_text(graphviz)
