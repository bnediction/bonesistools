#!/usr/bin/env python

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import pytest
from scipy.sparse import csr_matrix

import bonesistools as bt


def _add_paga_edges(adata):
    adata.uns["paga_edges"] = csr_matrix([[0.0, 0.2], [0.4, 0.0]])


@pytest.mark.parametrize("use_rep", ["X_paga_2d", "X_pca"])
def test_draw_paga_reuses_axes_with_2d_and_3d_representations(mini_adata, use_rep):
    _add_paga_edges(mini_adata)
    mini_adata.obsm["X_paga_2d"] = mini_adata.obsm["X_pca"][:, :2].copy()

    fig, ax = plt.subplots()
    returned_ax = bt.sct.pl.draw_paga(
        mini_adata,
        obs="cluster",
        use_rep=use_rep,
        edges="paga_edges",
        threshold=0.1,
        ax=ax,
        with_labels=True,
        node_color={"A": "red", "B": "blue"},
    )

    assert returned_ax is ax
    assert isinstance(returned_ax, Axes)
    plt.close(fig)


def test_draw_paga_writes_outfile(mini_adata, tmp_path):
    _add_paga_edges(mini_adata)

    outfile = tmp_path / "paga.png"
    result = bt.sct.pl.draw_paga(
        mini_adata,
        obs="cluster",
        use_rep="X_pca",
        edges="paga_edges",
        outfile=outfile,
    )

    assert result is None
    assert outfile.exists()
