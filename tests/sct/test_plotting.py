import bonesistools as bt

from matplotlib.axes import Axes
from matplotlib.figure import Figure

bt.sct.pl.set_default_params()

ADATA = bt.sct.datasets.nestorowa()


def test_embedding_plot_with_test_representation():
    adata = ADATA.copy()

    adata.obsm["X_test"] = adata.X[:, :2].copy()
    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    fig, ax = bt.sct.pl.embedding_plot(
        adata,
        obs="cluster",
        use_rep="X_test",
        xlabel="X1",
        ylabel="X2",
        figwidth=6,
        s=2,
        alpha=1,
        add_legend=True,
        lgd_params={
            "title": "clusters",
            "ncol": 1,
            "markerscale": 5,
            "frameon": True,
            "edgecolor": bt.sct.pl.get_color("black"),
            "shadow": False,
        },
        n_components=2,
        background_visible=False,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_kde_plot_returns_matplotlib_objects():
    adata = ADATA.copy()

    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    gene = adata.var_names[0]

    fig, ax = bt.sct.pl.kde_plot(
        adata,
        gene=gene,
        xlabel="count",
        ylabel="density",
        show_legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_ecdf_plot_returns_matplotlib_objects():
    adata = ADATA.copy()

    adata.obs["cluster"] = adata.obs["clusters"].astype("category")

    gene = adata.var_names[0]

    fig, ax = bt.sct.pl.ecdf_plot(
        adata,
        gene=gene,
        xlabel="count",
        ylabel="ECDF",
        show_legend=True,
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
