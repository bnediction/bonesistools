import bonesistools as bt

ADATA = bt.sct.datasets.nestorowa()


def test_mitochondrial_and_ribosomal_gene_classification():

    adata = ADATA.copy()

    bt.sct.tl.mitochondrial_genes(
        adata,
        index_type="name",
        key="mt",
        axis=1,
        copy=False,
    )

    bt.sct.tl.ribosomal_genes(
        adata,
        index_type="name",
        key="rps",
        axis=1,
        copy=False,
    )

    assert "mt" in adata.var
    assert "rps" in adata.var
    assert adata.var["mt"].dtype == bool
    assert adata.var["rps"].dtype == bool


def test_filter_var_removes_genes():

    adata = ADATA.copy()
    initial_n_vars = adata.n_vars

    adata.var["test_filter"] = False
    adata.var.iloc[:10, adata.var.columns.get_loc("test_filter")] = True

    bt.sct.pp.filter_var(
        adata,
        "test_filter",
        lambda x: x,
    )

    assert adata.n_vars == 10
    assert adata.n_vars < initial_n_vars


def test_filter_obs_removes_cells():

    adata = ADATA.copy()
    initial_n_obs = adata.n_obs

    adata.obs["test_filter"] = False
    adata.obs.iloc[:20, adata.obs.columns.get_loc("test_filter")] = True

    bt.sct.pp.filter_obs(
        adata,
        "test_filter",
        lambda x: x,
    )

    assert adata.n_obs == 20
    assert adata.n_obs < initial_n_obs
