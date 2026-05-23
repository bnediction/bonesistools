#!/usr/bin/env python

try:
    from collections import Sequence as Seq
except ImportError:
    from collections.abc import Sequence as Seq

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from typing import Union, Optional, Iterable, Sequence, Callable
from .._typing import anndata_checker

import pandas as pd
import numpy as np

from anndata import AnnData
from numpy import log2

from ._conversion import anndata_to_dataframe

_CorrMethod = Literal["benjamini-hochberg", "bonferroni"]
_Alternatives = Literal["two-sided", "less", "greater"]


@anndata_checker
def calculate_logfoldchanges(
    adata: AnnData,
    groupby: str,
    layer: Optional[str] = None,
    column_name: str = "logfoldchanges",
    is_log: bool = False,
    cluster_rebalancing: bool = False,
    filter_logfoldchanges: Optional[Callable] = None,
) -> pd.DataFrame:
    """
    Compute log2 fold changes between each group and the remaining observations.

    The computation follows the fold-change convention used by Seurat's
    FindAllMarkers rather than Scanpy's `rank_gene_groups` output.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation key in `adata.obs` defining groups.
    layer: str, optional
        Layer to use instead of `adata.X`.
    column_name: str (default: "logfoldchanges")
        Column name in output dataframe storing log2 fold-change values.
    is_log: bool (default: False)
        Whether values are already log1p-transformed.
    cluster_rebalancing: bool (default: False)
        Whether to average groups before computing fold changes, so groups
        contribute equally regardless of size.
    filter_logfoldchanges: Callable, optional
        Callable used to filter rows from the resulting log-fold-change table.

    Returns
    -------
    DataFrame
        Table with `group`, `names` and `<column_name>` columns.

    See Also
    --------
    Get more information about the difference between log2 fold-changes derived with Seurat and Scanpy here:
    <https://www.biostars.org/p/453129/>

    Raises
    ------
    TypeError
        If `filter_logfoldchanges` is specified but is not callable.
    """

    def compute_logfc(mean_in, mean_out, cluster):

        __df = pd.DataFrame(log2(mean_in) - log2(mean_out), columns=[column_name])
        __df.reset_index(names="names", inplace=True)
        __df.insert(0, "group", cluster)
        return __df

    logfoldchanges = []
    counts_df = anndata_to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log)

    if cluster_rebalancing:
        mean_counts_df = counts_df.groupby(by=groupby, sort=True).mean()

        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = mean_counts_df.loc[cluster]
            _mean_out = mean_counts_df.drop(index=cluster, inplace=False).mean()

            logfoldchanges.append(compute_logfc(_mean_in, _mean_out, cluster))

    else:
        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = counts_df.loc[
                counts_df[groupby] == cluster,
                counts_df.columns != groupby,
            ].mean()

            _mean_out = counts_df.loc[
                counts_df[groupby] != cluster,
                counts_df.columns != groupby,
            ].mean()

            logfoldchanges.append(compute_logfc(_mean_in, _mean_out, cluster))

    if len(logfoldchanges) == 0:
        return pd.DataFrame(columns=["group", "names", column_name])

    logfoldchanges_df = pd.concat(logfoldchanges, ignore_index=True)

    if filter_logfoldchanges is not None:
        if not callable(filter_logfoldchanges):
            raise TypeError(
                f"unsupported argument type for 'filter_logfoldchanges': "
                f"expected callable object but received {type(filter_logfoldchanges)}"
            )
        else:
            logfoldchanges_df = logfoldchanges_df.loc[
                filter_logfoldchanges(logfoldchanges_df[column_name].values)
            ]

    return logfoldchanges_df.reset_index(drop=True)


def hypergeometric_test(
    adata: AnnData,
    signature: Sequence[str],
    markers: Sequence[str],
) -> float:
    """
    Test marker/signature overlap with a hypergeometric survival function.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    signature: Sequence[str]
        Signature genes for a cell type.
    markers: Sequence[str]
        Marker genes for a cluster.

    Returns
    -------
    float
        Hypergeometric survival-function p-value.

    Raises
    ------
    TypeError
        If `adata` is not an AnnData object.
    """

    from scipy.stats import hypergeom

    if not isinstance(adata, AnnData):
        raise TypeError(
            f"unsupported argument type for 'adata': "
            f"expected {AnnData} but received {type(adata)}"
        )

    background = set(adata.var.index)
    if not isinstance(signature, set):
        signature = set(signature)
    if not isinstance(markers, set):
        markers = set(markers)
    marked_genes = markers.intersection(signature)

    N = len(background)  # population size
    K = len(signature)  # number of success states
    n = len(markers)  # number of draws
    k = len(marked_genes)  # number of observed successes (matching genes)

    return hypergeom.sf(k=k, M=N, n=K, N=n, loc=1)


def smirnov_tests(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    layer: Optional[str] = None,
    alternative: _Alternatives = "two-sided",
    corr_method: _CorrMethod = "benjamini-hochberg",
    pval_cutoff: Optional[float] = None,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Compare whether a subsample and a reference sample have the same distribution
    using a two-sample Kolmogorov-Smirnov test.

    Expects logarithmized data.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation key in `adata.obs` defining groups.
    groups: 'all' or sequence of str (default: "all")
        Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparisons
        shall be restricted, or 'all' for performing comparisons for all groups.
    reference: 'rest' | str
        Name of the group used as reference.
        If 'rest', compare each group to all observations not in the group.
    layer: str, optional
        Layer to use instead of `adata.X`.
    alternative: 'two-sided' | 'less' | 'greater' (default: two-sided)
        Defines the null and alternative hypotheses.
    corr_method: 'benjamini-hochberg' | 'bonferroni'
        Correction method for computing adjusted p-values.
    pval_cutoff: float (optional, default: None)
        Return only adjusted p-values below the cutoff.
    key_added: str, optional
        Key in `adata.uns` where information is saved.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    AnnData or None
        Copy of `adata` with results if `copy=True`; otherwise None after
        modifying `adata`.

    Notes
    -----
    Results are stored in `adata.uns[key_added]` with the following values:
    `group`, `names`, `statistics`, `locations`, `signs`, `pvals` and
    `pvals_adj`.

    Raises
    ------
    TypeError
        If `groups` is neither `"all"` nor a non-string sequence.
    ValueError
        If `reference` is neither `"rest"` nor one of the categories in
        `adata.obs[groupby]`.
    """

    from scipy.sparse import issparse
    from scipy.stats import kstest

    adata = adata.copy() if copy else adata

    if groups == "all":
        groups = sorted(list(adata.obs[groupby].cat.categories))
    elif isinstance(groups, Seq) and not isinstance(groups, str):
        groups = list(groups)
    else:
        raise TypeError(
            f"unsupported argument type for 'groups': "
            f"expected sequence but received {type(groups)}"
        )

    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        cats = sorted(list(adata.obs[groupby].cat.categories))
        raise ValueError(
            f"invalid argument value for 'reference': "
            f"expected 'rest' or one of {cats} but received {reference!r}"
        )

    if key_added is None:
        key_added = "smirnov_tests"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {
        "groupby": groupby,
        "reference": reference,
        "layer": layer,
        "corr_method": corr_method,
    }

    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if issparse(X):
        X.eliminate_zeros()

    def get_gene_values(obs_mask, gene_name):
        obs_indices = np.where(np.asarray(obs_mask))[0]
        var_indices = np.where(np.asarray(adata.var.index == gene_name))[0]
        if issparse(X):
            values = X[obs_indices][:, var_indices].toarray()
        else:
            values = X[np.ix_(obs_indices, var_indices)]
        return np.asarray(values).reshape(-1)

    df = pd.DataFrame(
        columns=["group", "names", "statistics", "locations", "signs", "pvals"]
    )
    index = 0

    for name in adata.var_names:
        if reference != "rest":
            ref_sample = get_gene_values(adata.obs[groupby] == reference, name)
        for group in groups:
            if reference == "rest":
                ref_sample = get_gene_values(adata.obs[groupby] != group, name)
            sample = get_gene_values(adata.obs[groupby] == group, name)
            ks = kstest(sample, ref_sample, alternative=alternative)
            df.loc[index] = [
                group,
                name,
                ks.statistic,
                ks.statistic_location,
                ks.statistic_sign,
                ks.pvalue,
            ]
            index += 1

    df.sort_values(by=["statistics"], ascending=False, inplace=True, ignore_index=True)

    if corr_method == "benjamini-hochberg":
        pvals = df["pvals"].to_numpy(dtype=float)
        if pvals.size:
            order = np.argsort(pvals)
            ranks = np.arange(1, pvals.size + 1)
            adjusted = pvals[order] * pvals.size / ranks
            adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
            pvals_adj = np.empty_like(adjusted)
            pvals_adj[order] = np.minimum(adjusted, 1.0)
            df["pvals_adj"] = pvals_adj
        else:
            df["pvals_adj"] = []
    elif corr_method == "bonferroni":
        n_genes = adata.var.shape[0]
        df["pvals_adj"] = np.minimum(df["pvals"] * n_genes, 1.0)

    if pval_cutoff is not None:
        df = df[df["pvals_adj"] < pval_cutoff]

    if copy:
        return df
    else:
        adata.uns[key_added]["results"] = df
