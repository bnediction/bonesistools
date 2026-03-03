#!/usr/bin/env python

try:
    from collections import Sequence as Seq
except:
    from _collections_abc import Sequence as Seq

from typing import (
    Union,
    Optional,
    Iterable,
    Sequence,
    Literal,
    Callable
)
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
    filter_logfoldchanges: Optional[Callable] = None
) -> pd.DataFrame:
    """
    Log2 fold-change is a metric translating how much the transcript's expression
    has changed between cells in and out of a cluster. The reported values are based
    on a logarithmic scale to base 2 with respect to the fold change ratios.
    According to <https://www.biostars.org/p/453129/>, computed log2 fold changes
    are different between FindAllMarkers (package Seurat) and rank_gene_groups
    (module Scanpy) functions. As mentionned, results derived from Scanpy are inconsistent.
    This function computes it in the right way, with identical results to Seurat by keeping default options.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    groupby: str
        Any key in 'adata.obs' corresponding to defined clusters or groups.
    layer: str (optional, default: None)
        Any key in 'adata.layers'.
        If not specify, log2 fold changes are derived from 'adata.X'.
    column_name: str (default: 'logfoldchanges')
        Column name in output dataframe storing log2 fold-change values.
    is_log: bool (default: False)
        Specify whether count matrix is logarithmized or not.
    cluster_rebalancing: bool (default: False)
        If no cluster rebalancing, cells are equally-weighted.
        Otherwise, cells are weighted with cluster size such as clusters are equally-weighted.
        It means that cells in small cluster have a greater weight than other cells in order
        to correct cluster size effects.
    filter_logfoldchanges: Function (optional, default: None)
        Function filtering results with respect log2 fold-change values.
    
    Returns
    -------
    Return DataFrame storing following values:
    - **group**: group names.
    - **names**: gene names or gene ids.
    - **<column_name>**: log2 fold-change values.
    
    See also
    --------
    Get more information about the difference between log2 fold-changes derived with Seurat and Scanpy here:
    <https://www.biostars.org/p/453129/>
    """

    def compute_logfc(mean_in, mean_out, cluster):

        __df = pd.DataFrame(log2(mean_in) - log2(mean_out), columns=[column_name])
        __df.reset_index(names="names", inplace=True)
        __df.insert(0, "group", cluster)
        return __df
    
    logfoldchanges_df = pd.DataFrame(columns=["group","names",column_name])
    counts_df = anndata_to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log)

    if cluster_rebalancing:
        mean_counts_df = counts_df.groupby(by=groupby, sort=True).mean()
        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = mean_counts_df.loc[cluster]
            _mean_out = mean_counts_df.drop(index=cluster, inplace=False).mean()
            _logfoldchanges_df = compute_logfc(_mean_in, _mean_out, cluster)
            logfoldchanges_df = pd.concat([logfoldchanges_df, _logfoldchanges_df.copy()])
    else:
        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = counts_df.loc[counts_df[groupby] == cluster, counts_df.columns != groupby].mean()
            _mean_out = counts_df.loc[counts_df[groupby] != cluster, counts_df.columns != groupby].mean()
            _logfoldchanges_df = compute_logfc(_mean_in, _mean_out, cluster)
            logfoldchanges_df = pd.concat([logfoldchanges_df, _logfoldchanges_df.copy()])
            del _logfoldchanges_df
    
    if filter_logfoldchanges is not None:
        if not callable(filter_logfoldchanges):
            raise TypeError(f"unsupported argument type for 'filter_logfoldchanges': expected callable object")
        else:
            logfoldchanges_df = logfoldchanges_df.loc[filter_logfoldchanges(logfoldchanges_df[column_name].values)]

    return logfoldchanges_df.reset_index(drop=True)

def hypergeometric_test(
    adata: AnnData,
    signature: Sequence[str],
    markers: Sequence[str],
) -> float:
    """
    Estimates the p-value (or survival function) of an hypergeometric distribution in order
    to test whether marker genes significantly match signature genes.
    Given a population size N and a number of sccess states K, it describes the probability
    of having at least k successes in n draws, without replacement, where:
    - N is the number of genes in anndata,
    - K is the number of signature genes,
    - n is the number of markers,
    - k is the number of gene matching both signature genes and markers.
    Smaller the p-value, higher the probability that genes of the given
    cluster comes from the cell-type associated to the given signature.    

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    signature: Sequence[str]
        Set of signature genes in a given cell-type.
        A signature is a set of over-expressed genes in a cell-type.
    markers: Sequence[str]
        Set of markers (genes) in a given cluster.
        A marker set is a set of over-expressed genes in a cluster.
    
    Returns
    -------
    Return the p-value.
    """

    from scipy.stats import hypergeom

    if not isinstance(adata, AnnData):
        raise TypeError(f"unsupported argument type for 'adata': expected {AnnData} but received {type(adata)}")
    
    background = set(adata.var.index)
    if not isinstance(signature, set):
        signature = set(signature)
    if not isinstance(markers, set):
        markers = set(markers)
    marked_genes = markers.intersection(signature)

    N = len(background)         # population size
    K = len(signature)          # number of success states
    n = len(markers)            # number of draws
    k = len(marked_genes)       # number of observed successes (matching genes)
    
    return hypergeom.sf(
        k=k,
        M=N,
        n=K,
        N=n,
        loc=1
    )

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
    copy: bool = False
) -> Optional[AnnData]:
    """
    Compare whether a subsample and a reference sample have the same distribution
    using a two-sample Kolmogorov-Smirnov test.

    Expects logarithmized data.

    Parameters
    ----------
    adata: ad.AnnData
        Unimodal annotated data matrix.
    groupby: str
        Any key in 'adata.obs' corresponding to defined groups to consider.
    groups: 'all' | Sequence (default: 'all')
        Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparisons
        shall be restricted, or 'all' for performing comparisons for all groups.
    reference: 'rest' | str
        Name of the group used as reference.
        If 'rest', compare each group to all observations not in the group.
    layer: str (optional, default: None)
        Any key in 'adata.layers' whose value will be used to perform tests on.
        If not specify, tests are derived from 'adata.X'.
    alternative: 'two-sided' | 'less' | 'greater' (default: two-sided)
        Defines the null and alternative hypotheses.
    corr_method: 'benjamini-hochberg' | 'bonferroni'
        Correction method for computing adjusted p-values.
    pval_cutoff: float (optional, default: None)
        Return only adjusted p-values below the cutoff.
    key_added
        Key in 'adata.uns' where information is saved to.
    copy: bool (default: False)
        Return a copy instead of updating 'adata' object.
    
    Returns
    -------
    Return DataFrame ordered by ks statistic storing following values:
    - **group**: group names.
    - **names**: gene names or gene ids.
    - **statistics**: Kolmogorov-Smirnov test statistic.
    - **locations**: value from empirical cumulative distribution functions corresponding
        with the KS test statistic.
    - **signs**: +1 if the empirical distribution function of the subsample
        exceeds the empirical distribution function of the reference sample
        at 'location', otherwise -1.
    - **pvals**: two-tailed p-value.
    - **pvals_adj**: two-tailed p-value corrected by benjamini-hochberg or bonferroni.
    """

    from scipy.sparse import issparse
    from scipy.stats import kstest

    adata = adata.copy() if copy else adata

    if groups == "all":
        groups = sorted(list(adata.obs[groupby].cat.categories))
    elif isinstance(groups, Seq) and not isinstance(groups, str):
        groups = list(groups)
    else:
        raise TypeError(f"unsupported argument type for 'groups': expected sequence but received {type(groups)}")

    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        cats = sorted(list(adata.obs[groupby].cat.categories))
        raise ValueError(f"reference '{reference}' not found in adata.obs['groupby'] (avalaible values: {cats}).")

    if key_added is None:
        key_added = "smirnov_tests"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {
        "groupby": groupby,
        "reference": reference,
        "layer": layer,
        "corr_method": corr_method
    }

    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if issparse(X):
        X.eliminate_zeros()

    df = pd.DataFrame(
        columns=["group", "names", "statistics", "locations", "signs", "pvals"]
    )
    index = 0

    for name in adata.var_names:
        if reference != "rest":
            ref_sample = np.asarray(X[adata.obs[groupby] == reference, adata.var.index == name].to_numpy()).reshape(-1)
        for group in groups:
            if reference == "rest":
                ref_sample = np.asarray(X[adata.obs[groupby] != group, adata.var.index == name]).reshape(-1)
            sample = np.asarray(X[adata.obs[groupby] == group, adata.var.index == name]).reshape(-1)
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
    
    df.sort_values(
        by=["statistics"],
        ascending=False,
        inplace=True,
        ignore_index=True
    )

    if corr_method == "benjamini-hochberg":
        from statsmodels.stats.multitest import multipletests
        _, df["pvals_adj"], _, _ = multipletests(
            df["pvals"], alpha=0.05, method="fdr_bh"
        )
    elif corr_method == "bonferroni":
        n_genes = adata.var.shape[0]
        df["pvals_adj"] = np.minimum(df["pvals"] * n_genes, 1.0)

    if pval_cutoff is not None:
        df = df[df["pvals_adj"] < pval_cutoff]

    if copy:
        return df
    else:
        adata.uns[key_added]["results"] = df
    