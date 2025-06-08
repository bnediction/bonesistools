#!/usr/bin/env python

from typing import (
    Optional,
    Sequence,
    Callable
)
from .._typing import anndata_checker

import pandas as pd
from anndata import AnnData
from numpy import log2

import scipy

from ._conversion import anndata_to_dataframe

@anndata_checker
def calculate_logfoldchanges(
    adata: AnnData,
    groupby: str,
    layer: Optional[str] = None,
    is_log: bool =False,
    cluster_rebalancing: bool = False
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
        Any key in anndata.obs corresponding to defined clusters or groups.
    layer: str (optional, default: None)
        Any key in anndata.layer.
        If not specify, log2 fold changes are derived from anndata.X.
    is_log: bool (default: False)
        Specify whether counts are logarithmized or not.
    cluster_rebalancing: bool (default: False)
        If no cluster rebalancing, cells are equally-weighted.
        Otherwise, cells are weighted with cluster size such as clusters are equally-weighted.
        It means that cells in small cluster have a greater weight than other cells in order
        to correct cluster size effects.
    
    Returns
    -------
    Return DataFrame.
    DataFrame stores computed log2 fold-changes.
    
    See also
    --------
    Get more information about the difference between log2 fold-changes derived with Seurat and Scanpy here:
    <https://www.biostars.org/p/453129/>
    """

    def compute_logfc(mean_in, mean_out, cluster):

        __df = pd.DataFrame(log2(mean_in) - log2(mean_out), columns=["logfoldchanges"])
        __df.reset_index(names="names", inplace=True)
        __df.insert(0, "group", cluster)
        return __df
    
    logfoldchanges_df = pd.DataFrame(columns=["group","names","logfoldchanges"])
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

    return logfoldchanges_df.reset_index(drop=True)

def update_logfoldchanges(
    df: pd.DataFrame,
    adata: AnnData,
    groupby: str,
    layer: str,
    is_log: Optional[bool] = True,
    cluster_rebalancing: Optional[bool] = False,
    filter_logfoldchanges: Optional[Callable]=None,
) -> pd.DataFrame:

    logfoldchanges_df = calculate_logfoldchanges(
        adata,
        groupby=groupby,
        layer=layer,
        is_log=is_log,
        cluster_rebalancing=cluster_rebalancing
    )
    df = df.loc[:, df.columns != "logfoldchanges"]

    if filter_logfoldchanges is not None:
        if not callable(filter_logfoldchanges):
            raise TypeError(f"unsupported argument type for 'filter_logfoldchanges': expected callable object")
        else:
            logfoldchanges_df = logfoldchanges_df.loc[filter_logfoldchanges(logfoldchanges_df["logfoldchanges"].values)]

    df = pd.merge(
        df,
        logfoldchanges_df,
        left_on=["names", "group"],
        right_on=["names", "group"],
        how="inner"
    )
    return df

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
    
    return scipy.stats.hypergeom.sf(
        k=k,
        M=N,
        n=K,
        N=n,
        loc=1
    )
