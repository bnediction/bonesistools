#!/usr/bin/env python

from typing import Optional, Sequence
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
    is_log: Optional[bool] =False,
    cluster_rebalancing: Optional[bool] = False
) -> pd.DataFrame:
    """Log2 fold-change is a metric translating how much the transcript's expression
    has changed between cells in and out of a cluster. The reported values are based
    on a logarithmic scale to base 2 with respect to the fold change ratios.
    According to <https://www.biostars.org/p/453129/>, computed log2 fold changes
    are different between FindAllMarkers (package Seurat) and rank_gene_groups
    (module Scanpy) functions. As mentionned, results derived from Scanpy are inconsistent.
    This current function 'calculate_logfoldchanges' computes it in the right way, with identical
    results to Seurat by keeping default options.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        Any key in anndata.obs corresponding to defined clusters or groups.
    layer
        Any key in anndata.layer.
        If not specify, log2 fold changes are derived from anndata.X.
    is_log
        Boolean value specifying if the counts is logarithmized.
    cluster_rebalancing
        If no cluster rebalancing, cells are equally-weighted.
        Otherwise, cells are weighted with cluster size such as clusters are equally-weighted.
        It means that cells in small cluster have a greater weight than other cells in order
        to correct cluster size effects.    
    """

    def compute_logfc(mean_in, mean_out, cluster):

        __df = pd.DataFrame(log2(mean_in) - log2(mean_out), columns=["logfoldchanges"])
        __df.reset_index(names="names", inplace=True)
        __df.insert(0, "group", cluster)
        return __df
    
    logfc_df = pd.DataFrame(columns=["group","names","logfoldchanges"])
    counts_df = anndata_to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log)

    if cluster_rebalancing:
        mean_counts_df = counts_df.groupby(by=groupby, sort=True).mean()
        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = mean_counts_df.loc[cluster]
            _mean_out = mean_counts_df.drop(index=cluster, inplace=False).mean()
            _logfc_df = compute_logfc(_mean_in, _mean_out, cluster)
            logfc_df = pd.concat([logfc_df, _logfc_df.copy()])
    else:
        for cluster in sorted(adata.obs[groupby].unique().dropna()):
            _mean_in = counts_df.loc[counts_df[groupby] == cluster, counts_df.columns != groupby].mean()
            _mean_out = counts_df.loc[counts_df[groupby] != cluster, counts_df.columns != groupby].mean()
            _logfc_df = compute_logfc(_mean_in, _mean_out, cluster)
            logfc_df = pd.concat([logfc_df, _logfc_df.copy()])
            del _logfc_df

    return logfc_df.reset_index(drop=True)

def update_logfoldchanges(
    df: pd.DataFrame,
    adata: AnnData,
    groupby: str,
    layer: str,
    is_log: Optional[bool] = True,
    cluster_rebalancing: Optional[bool] = False,
    threshold: Optional[float] = None
) -> pd.DataFrame:

    logfc_df = calculate_logfoldchanges(
        adata,
        groupby=groupby,
        layer=layer,
        is_log=is_log,
        cluster_rebalancing=cluster_rebalancing
    )
    df = df.loc[:, df.columns != "logfoldchanges"]
    if threshold:
        logfc_df = logfc_df.loc[logfc_df["logfoldchanges"] > threshold]
    df = pd.merge(
        df,
        logfc_df,
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
    """Estimates the p-value (or survival function) of an hypergeometric
    distribution in order to test whether marker genes significantly
    match signature genes.
    Given a population size N and a number of success states K,
    it describes the probability of having at least k successes
    in n draws, without replacement, where:
    - N is the number of genes in anndata,
    - K is the number of signature genes,
    - n is the number of markers,
    - k is the number of gene matching both signature genes and markers.
    Smaller the p-value, higher the probability that genes of the given
    cluster comes from the cell-type associated to the given signature.    

    Parameters
    ----------
    adata
        Annotated data matrix.
    signature
        Set of signature genes in a given cell-type.
        A signature is a set of over-expressed genes in a cell-type.
    markers
        Set of markers (genes) in a given cluster.
        A marker set is a set of over-expressed genes in a cluster.
    """

    if not isinstance(adata, AnnData):
        raise TypeError(f"Argument 'adata' must be of type {type(AnnData)}, not {type(adata)}")
    
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
