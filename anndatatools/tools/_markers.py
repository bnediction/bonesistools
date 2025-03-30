#!/usr/bin/env python

from typing import Optional, Sequence, Any
from .._typing import adata_checker

import pandas as pd
from anndata import AnnData
from numpy import log2

import scipy

from ._conversion import anndata_to_dataframe

@adata_checker
def extract_rank_genes_groups(
    adata: AnnData,
    logfc_keeping: Optional[bool] = None
) -> pd.DataFrame:
    """Extracts information in adata.uns['rank_genes_groups'] in a comprehensible way.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    keep_logfoldchanges
        Specify if dataframe columns contain log2_fold_changes computed with Scanpy.
        Since these values are inconsistent (<https://www.biostars.org/p/453129/>),
        one does prefer recompute consistent log2_fold_changes.
    
    Returns
    -------
    Dataframe with information related to gene rankings. Need to use `scanpy.tl.rank_genes_groups` function
    on anndata object before.
    """

    if "rank_genes_groups" in adata.uns.keys():
        markers_uns = adata.uns["rank_genes_groups"]
    else:
        raise ValueError("adata.uns does not contain key 'rank_genes_groups'.\
            Please use `scanpy.tl.rank_genes_groups` function before, aborting")
    
    groupby = markers_uns["params"]["groupby"]

    markers_d = {key: list() for key in ["genes", "clusters", "pvals", "adj_pvals", "scores", "log_fc"]}

    for cluster in sorted(adata.obs[groupby].unique()):
        markers_d["genes"].extend(markers_uns["names"][cluster])
        markers_d["clusters"].extend([cluster] * adata.n_vars)
        markers_d["pvals"].extend(markers_uns["pvals"][cluster])
        markers_d["adj_pvals"].extend(markers_uns["pvals_adj"][cluster])
        markers_d["scores"].extend(markers_uns["scores"][cluster])
        if logfc_keeping is True:
            markers_d["log_fc"].extend(markers_uns["logfoldchanges"][cluster])
        else:
            markers_d["log_fc"].extend([float("nan")] * adata.n_vars)

    return pd.DataFrame.from_dict(markers_d, orient="columns")

@adata_checker
def log_fold_changes(
    adata: AnnData,
    groupby: str,
    layer: Optional[str] = None,
    is_log: Optional[bool] =False,
    cluster_rebalancing: Optional[bool] = False
) -> pd.DataFrame:
    """Log2 fold change is a metric translating how much the transcript's expression
    has changed between cells in and out of a cluster. The reported values are based
    on a logarithmic scale to base 2 with respect to the fold change ratios.
    According to <https://www.biostars.org/p/453129/>, computed log2 fold changes
    are different between FindAllMarkers (package Seurat) and rank_gene_groups
    (module Scanpy) functions. As mentionned, results derived from Scanpy are inconsistent.
    This current function `log_fold_changes` computes it in the right way, with identical
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

        __df = pd.DataFrame(log2(mean_in) - log2(mean_out), columns=["log_fc"])
        __df.reset_index(names="genes", inplace=True)
        __df.insert(0, "clusters", cluster)
        return __df
    
    logfc_df = pd.DataFrame(columns=["clusters","genes","log_fc"])
    counts_df = anndata_to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log)

    if cluster_rebalancing:
        mean_counts_df = counts_df.groupby(by=groupby, sort=True).mean()
        for cluster in sorted(pd.unique(adata.obs[groupby])):
            _mean_in = mean_counts_df.loc[cluster]
            _mean_out = mean_counts_df.drop(index=cluster, inplace=False).mean()
            _logfc_df = compute_logfc(_mean_in, _mean_out, cluster)
            logfc_df = pd.concat([logfc_df, _logfc_df.copy()])
    else:
        for cluster in sorted(pd.unique(adata.obs[groupby])):
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

    logfc_df = log_fold_changes(
        adata,
        groupby=groupby,
        layer=layer,
        is_log=is_log,
        cluster_rebalancing=cluster_rebalancing
    )
    df = df.loc[:, df.columns != "log_fc"]
    if threshold:
        logfc_df = logfc_df.loc[logfc_df["log_fc"] > threshold]
    df = pd.merge(
        df,
        logfc_df,
        left_on=["genes", "clusters"],
        right_on=["genes", "clusters"],
        how="inner"
    )
    return df

def hypergeometric_test(
    adata: AnnData,
    signature: Sequence[str],
    markers: Sequence[str],
) -> float:
    """Computes the p-value (or survival function) of an hypergeometric
    distribution using scRNA-seq data in order to test whether marker genes
    significantly match signature genes.
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
        A signature is a set of overexpressed genes in a cell-type.
    markers
        Set of markers (genes) in a given cluster.
        A marker set is a set of overexpressed genes in a cluster.
    """

    if not isinstance(adata, AnnData):
        raise TypeError(f"Argument `adata` must be of type {type(AnnData)}, not {type(adata)}")
    
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
    
    return scipy.stats.hypergeom.sf(k = k, M = N, n = K, N = n, loc = 1)

@adata_checker
def multiple_hypergeometric_test(
    adata: AnnData,
    signatures: dict,
    markers: pd.DataFrame,
    cluster: str,
) -> dict:

    _markers = markers[markers["clusters"] == cluster]["genes"]
    return {cell_type: hypergeometric_test(adata, signature, _markers) for cell_type, signature in signatures.items()}

@adata_checker
def get_info(
    adata: AnnData,
    signatures: dict,
    markers: pd.DataFrame,
    groupby: str = "cluster",
    by: Optional[Any] = None,
) -> dict:

    columns = ["genes", "clusters", "pvals", "adj_pvals", "scores", "log_fc"]
    for idx, column in enumerate(columns):
        if not column == markers.columns[idx]:
            raise ValueError("`markers` dataframe must contain specific rows with the specific order\
                `genes`, `clusters`, `pvals`, `adj_pvals`, `scores`, `log_fc`")

    if by:
        group_ad = adata[adata.obs[groupby] == by]
        group_info_d = dict()
        group_info_d["n_cells"] = group_ad.n_obs
        group_info_d["proportion_cells"] = round(group_ad.n_obs / adata.n_obs, ndigits=6)
        proportion_phases = group_ad.obs["pypairs_max_class"].value_counts() / group_ad.n_obs
        group_info_d.update({phase: round(proportion_phases[phase], ndigits=6) for phase in sorted(proportion_phases.index)})
        group_info_d["median_expressed_genes"] = group_ad.obs["n_genes_by_counts"].median()
        group_info_d["median_total_counts"] = group_ad.obs["total_counts"].median()
        group_info_d["median_proportion_mito"] = f"{group_ad.obs['pct_counts_mitochondrion'].median():.4f}%"
        pvalues_d = multiple_hypergeometric_test(group_ad, signatures, markers, cluster=by)
        group_info_d.update({cell_type: round(pvalue, ndigits=6) for cell_type, pvalue in pvalues_d.items()})
        return group_info_d
    else:
        return {group: get_info(adata, signatures, markers, groupby=groupby, by=group) for group in sorted(adata.obs[groupby].unique())}
