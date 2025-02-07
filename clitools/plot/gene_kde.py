#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Optional, Tuple

import os, argparse
from pathlib import Path

import pandas as pd
import anndata as ad
import scanpy as sc
import anndatatools as adt

import scipy

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
from anndatatools.plotting import (
    fig,
    color
)

class GeneKde(object):

    def __init__(self, adata: ad.AnnData):

        if isinstance(adata, ad.AnnData):
            self.adata = adata
        else:
            raise TypeError(f"unsupported type for instantiate '{type(self)}': '{type(adata)}' instead of '{type(ad.AnnData)}'")

    def plot(
        self,
        gene: str,
        layer: Optional[str]=None,
        obs: Optional[str]=None
    ) -> Tuple[mpl.figure.Figure, mpl.axes._axes.Axes]:

        counts = self.adata[:,gene].layers[layer] if layer else self.adata[:,gene].X
        if scipy.sparse.issparse(counts):
            counts = pd.Series(counts.toarray().squeeze(), index=self.adata.obs.index, name="counting")
        if obs:
            counts = pd.concat([counts, self.adata.obs[obs].astype("category")], axis=1)
            colors = color.COLORS[0:len(counts[obs].cat.categories)]
        
        fig, ax = plt.subplots()
        sns.kdeplot(
            data=counts["counting"],
            ax=ax,
            color=color.gray if obs else color.red,
            fill=True,
            label="all"
        )
        if obs:
            for _cluster, _color in zip(sorted(counts["clusters"].cat.categories), colors):
                _counts = counts.loc[counts["clusters"] == _cluster]["counting"]
                sns.kdeplot(
                    data=_counts,
                    ax=ax,
                    color=_color,
                    fill=False,
                    label=_cluster
                )
        if min(counts["counting"]) == 0:
            plt.xlim(min(counts["counting"]),max(counts["counting"])*1.1)
        if obs:
            ax.legend(loc='upper right')

        return fig, ax

parser = argparse.ArgumentParser(
    prog="Single-cell gene kde",
    description="""Compute gene-related kernel density distribution from single-cell sequencing data.""",
    usage="python gene_distribution.py <FILE> <PATH> [<args>]"
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="counting file (h5ad format)"
)

parser.add_argument(
    "outpath",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="output path"
)

parser.add_argument(
    "--genes",
    dest="genes",
    type=str,
    required=True,
    nargs="+",
    metavar="LITERAL",
    help="genes of interest"
)

parser.add_argument(
    "--layer",
    dest="layer",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="layer used (if not specified, consider adata.X)"
)

parser.add_argument(
    "--obs",
    dest="obs",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="plot kde for each distinct group in adata.obs (default: None)"
)

args = parser.parse_args()

if not args.outpath.exists():
    os.makedirs(args.outpath)

adata = sc.read_h5ad(args.infile)

gene_kde = GeneKde(adata)

gene = args.genes[0]
for gene in args.genes:
    fig, ax = gene_kde.plot(gene, layer=args.layer, obs=args.obs)
    plt.savefig(Path(f"{args.outpath}/{gene}.pdf"))
    plt.close()
