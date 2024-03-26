#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import os, argparse
from pathlib import Path

import numpy as np, random, math

import pandas as pd, scanpy as sc, rdata
from pypairs import pairs

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from anndatatools.plotting import (
    fig,
    color
)

random.seed(1000)

pd.DataFrame.iteritems = pd.DataFrame.items

def median_absolute_deviation(x, consistency=False):
    """Compute the mean absolute deviation (MAD),
    i.e. the median of the absolute deviations from the median.
    If consistency is true, adjust by a factor for asymptotically normal consistency.
    Asymptotic normal consistency means that:
        E[MAD(X_1,...,X_n)] = sigma
    for X_i following a gaussian distribution N(mu, sigma^2).
    """
    constant = 1.4826 if consistency else 1
    return constant*np.median(np.absolute(x - np.median(x)))

def marker_pairs_converter(ensembl_to_symbol: dict, ensembl_marker_pairs):
    """convert marker pairs from ensembl id into symbol id"""
    symbol_marker_pairs = dict()
    for cycle, pairs in ensembl_marker_pairs.items():
        cycle_pairs = list()
        for _, (first, second) in pairs.iterrows():
            if first in ensembl_to_symbol.keys() and second in ensembl_to_symbol.keys():
                cycle_pairs.append([ensembl_to_symbol[first], ensembl_to_symbol[second]])
        symbol_marker_pairs[cycle] = cycle_pairs
    return symbol_marker_pairs

parser = argparse.ArgumentParser(
    prog="Cell filtering and cell cycle phases assignement on sc-RNAseq data",
    description="""From sc-rnaSeq data recorded in the hdf5 format (<filename>.h5ad),
    filter low-quality cells and assign cell cycle phases using marker pairs.
    Some quality metrics are computed before and after filtering.""",
    usage="python cluster.py [-h] --i <path> [<args>]"
)

parser.add_argument(
    "-i", "--infile",
    dest="count_infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad file (including file)"
)

parser.add_argument(
    "-c", "--cycle","--marker",
    dest="marker_infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .rds cell cycle phase markers file (including file)"
)

parser.add_argument(
    "-o", "--outpath",
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=Path("./").resolve(),
    metavar="PATH",
    help="output path"
)

parser.add_argument(
    "-t", "--mitochondrial_threshold",
    dest="mitochondrial_threshold",
    type=float,
    required=False,
    default=100,
    metavar="FLOAT",
    help="maximum mitochondria threshold in a cell (value between 0 and 100)"
)

parser.add_argument(
    "-u", "--upper-mad",
    dest="upper_mad",
    type=float,
    required=False,
    default=2,
    metavar="FLOAT",
    help="factor removing cells for which their count are higher than this factor*MAD with respect to the median count"
)

parser.add_argument(
    "-l", "--lower-mad",
    dest="lower_mad",
    type=float,
    required=False,
    default=2,
    metavar="FLOAT",
    help="factor removing cells for which their count are lower than this factor*MAD with respect to the median count"
)

parser.add_argument(
    "-m", "--consistency-mad",
    dest="consistency_mad",
    required=False,
    action="store_true",
    help="median absolute deviation (MAD) is refactorised for asymptotically normal consistency"
)

args = parser.parse_args()

data_outpath = Path(f"{args.outpath}/tables")
fig_outpath = Path(f"{args.outpath}/figures")

if not data_outpath.exists():
    os.makedirs(data_outpath)
if not fig_outpath.exists():
    os.makedirs(fig_outpath)

print(f"Loading data...")

adata = sc.read_h5ad(Path(f"{args.count_infile}").resolve())
_s_ante_filter = adata.shape

ensembl_to_symbol = dict()
for _, row in adata.var.iterrows():
    ensembl_to_symbol[row["gene_ids"]] = row["symbol"]

print(f"Assigning cell cycle phases...")

parser = rdata.parser.parse_file(args.marker_infile)
marker_pairs = rdata.conversion.convert(parser)
marker_pairs = marker_pairs_converter(ensembl_to_symbol, marker_pairs)
scores = pairs.cyclone(adata, marker_pairs)
adata.obs.rename(columns={
    "pypairs_G1": "G1_score",
    "pypairs_S": "S_score",
    "pypairs_G2M": "G2M_score"
}, inplace=True)

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.scatter(adata.obs.G1_score, adata.obs.G2M_score, s=30, facecolors=color.white, edgecolors=color.blue, alpha=1)
ax.set_xlabel(r"score $\mathrm{G_{1}}$")
ax.set_ylabel(r"score $\mathrm{G_{2}/M}$")
plt.sca(ax)
ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
plt.savefig(f"{fig_outpath}/cell-cycle-phases-assignment")
plt.close()

adata.var_names_make_unique()

adata.var["mitochondrion"] = adata.var_names.str.startswith("mt-")          # annotate the group of mitochondrial genes
adata.var["ribosome"] = adata.var_names.str.startswith(("Rps","Rpl","Mrp")) # annotate the group of ribosomal genes

print("Computing Violin plot before cell filtering...")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True, qc_vars=["mitochondrion","ribosome"])
ax = sc.pl.violin(
    adata=adata,
    keys=["n_genes_by_counts", "total_counts", "pct_counts_mitochondrion", "pct_counts_ribosome"],
    jitter=0.4,
    multi_panel=True,
    stripplot=False,
    show=False,
    save=False
)
ax.axes[0,0].set_title(r"gene number")
ax.axes[0,1].set_title(r"gene counts")
ax.axes[0,2].set_title(r"mitochondrion proportion")
ax.axes[0,3].set_title(r"ribosome proportion")
plt.savefig(f"{fig_outpath}/violin-plot-on-UMI-before-filtering.png")

print(f"Filtering low-quality cells...")

min_counts_threshold = np.exp(np.median(np.log(adata.obs.total_counts)) \
    - args.lower_mad*median_absolute_deviation(np.log(adata.obs.total_counts),consistency=args.consistency_mad))
max_counts_threshold = np.exp(np.median(np.log(adata.obs.total_counts)) \
    + args.upper_mad*median_absolute_deviation(np.log(adata.obs.total_counts),consistency=args.consistency_mad))

_ylim = [0, round(math.ceil(max(adata.obs.total_counts)+1000),-3)]
fig, ax = plt.subplots(nrows=1, ncols=2)
sc.pl.violin(
    adata=adata,
    keys="total_counts",
    jitter=0.4,
    multi_panel=None,
    stripplot=False,
    ax=ax[0],
    show=False,
    save=False,
)
[ax[0].axhline(threshold, linewidth=1.5, linestyle='--', color=color.red) for threshold in [min_counts_threshold, max_counts_threshold]]
ax[0].set_ylim(_ylim)
ax[0].set(title="before cell filtering")

sc.pp.filter_cells(adata, min_counts=min_counts_threshold, inplace=True)
sc.pp.filter_cells(adata, max_counts=max_counts_threshold, inplace=True)
adata = adata[adata.obs.pct_counts_mitochondrion < args.mitochondrial_threshold, :]
_s_post_filter = adata.shape

sc.pl.violin(
    adata=adata,
    keys="total_counts",
    jitter=0.4,
    multi_panel=None,
    stripplot=False,
    ax=ax[1],
    show=False,
    save=False,
)
ax[1].axhline(min_counts_threshold, linewidth=1.5, linestyle='--', color=color.red)
ax[1].axhline(max_counts_threshold, linewidth=1.5, linestyle='--', color=color.red)
ax[1].set_ylim(_ylim)
ax[1].set(title="after cell filtering")
plt.savefig(f"{fig_outpath}/violin-plot-on-barcode-counts.png")

fig, ax = plt.subplots(nrows=1, ncols=1)
ax = adata.obs.pypairs_cc_prediction.value_counts().plot.bar(rot=0)
ax.set(xlabel="cell cycle phases")
plt.savefig(f"{fig_outpath}/assigned-cell-cycle-phases-counting")

print("Computing Violin plot after cell filtering...")

ax = sc.pl.violin(
    adata=adata,
    keys=["n_genes_by_counts", "total_counts", "pct_counts_mitochondrion", "pct_counts_ribosome"],
    jitter=0.4,
    multi_panel=True,
    stripplot=False,
    show=False,
    save=False
)
ax.axes[0,0].set_title(r"gene number")
ax.axes[0,1].set_title(r"gene counts")
ax.axes[0,2].set_title(r"mitochondrion proportion")
ax.axes[0,3].set_title(r"ribosome proportion")
plt.savefig(f"{fig_outpath}/violin-plot-on-UMI-after-filtering2.png")

print("Saving data...")

adata.write_h5ad(filename=f"{data_outpath}/counts.h5ad")

print(f"Dimension of the couting matrix before filtering: {_s_ante_filter}")
print(f"Dimension of the couting matrix after filtering: {_s_post_filter}")