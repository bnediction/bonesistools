#!/usr/bin/python3

import warnings
warnings.filterwarnings("ignore")

import random
random.seed(100)

import os, argparse
from pathlib import Path
from utils.argtype import Store_prefix

import numpy as np, math

import pandas as pd, anndata as ad, scanpy as sc, json
import anndatatools as adt

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from anndatatools.plotting import (
    fig,
    color,
    color_cycle
)

parser = argparse.ArgumentParser(
    prog="Clusterization of sc-RNAseq data",
    description="""From one-condition sc-rnaSeq data recorded in the hdf5 format (<filename>.h5ad),
    perform dimension reduction using PCA, create clusters using leiden algorithm,
    run t-SNE and UMAP algorithm, search for gene markers and compare markers and
    signatures in order to provide useful information about potential cell-types
    of each cluster.""",
    usage="python cluster.py -i <path> -s <path> [<args>]"
)

parser.add_argument(
    "-i", "--infile",
    dest="infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad file (including file)"
)

parser.add_argument(
    "-s", "--signatures",
    dest="signatures",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .json signatures file (including file)"
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
    "-p", "--prefix",
    dest="prefix",
    action=Store_prefix,
    required=False,
    default="",
    help="prefix for each saving file"
)

parser.add_argument(
    "-k", "--k-neighbors",
    dest="k_neighbors",
    type=int,
    required=False,
    default=20,
    metavar="PATH",
    help="number of closest neighbors computed when computing KNN graph (default: 20)"
)

parser.add_argument(
    "-ng", "--neighborhood-graph",
    dest="neighborhood_graph",
    type=str,
    required=False,
    default="knn",
    choices=["knn","snn"],
    metavar="[knn | snn]",
    help="neighborhood graph used by Leiden clustering algorithm (default: knn)"
)

parser.add_argument(
    "-n", "--dimensions",
    dest="n_dimensions",
    type=int,
    required=False,
    default=15,
    metavar="INT",
    help="number of principal components taken into account for clustering cells and running t-SNE/UMAP (default: 15)"
)

parser.add_argument(
    "-r", "--resolution",
    dest="resolution",
    type=float,
    required=False,
    default=0.6,
    metavar="FLOAT",
    help="parameter value controlling the coarseness of the clustering when using Leiden algorithm (default: 0.6)"
)

parser.add_argument(
    "-l", "--logfc-threshold",
    dest="logfc_threshold",
    type=float,
    required=False,
    default=0.25,
    metavar="FLOAT",
    help="threshold describing the minimum log2 fold-changes for being a gene marker (default: 0.25)"
)

parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    required=False,
    default=False,
    action="store_true",
    help="get summarizing information about cluster in stdout"
)

args = parser.parse_args()

data_outpath = Path(f"{args.outpath}/tables")
fig_outpath = Path(f"{args.outpath}/figures")

if not data_outpath.exists():
    os.makedirs(data_outpath)
if not fig_outpath.exists():
    os.makedirs(fig_outpath)

print(f"Loading data...")

adata = sc.read_h5ad(args.infile)

n_comps = 50 if args.n_dimensions <= 15 else args.n_dimensions

resolutions = [0.6,0.8,1,1.2]

color_d = {
    "G1": color.blue,
    "G2M": color.red,
    "S": color.green
}
phase = adata.obs["pypairs_cc_prediction"]

print(f"Running principal component analysis (PCA)...")

adata.X = adata.layers["correct"]
sc.tl.pca(
    adata,
    zero_center=False,
    n_comps=n_comps,
    use_highly_variable=True,
    copy=False
)

pc1 = adata.obsm["X_pca"][:,0]
pc2 = adata.obsm["X_pca"][:,1]
fig, ax = plt.subplots(nrows=1, ncols=1)
for p in np.unique(phase):
    idx = np.where(phase == p)[0]
    ax.scatter(pc1[idx], pc2[idx], s=5, facecolors=color_d[p], edgecolors="none", alpha=1, label=p)
ax.set_xlabel(r"$\mathrm{PC_{1}}$")
ax.set_ylabel(r"$\mathrm{PC_{2}}$")
ax.legend(markerscale=2, edgecolor=color.black)
plt.sca(ax)
ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
plt.savefig(f"{fig_outpath}/{args.prefix}principal-component-analysis")

print(f"Clustering...")

knn_key = "knn"
snn_key = "snn"
sc.pp.neighbors(adata, n_neighbors=args.k_neighbors, n_pcs=args.n_dimensions, key_added=knn_key, copy=False)
adt.tl.shared_neighbors(adata, knn_key=knn_key, snn_key=snn_key, prune_snn = 1/15, copy=False)

for resolution in resolutions:
    sc.tl.leiden(adata, resolution=resolution, neighbors_key=knn_key, key_added=f"leiden_{resolution}")

if args.resolution in resolutions and args.neighborhood_graph == "knn":
    adata.obs["leiden"] = adata.obs[f"leiden_{args.resolution}"]
elif args.neighborhood_graph == "knn":
    sc.tl.leiden(adata, resolution=args.resolution, neighbors_key=knn_key, key_added=f"leiden")
else:
    obsp = adata.uns[snn_key]["similarities_key"]
    sc.tl.leiden(adata, resolution=resolution, adjacency=adata.obsp[obsp].copy(), key_added=f"leiden")

print(f"Running t-SNE...")

sc.tl.tsne(adata, n_pcs=args.n_dimensions, learning_rate=1000)

tsne1 = adata.obsm["X_tsne"][:,0]
tsne2 = adata.obsm["X_tsne"][:,1]

fig, axes = plt.subplots(nrows=2, ncols=2)
fig.set_figheight(8)
fig.set_figwidth(8)
for i, resolution in enumerate(resolutions):
    for _cluster, _color in zip(np.unique(adata.obs[f"leiden_{resolution}"]), color_cycle):
        idx = np.where(adata.obs[f"leiden_{resolution}"] == _cluster)[0]
        axv, axh = [math.floor(i/2), i%2]
        axes[axv, axh].scatter(tsne1[idx], tsne2[idx], s=2, facecolors=_color, edgecolors="none", alpha=1, label=_cluster)
        axes[axv, axh].title.set_text(f"resolution: {resolution}")
        if axv == 1:
            axes[axv, axh].set_xlabel(r"$t$-$\mathrm{SNE_{1}}$")
        if axh == 0:
            axes[axv, axh].set_ylabel(r"$t$-$\mathrm{SNE_{2}}$")
plt.savefig(f"{fig_outpath}/{args.prefix}tsne_clusters")

print(f"Running uniform manifold approximation and projection (UMAP)...")

sc.tl.umap(adata, neighbors_key=knn_key, n_components=2)

umap1 = adata.obsm["X_umap"][:,0]
umap2 = adata.obsm["X_umap"][:,1]

fig, axes = plt.subplots(nrows=2, ncols=2)
fig.set_figheight(8)
fig.set_figwidth(8)
for i, resolution in enumerate(resolutions):
    for _cluster, _color in zip(np.unique(adata.obs[f"leiden_{resolution}"]), color_cycle):
        idx = np.where(adata.obs[f"leiden_{resolution}"] == _cluster)[0]
        axv, axh = [math.floor(i/2), i%2]
        axes[axv, axh].scatter(umap1[idx], umap2[idx], s=2, facecolors=_color, edgecolors="none", alpha=1, label=_cluster)
        axes[axv, axh].title.set_text(f"resolution: {resolution}")
        if axv == 1:
            axes[axv, axh].set_xlabel(r"$\mathrm{UMAP_{1}}$")
        if axh == 0:
            axes[axv, axh].set_ylabel(r"$\mathrm{UMAP_{2}}$")
plt.savefig(f"{fig_outpath}/{args.prefix}umap_clusters")

fig, ax = plt.subplots(nrows=1, ncols=1)
for p in np.unique(phase):
    idx = np.where(phase == p)[0]
    ax.scatter(umap1[idx], umap2[idx], s=2, facecolors=color_d[p], edgecolors="none", alpha=1, label=p)
ax.set_xlabel(r"$\mathrm{UMAP_{1}}$")
ax.set_ylabel(r"$\mathrm{UMAP_{2}}$")
ax.legend(markerscale=5, edgecolor=color.black)
plt.sca(ax)
ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
plt.savefig(f"{fig_outpath}/{args.prefix}umap_phases")

for metric in ["total_counts", "pct_counts_mitochondrion"]:
    fig, ax = plt.subplots(nrows=1, ncols=1)
    if metric == "total_counts":
        cmap = "Greens"
        label = r"$\# \mathrm{read\ counts}$"
    elif metric == "pct_counts_mitochondrion":
        cmap = "Blues"
        label = r"$\frac{\# \mathrm{mitochondrion\ counts}}{\# \mathrm{read\ counts}}$"
    mapping = ax.scatter(umap1, umap2, s=2, c=adata.obs[metric], cmap=cmap, alpha=1)
    cbar = fig.colorbar(mapping)
    cbar.set_label(label, loc="center", labelpad=5)
    ax.set_xlabel(r"$\mathrm{UMAP_{1}}$")
    ax.set_ylabel(r"$\mathrm{UMAP_{2}}$")
    plt.sca(ax)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    plt.savefig(f"{fig_outpath}/{args.prefix}umap_{metric}")

print(f"Marker analysis...")

layer = "log-normalize"
groupby = "leiden"

sc.tl.rank_genes_groups(
    adata,
    layer=layer,
    use_raw=False,
    groupby=groupby,
    reference="rest",
    method="wilcoxon",
    tie_correct=True,
    corr_method="bonferroni"
)
markers_df = adt.tl.extract_rank_genes_groups(
    adata,
    logfc_keeping=False
)
markers_df = markers_df.loc[markers_df["adj_pvals"] < 0.05]
markers_df = adt.tl.update_logfoldchanges(
    markers_df,
    adata,
    groupby,
    layer,
    threshold=args.logfc_threshold
)

print(f"Signature analysis...")

with open(args.signatures, "r") as signatures_f:
    signatures_d = json.load(signatures_f)

valid_gene_names = set(adata.var_names)
for cell_type, signature in signatures_d.items():
    signatures_d[cell_type] = {gene for gene in signature if gene in valid_gene_names}
signatures_d = {cell_type: signature for cell_type, signature in signatures_d.items() if signature}
del valid_gene_names

layer="log-normalize"
adata.X = adata.layers[layer]
for cell_type, signature in signatures_d.items():
    sc.tl.score_genes(
        adata,
        gene_list=signature,
        gene_pool=None,
        n_bins=25,
        ctrl_size=100,
        score_name=cell_type,
        random_state=0,
        copy=False,
        use_raw=False
    )

print("Summarizing clusters...")

clust_info_d = adt.tl.get_info(adata, signatures_d, markers_df, groupby=groupby)
clust_info_df = pd.DataFrame.from_dict(clust_info_d, orient="index")

print("Saving data...")

adata.write_h5ad(filename=f"{data_outpath}/{args.prefix}counts.h5ad", compression="gzip")
markers_df.to_csv(f"{data_outpath}/{args.prefix}markers.csv", sep=",", index=False)
clust_info_df.to_csv(f"{data_outpath}/{args.prefix}cluster_info.csv", sep=",", index=True)

if args.verbose:
    print(clust_info_df.transpose())
