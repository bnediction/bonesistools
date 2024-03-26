#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import os, argparse
import pickle
from pathlib import Path
from utils.argtype import Store_prefix, Range
from utils.stdout import disable_print

import anndata as ad, anndatatools as adt, stream as st

from numpy import nan
from scipy.sparse import issparse
from pandas.api.types import is_float_dtype
from networkx.classes.graph import Graph
from rpy2.rinterface import ListSexpVector

import matplotlib.pyplot as plt
from anndatatools.plotting import color

parser = argparse.ArgumentParser(
    prog="pseudotime computation of sc-RNAseq data",
    description="""From sc-RNAseq data recorded in hdf5 format, compute pseudotime based on STREAM method \
    (see Chen et al. (2019): <https://www.nature.com/articles/s41467-019-09670-4>).""",
    usage=""""python pseudotime.py [-h] -i <path> [<args>]"""
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
    "-o", "--outpath",
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=Path("./pseudotime").resolve(),
    metavar="PATH",
    help="output path (default: ./pseudotime)"
)

parser.add_argument(
    "-p", "--prefix",
    dest="prefix",
    action=Store_prefix,
    required=False,
    default="",
    metavar="LITERAL",
    help="prefix for each saving file"
)

parser.add_argument(
    "--st", "--save-tables",
    dest="save_tables",
    required=False,
    action="store_true",
    help="save the anndata object"
)

parser.add_argument(
    "-e", "--extension",
    dest="extension",
    type=str,
    required=False,
    choices=["h5ad", "pkl", "both"],
    default="h5ad",
    metavar="[h5ad | pkl | both]",
    help="""file extension of saved anndata object (default: h5ad).
    if `both` is specified, save both anndata object in h5ad and pkl format."""
)

parser.add_argument(
    "-s", "--use-stream-embedding",
    dest="use_stream_embedding",
    required=False,
    action="store_true",
    help="""compute embedding components using classic stream data preprocessing.
    Trajectories are then inferred using the output array given by stream embedding.
    Cannot be used with --use-stream-embedding argument.
    If neither --use-stream-embedding argument nor --obsm argument is used,
    searches pre-computed embedding components for inferring trajectories:
    .obsm[`X_umap`] and .obsm[`X_scanorama`]."""
)

parser.add_argument(
    "-m", "--method",
    dest="method",
    type=str,
    required=False,
    choices=["se", "mlle", "umap", "pca"],
    metavar="[se | mlle | umap | pca]",
    help="method used for dimension reduction (used only if --use-stream-embedding)"
)

parser.add_argument(
    "--obsm",
    dest="obsm",
    required=False,
    default=None,
    metavar="LITERAL",
    help="""ndarray name stored in .obsm[`obsm`] used for trajectory inference.
    Cannot be used with --use-stream-embedding argument.
    If neither --obsm argument nor --use-stream-embedding argument is used,
    searches pre-computed embedding components for inferring trajectories:
    .obsm[`X_umap`] and .obsm[`X_scanorama`]."""
)

parser.add_argument(
    "-l", "--layer",
    dest="layer",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="layer used for dimension reduction (used only if --use-stream-embedding)"
)

parser.add_argument(
    "--hvg",
    dest="hvg",
    required=False,
    action="store_true",
    help="select the most variable genes for dimension reduction (used only if --use-stream-embedding)"
)

parser.add_argument(
    "-d", "--dimensions",
    dest="n_dimensions",
    type=int,
    required=False,
    default=3,
    metavar="INT",
    help="number of components to keep (used only if --use-stream-embedding, default: 3)"
)

parser.add_argument(
    "-j", "--jobs",
    dest="n_jobs",
    type=int,
    required=False,
    default=1,
    metavar="INT",
    help="number of parallel jobs to run"
)

parser.add_argument(
    "-n", "--cluster-number",
    dest="n_clusters",
    type=int,
    required=False,
    default=5,
    metavar="INT",
    help="number of clusters to compute for elastic principal graph (default: 5)"
)

parser.add_argument(
    "-g", "--groups",
    dest="groups",
    type=str,
    required=False,
    nargs="+",
    metavar="LITERAL",
    help="clusters retrieving from adata.obs[`cluster`] used for cluster-related trajectory plotting"
)

parser.add_argument(
    "--lambda", "--epg-lambda",
    dest="epg_lambda",
    type=float,
    required=False,
    default=0.05,
    metavar="FLOAT",
    help="lambda parameter used to compute the elastic energy (default: 0.05)"
)

parser.add_argument(
    "--mu", "--epg-mu",
    dest="epg_mu",
    type=float,
    required=False,
    default=0.05,
    metavar="FLOAT",
    help="mu parameter used to compute the elastic energy (default: 0.05)"
)

parser.add_argument(
    "--alpha", "--epg-alpha",
    dest="epg_alpha",
    type=float,
    required=False,
    default=0.01,
    metavar="FLOAT",
    help="alpha parameter of the penalized elastic energy (default: 0.01)"
)

parser.add_argument(
    "--extend-leaf-nodes",
    dest="extend_leaf_nodes",
    required=False,
    action="store_true",
    help="attach new node to each leaf node in order to connect the border of data point cloud to nodes"
)

parser.add_argument(
    "--extend-mode",
    dest="extend_mode",
    type=str,
    required=False,
    choices=["QuantDists","QuantCentroid","WeigthedCentroid"],
    default="QuantDists",
    metavar="[QuantDists | QuantCentroid | WeigthedCentroid]",
    help="mode used to extend the leaves (used only if --extend-leaf-nodes, default: QuantDists)"
)

parser.add_argument(
    "--extend-parameter",
    dest="extend_parameter",
    type=float,
    action=Range,
    min=0,
    max=1,
    required=False,
    default=0.5,
    help="parameter value used to extend the leaves (used only if --extend-leaf-nodes, default: 0.5)"
)

parser.add_argument(
    "--prune-graph",
    dest="prune_graph",
    required=False,
    action="store_true",
    help="Prune the learnt elastic principal graph by filtering out trivial branches"
)

parser.add_argument(
    "--collapse-mode",
    dest="collapse_mode",
    type=str,
    required=False,
    choices=["PointNumber", "PointNumber_Extrema", "PointNumber_Leaves", "EdgesNumber", "EdgesLength"],
    default="PointNumber",
    metavar="[PointNumber | PointNumber_Extrema | PointNumber_Leaves | EdgesNumber | EdgesLength]",
    help="mode used to prune the graph (used only if --prune-graph, default: PointNumber)"
)

parser.add_argument(
    "--collapse-parameter",
    dest="collapse_parameter",
    type=float,
    required=False,
    default=5,
    metavar="FLOAT",
    help="parameter value used to prune the graph (used only if --prune-graph, default: 5)"
)

parser.add_argument(
    "--add-legend",
    dest="legend",
    required=False,
    action="store_true",
    help="add legend to figures"
)

parser.add_argument(
    "--add-graph",
    dest="graph",
    required=False,
    action="store_true",
    help="add elastic principal graph to figures"
)

parser.add_argument(
    "--plot-3d",
    dest="plot_3d",
    required=False,
    action="store_true",
    help="plot figures in three dimensions"
)

args = parser.parse_args()

data_outpath = Path(f"{args.outpath}/tables")
fig_outpath = Path(f"{args.outpath}/figures")

if not data_outpath.exists():
    os.makedirs(data_outpath)
if not fig_outpath.exists():
    os.makedirs(fig_outpath)

groups = set(args.groups)

print("Loading data...")

adata = ad.read_h5ad(args.infile)
adata.obs_names_make_unique()
adata.uns["workdir"] = args.outpath

if args.use_stream_embedding is True and args.obsm is not None:
    raise argparse.ArgumentError("--use-stream-embedding and --obsm arguments cannot be used simultaneously.")
elif args.use_stream_embedding is True:
    print("Computing embedding components using stream...")
    with disable_print():
        adata.X = adata.layers[args.layer].toarray() if issparse(adata.layers[args.layer]) else adata.layers[args.layer]
        if args.hvg:
            st.select_variable_genes(
                adata,
                loess_frac=0.02
            )
        st.select_top_principal_components(
            adata,
            first_pc=True,
            n_pc=15,
            feature="var_genes" if args.hvg else None
        )
        st.dimension_reduction(
            adata,
            method=args.method,
            feature="top_pcs",
            n_components=args.n_dimensions,
            n_neighbors=50,
            n_jobs=args.jobs
        )
        adata.uns["dr"] = f"{args.method}"
else:
    if args.obsm:
        adata.uns["dr"] = args.obsm
    if "X_umap" in adata.obsm.keys():
        adata.uns["dr"] = "X_umap"
    elif "X_scanorama" in adata.obsm.keys():
        adata.uns["dr"] = "X_scanorama"
    else:
        raise ValueError("Integrated components (`X_umap` or `X_scanorama`) in adata.obsm not found.")
    adata.obsm["X_dr"] = adata.obsm[adata.uns["dr"]].copy()

for group in groups:
    try:
        adata.obs[group] = adata.obs[group].astype(object)
    except:
        pass

print("Computing elastic principal graph...")

with disable_print():
    st.seed_elastic_principal_graph(
        adata,
        clustering="kmeans",
        n_clusters=args.n_clusters
    )
    st.elastic_principal_graph(
        adata,
        epg_alpha=args.epg_alpha,
        epg_mu=args.epg_mu,
        epg_lambda=args.epg_lambda
    )
    if args.extend_leaf_nodes is True:
        st.extend_elastic_principal_graph(
            adata,
            epg_ext_mode=args.extend_mode,
            epg_ext_par=args.extend_parameter
        )
    if args.prune_graph is True:
        st.prune_elastic_principal_graph(
            adata,
            epg_collapse_mode = args.collapse_mode,
            epg_collapse_par = args.collapse_parameter,
            epg_n_processes=args.n_jobs
        )

adata.obs["node_clusters"] = nan
adata.obs["node_clusters"] = adata.obs["node_clusters"].astype(str)

nodes_mapping = dict()
for key, value in adata.uns["flat_tree"]._node.items():
    nodes_mapping[key] = value["label"]

node_clusters = dict()
for node in nodes_mapping.keys():
    _true = adata.obs["node"] == node
    adata.obs["node_clusters"][_true] = str(nodes_mapping[node])

groups = groups.union({"kmeans", "node_clusters"})

print("Plotting trajectories...")

for _group in groups:
    fig, ax = adt.pl.embedding_plot(
        adata,
        obs=_group,
        obsm=adata.uns["dr"],
        colors=[color.blue]*(len(nodes_mapping)) + [color.lightgray] if _group == "node_clusters" else None,
        xlabel=r"$\mathrm{UMAP_{1}}$" if adata.uns["dr"] == "X_umap" else r"$\mathrm{x_{1}^{\mathrm{scanorama}}}$",
        ylabel=r"$\mathrm{UMAP_{2}}$" if adata.uns["dr"] == "X_umap" else r"$\mathrm{x_{2}^{\mathrm{scanorama}}}$",
        zlabel=r"$\mathrm{UMAP_{3}}$" if adata.uns["dr"] == "X_umap" else r"$\mathrm{x_{3}^{\mathrm{scanorama}}}$",
        add_graph=args.graph,
        add_text=True if not is_float_dtype(adata.obs[_group]) else False,
        add_legend=args.legend if _group != "node_clusters" else False,
        figwidth=6 if args.legend else 5,
        s=2,
        alpha=0.4 if not is_float_dtype(adata.obs[_group]) else 0.7,
        lgd_params={
            "title":"clusters" if _group != "condition" else "conditions",
            "labels":[string.replace("cluster ","") for string in sorted(adata.obs[_group].unique())],
            "ncol":1,
            "markerscale":5,
            "frameon":True,
            "shadow":False
        } if not is_float_dtype(adata.obs[_group]) else None,
        text={
            "fontsize":14,
            "fontweight":"extra bold"
        },
        n_components = 3 if args.plot_3d is True else 2,
        background_visible=False
    )
    adt.pl.set_default(ax)
    plt.savefig(f"{fig_outpath}/{args.prefix}{_group}_{adata.uns['dr'].split('_')[-1].lower()}_trajectory_plot")
    if args.plot_3d is True:
        pickle.dump(fig, open(Path(f"{fig_outpath}/{args.prefix}{_group}_{adata.uns['dr'].split('_')[-1].lower()}_trajectory_plot.pkl"), "wb"))
    else:
        pass

if args.save_tables:

    print("Saving data...")

    if args.extension == "h5ad" or args.extension == "both":
        for key in list(adata.obs.keys()):
            if isinstance (adata.obs[key][0], tuple):
                del adata.obs[key]
        for key in list(adata.uns.keys()):
            if isinstance(adata.uns[key], (tuple, Path, Graph, ListSexpVector)):
                del adata.uns[key]
            if key.startswith("stream_S"):
                del adata.uns[key]
        adata.write_h5ad(filename=f"{data_outpath}/{args.prefix}stream.h5ad", compression="gzip")
    elif args.extension == "pkl" or args.extension == "both":
        st.write(adata, file_name=f"{data_outpath}/{args.prefix}stream.h5ad.pkl")
