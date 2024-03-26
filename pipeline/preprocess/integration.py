#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import os, argparse
import pickle
from pathlib import Path
from utils.stdout import Section
from utils.argtype import Store_prefix

from typing import Sequence

from collections import OrderedDict as odict

import anndata as ad, anndatatools as adt
import scanpy as sc
import scanorama

import matplotlib.pyplot as plt
from anndatatools.plotting import color

@adt._adata_arg_checking
def clean_adata(
    adata: ad.AnnData,
    obs: Sequence[str] = None,
    var: Sequence[str] = None,
    copy: bool = False
    ) -> None:
        
        adata = adata.copy() if copy else adata

        if obs:
            for _obs in obs:
                if _obs in adata.obs.columns:
                    del adata.obs[_obs]
        if var:
            for _var in var:
                if _var in adata.var.columns:
                    del adata.var[_var]
        if "pca" in adata.uns.keys():
            del adata.uns["pca"]
        if "neighbors" in adata.uns.keys():
            del adata.uns["neighbors"]
        if "leiden" in adata.uns.keys():
            del adata.uns["leiden"]
        if "tsne" in adata.uns.keys():
            del adata.uns["tsne"]
        if "umap" in adata.uns.keys():
            del adata.uns["umap"]
        if "hvg" in adata.uns.keys():
            try:
                del adata.uns["hvg"]
                del adata.var["highly_variable"]
                del adata.var["highly_variable_rank"]
            except:
                pass
        del adata.obsm, adata.obsp, adata.varm, adata.varp

        return adata if copy else None

parser = argparse.ArgumentParser(
    prog="Integration and clustering on sc-RNAseq data",
    description="""From two samples of sc-rnaSeq data recorded in the hdf5 format (<filename>.h5ad),
    perform integration on embedding dimensions, create clusters using leiden algorithm,
    and run UMAP algorithm. This programm allows to search cell evolutions between two experiments""",
    usage="""python cluster.py [-h] --i1 <path> --i2 <path> [<args>]"""
)

parser.add_argument(
    "--i1", "--infile1", "--infile-ref",
    dest="infile_ref",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad file (including file) considered as reference sample"
)

parser.add_argument(
    "--i2", "--infile2", "--infile-interest",
    dest="infile_interest",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad file (including file) considered as sample to integrate"
)

parser.add_argument(
    "-o", "--outpath",
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=Path("./").resolve(),
    metavar="PATH",
    help="output path (default: ./)"
)

parser.add_argument(
    "--prefix",
    dest="prefix",
    action=Store_prefix,
    required=False,
    default="",
    help="prefix for each saving file"
)

parser.add_argument(
    "-l", "--label",
    dest="label",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="label used in `adata.obs` for characterizing sample, useful for plotting"
)

parser.add_argument(
    "--method",
    dest="method",
    type=str,
    required=False,
    default="all",
    choices=["all", "ingest", "bbknn", "scanorama"],
    metavar="[all | ingest | bbknn | scanorama]",
    help="integration method to use (default: all)"
)

parser.add_argument(
    "-p", "--dim-pca",
    dest="dim_pca",
    type=int,
    required=False,
    default=50,
    metavar="INT",
    help="number of principal components (default: 50)"
)

parser.add_argument(
    "-c", "--dim-clustering",
    dest="dim_clustering",
    type=int,
    required=False,
    default=15,
    metavar="INT",
    help="number of principal components taken into account for clustering cells (default: 15)"
)

parser.add_argument(
    "-e", "--dim-integration",
    dest="dim_integration",
    type=int,
    required=False,
    default=2,
    metavar="INT",
    help="number of embedding dimensions (default: 2)"
)

parser.add_argument(
    "-z", "--zero-center",
    dest="zero_center",
    required=False,
    action="store_true",
    help="compute standard PCA from covariance matrix if `True`, otherwise omit zero-centering variables"
)

parser.add_argument(
    "--hvg",
    dest="hvg",
    required=False,
    action="store_true",
    help="select the most variable genes for dimension reduction with pca"
)

parser.add_argument(
    "-m", "--metric",
    dest="metric",
    type=str,
    required=False,
    default="euclidean",
    metavar="METRIC",
    help="metric used for knn and bbknn-based integration algorithms (default: euclidean)"
)

parser.add_argument(
    "-k", "--k-neighbors",
    dest="k_neighbors",
    type=int,
    required=False,
    default=20,
    metavar="INT",
    help="number of closest neighbors computed when performing KNN graph (default: 20)"
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
    "--add-legend",
    dest="legend",
    required=False,
    action="store_true",
    help="add legend to figures"
)

parser.add_argument(
    "--plot-3d",
    dest="plot_3d",
    required=False,
    action="store_true",
    help="plot figures in three dimensions"
)

parser.add_argument(
    "-j", "--jobs",
    dest="n_jobs",
    type=int,
    required=False,
    default=1,
    metavar="INT",
    help="number of process to use (default: 1)"
)

parser.add_argument(
    "-s", "--seed",
    dest="seed",
    type=int,
    required=False,
    default=None,
    metavar="FLOAT",
    help="random number generator"
)

parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    required=False,
    action="store_true",
    help="display information about running programm"
)

args = parser.parse_args()

if args.dim_pca < max(args.dim_clustering, args.dim_integration):
    raise argparse.ArgumentError(
        f"dim_pca ({args.dim_pca}) must be superior to dim_clustering ({args.dim_clustering}) and dim_integration ({args.dim_integration}), aborting."
    )

data_outpath = Path(f"{args.outpath}/tables")
fig_outpath = Path(f"{args.outpath}/figures")

if not data_outpath.exists():
    os.makedirs(data_outpath)
if not fig_outpath.exists():
    os.makedirs(fig_outpath)

default_seed = args.seed if args.seed else 10
section = Section(verbose = args.verbose)

print(f"Loading data...")

adata_d = odict()
adata_d["reference"] = sc.read_h5ad(args.infile_ref)
adata_d["interest"] = sc.read_h5ad(args.infile_interest)
if args.label:
    label = [adata.uns[args.label] for adata in adata_d.values()]
else:
    label = ["reference", "interest"]
valid_genes = list(set(adata_d["reference"].var.index).intersection(set(adata_d["interest"].var.index)))

for key in adata_d.keys():
    clean_adata(adata_d[key])
    adata_d[key].X = adata_d[key].layers["correct"]
    adata_d[key] = adata_d[key][:,valid_genes]
    sc.pp.highly_variable_genes(adata_d[key], layer="raw", flavor="seurat_v3", span=0.3, n_bins=20, n_top_genes=2000, inplace=True)

del valid_genes

if args.method=="ingest" or args.method=="all":

    print("Integration using ingest:")

    section("Computation of reference sample embedding components...", reset=True)
    sc.tl.pca(
        adata_d["reference"],
        zero_center=args.zero_center,
        n_comps=args.dim_pca,
        use_highly_variable=args.hvg,
        copy=False
    )
    sc.pp.neighbors(
        adata_d["reference"],
        n_neighbors=args.k_neighbors,
        n_pcs=args.dim_clustering,
        copy=False
    )
    sc.tl.umap(
        adata_d["reference"],
        n_components=args.dim_integration,
        random_state=default_seed
    )

    section("Integration of interest sample...")
    sc.tl.ingest(
        adata=adata_d["interest"],
        adata_ref=adata_d["reference"],
        obs=None,
        embedding_method=["pca", "umap"],
        n_jobs=args.n_jobs
    )
    try:
        adata = ad.concat(
            list(adata_d.values()),
            join="inner",
            label=args.label,
            keys=label,
            merge="same",
            uns_merge="same"
        )
    except:
        raise RuntimeError("Anndatas concatenation did not work, aborting")
    clean_adata(adata)
    sc.pp.neighbors(
        adata,
        n_neighbors=args.k_neighbors,
        use_rep="X_pca",
        n_pcs=args.dim_clustering,
        metric=args.metric,
        copy=False
    )
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added=f"leiden"
    )

    section("Plot of embedding components...")
    adt.pl.embedding_plot(
        adata,
        obs="condition",
        obsm="X_pca",
        xlabel=r"$\mathrm{PC_{1}}$",
        ylabel=r"$\mathrm{PC_{2}}$",
        outfile=Path(f"{fig_outpath}/ingest_pca"),
        add_legend=args.legend,
        s=2,
        alpha=1,
        lgd_params={
            "title":"conditions",
            "ncol":1,
            "markerscale":5,
            "frameon":True,
            "edgecolor":color.black,
            "shadow":False,
            "loc":"best"
        }
    )
    for cluster in ["condition", "leiden"]:
        fig, _ = adt.pl.embedding_plot(
            adata,
            obs=cluster,
            obsm="X_umap",
            xlabel=r"$\mathrm{UMAP_{1}}$",
            ylabel=r"$\mathrm{UMAP_{2}}$",
            zlabel=r"$\mathrm{UMAP_{3}}$",
            add_legend=args.legend,
            s=2,
            alpha=1,
            lgd_params={
                "title":"clusters" if cluster != "condition" else "conditions",
                "ncol":1,
                "markerscale":5,
                "frameon":True,
                "edgecolor":color.black,
                "shadow":False
            },
            n_components = 3 if args.dim_integration > 2 and args.plot_3d is True else 2,
            background_visible=False
        )
        plt.savefig(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}"))
        if args.dim_integration > 2 and args.plot_3d:
            pickle.dump(fig, open(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}.fig.pickle"), "wb"))

    section("Saving data...")
    adata.write_h5ad(filename=f"{data_outpath}/{args.prefix}ingest.h5ad", compression="gzip")
    del adata

if args.method=="bbknn" or args.method=="all":

    print("Integration using bbknn:")

    if "adata" not in globals():
        try:
            adata = ad.concat(
                list(adata_d.values()),
                join="inner",
                label=args.label,
                keys=label,
                merge="same",
                uns_merge="same"
            )
        except:
            raise RuntimeError("Anndatas concatenation did not work, aborting")
    clean_adata(
        adata,
        obs="leiden"
    )

    section("Computation of embedding components...", reset=True)
    sc.pp.highly_variable_genes(
        adata,
        layer="raw",
        flavor="seurat_v3",
        span=0.3,
        n_bins=20,
        n_top_genes=2000,
        inplace=True
    )
    sc.tl.pca(
        adata,
        zero_center=args.zero_center,
        n_comps=args.dim_pca,
        use_highly_variable=args.hvg,
        copy=False
    )

    section("Integration of embedding components...")
    sc.external.pp.bbknn(
        adata,
        batch_key=args.label,
        use_rep="X_pca",
        metric=args.metric,
        copy=False,
        neighbors_within_batch=args.k_neighbors,
        n_pcs=args.dim_clustering,
    )
    sc.tl.umap(
        adata,
        n_components=args.dim_integration,
        random_state=default_seed
    )
    sc.pp.neighbors(
        adata,
        n_neighbors=args.k_neighbors,
        use_rep="X_pca",
        n_pcs=args.dim_clustering,
        metric=args.metric,
        copy=False
    )
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added=f"leiden"
    )

    section("Plot of embedding components...")
    adt.pl.embedding_plot(
        adata,
        obs="condition",
        obsm="X_pca",
        xlabel=r"$\mathrm{PC_{1}}$",
        ylabel=r"$\mathrm{PC_{2}}$",
        outfile=Path(f"{fig_outpath}/bbknn_pca"),
        add_legend=args.legend,
        s=2,
        alpha=1,
        lgd_params={
            "title":"conditions",
            "ncol":1,
            "markerscale":5,
            "frameon":True,
            "edgecolor":color.black,
            "shadow":False,
        },
    )
    for cluster in ["condition", "leiden"]:
        fig, _ = adt.pl.embedding_plot(
            adata,
            obs=cluster,
            obsm="X_umap",
            xlabel=r"$\mathrm{UMAP_{1}}$",
            ylabel=r"$\mathrm{UMAP_{2}}$",
            zlabel=r"$\mathrm{UMAP_{3}}$",
            add_legend=args.legend,
            figwidth=6,
            s=2,
            alpha=1,
            lgd_params={
                "title":"clusters" if cluster != "condition" else "conditions",
                "ncol":1,
                "markerscale":5,
                "frameon":True,
                "edgecolor":color.black,
                "shadow":False
            },
            n_components = 3 if args.dim_integration > 2 and args.plot_3d is True else 2,
            background_visible=False
        )
        plt.savefig(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}"))
        if args.dim_integration > 2 and args.plot_3d:
            pickle.dump(fig, open(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}.fig.pickle"), "wb"))

    section("Saving data...")
    adata.write_h5ad(filename=f"{data_outpath}/{args.prefix}bbknn.h5ad", compression="gzip")
    del adata

if args.method=="scanorama" or args.method=="all":

    print("Integration using scanorama...")

    for key in adata_d.keys():
        clean_adata(adata_d[key])

    adata_l = list(adata_d.values())
    del adata_d

    section("Computation of integrated embedding components...", reset=True)
    adata_l = scanorama.correct_scanpy(
        adata_l,
        dimred=args.dim_pca,
        return_dimred=True
    )
    try:
        adata = ad.concat(
            adata_l,
            join="inner",
            label=args.label,
            keys=label,
            merge="same",
            uns_merge="same"
        )
        del adata_l
    except:
        raise RuntimeError("Anndatas concatenation did not work, aborting")
    sc.pp.neighbors(
        adata,
        n_neighbors=args.k_neighbors,
        use_rep="X_scanorama",
        n_pcs=args.dim_clustering,
        copy=False
    )
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added=f"leiden"
    )
    sc.tl.umap(
        adata,
        n_components=args.dim_integration,
        random_state=default_seed
    )

    section("Plot of embedding components...")
    adt.pl.embedding_plot(
        adata,
        obs="condition",
        obsm="X_pca",
        xlabel=r"$\mathrm{PC_{1}}$",
        ylabel=r"$\mathrm{PC_{2}}$",
        outfile=Path(f"{fig_outpath}/bbknn_pca"),
        add_legend=args.legend,
        s=2,
        alpha=1,
        lgd_params={
            "title":"conditions",
            "ncol":1,
            "markerscale":5,
            "frameon":True,
            "edgecolor":color.black,
            "shadow":False,
            "loc":"best"
        }
    )
    for cluster in ["condition", "leiden"]:
        fig, _ = adt.pl.embedding_plot(
            adata,
            obs=cluster,
            obsm="X_umap",
            xlabel=r"$\mathrm{UMAP_{1}}$",
            ylabel=r"$\mathrm{UMAP_{2}}$",
            zlabel=r"$\mathrm{UMAP_{3}}$",
            add_legend=args.legend,
            s=2,
            alpha=1,
            lgd_params={
                "title":"clusters" if cluster != "condition" else "conditions",
                "ncol":1,
                "markerscale":5,
                "frameon":True,
                "edgecolor":color.black,
                "shadow":False
            },
            n_components = 3 if args.dim_integration > 2 and args.plot_3d is True else 2,
            background_visible=False
        )
        plt.savefig(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}"))
        if args.dim_integration > 2 and args.plot_3d:
            pickle.dump(fig, open(Path(f"{fig_outpath}/{args.prefix}bbknn_umap_{cluster}.fig.pickle"), "wb"))

    section("Saving data...")
    adata.write_h5ad(filename=f"{data_outpath}/{args.prefix}scanorama.h5ad", compression="gzip")
