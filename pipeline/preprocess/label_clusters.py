#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import os, argparse
from pathlib import Path
from utils.argtype import Store_dict

import anndata as ad

parser = argparse.ArgumentParser(
    prog="labeling of clusters",
    description="""From sc-rnaSeq data recorded in the hdf5 format,
    rename labels using user-defined names.""",
    usage="""python cluster_labeling.py [-h] -i <PATH> -o <PATH> -c <LITERAL> -n <LITERAL=LITERAL [LITERAL=LITERAL ...]> [<args>]"""
)

parser.add_argument(
    "-i", "--infile",
    dest="infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad infile (including file)"
)

parser.add_argument(
    "-o", "--outfile",
    dest="outfile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad outfile (including file)"
)

parser.add_argument(
    "-c", "--column",
    dest="column",
    type=str,
    required=True,
    metavar="LITERAL",
    help="name of the column in adata.obs from which user want to redefine category names"
)

parser.add_argument(
    "-n", "--name",
    dest="labels",
    action=Store_dict,
    required=True,
    nargs="+",
    help="mapping between old and new names for labels"
)

parser.add_argument(
    "--obsm",
    dest="obsm",
    required=False,
    default=None,
    metavar="LITERAL",
    help="""ndarray name stored in .obsm[`obsm`] used for plotting figure"""
)

args = parser.parse_args()

if not Path(os.path.dirname(args.outfile)).exists():
    os.makedirs(Path(os.path.dirname(args.outfile)))

print(f"Loading data...")
adata = ad.read_h5ad(args.infile)

print(f"Rename categories...")
if args.column not in adata.obs:
    raise KeyError(f"adata.obsm[`{args.column}`] does not exist.")
elif not hasattr(adata.obs[args.column], "cat"):
    raise ValueError("values in adata.obs[`{args.column}`] are not derived from a Categorical type.")
else:
    adata.obs[args.column] = adata.obs[args.column].cat.rename_categories(args.labels)

print("Saving data...")
adata.write_h5ad(filename=args.outfile, compression="gzip")
