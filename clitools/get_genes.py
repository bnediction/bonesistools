#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import os, argparse
from pathlib import Path

from databases.genesyn import GeneSynonyms

import scanpy as sc

parser = argparse.ArgumentParser(
    prog="gene name list",
    description="""Retrieve gene names.""",
    usage="""python genename.py [-h] <FILE> <FILE> [--axis <obs|var> --standardization]"""
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="h5ad or loom file"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="txt file containing list of genes"
)

parser.add_argument(
    "--axis",
    dest="axis",
    choices=["obs", "var"],
    default="var",
    required=False,
    metavar="[obs|var]",
    help="axis corresponding to genes (default: obs)"
)

parser.add_argument(
    "--standardization",
    dest="standardization",
    required=False,
    action="store_true",
    help="standardize gene name with their NCBI reference name"
)

args = parser.parse_args()

if not Path(os.path.dirname(args.outfile)).exists():
    os.makedirs(os.path.dirname(args.outfile))

if str(args.infile).endswith(".h5ad"):
    adata = sc.read_h5ad(args.infile)
elif str(args.infile).endswith(".loom"):
    adata = sc.read_loom(args.infile)
else:
    raise IOError(f"incorrect format (h5ad or loom format expected)")

if args.axis == "obs":
    gene_list = list(adata.obs.index)
else:
    gene_list = list(adata.var.index)

if args.standardization:
    gene_list = GeneSynonyms()(gene_list)

with open(args.outfile, "w") as file:
    for gene in gene_list:
        file.write(f"{gene}\n")
