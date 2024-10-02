#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import os
import argparse
import re

import pandas as pd, scanpy as sc

import anndatatools as adt

parser = argparse.ArgumentParser(
    prog="Convert multi-omics counting file into h5ad format",
    description="""Convert sc-RNAseq data from 10X sparse matrix or loom format into hdf5 format. In case of 10X sparse matrix format, it is a directory containing three files:
    - matrix.mtx.gz (sparse matrix in the Market Exchange MEX format) -- also named coordinate list format, which corresponds to compressed reordered sparse counting data)
    - barcodes.tsv.gz (information about each cell)
    - features.tsv.gz (information about each gene)""",
    usage="python convert_to_h5ad.py [-h] <PATH|FILE> <FILE> [--sample-info <KEY=VALUE ...>]",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    "inpath",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="10X sparse matrix data directory or loom file"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="outfile in h5ad format"
)

parser.add_argument(
    "--sample-info",
    dest="sample_info",
    type=str,
    nargs="*",
    required=False,
    default=None,
    metavar="KEY=VALUE",
    help="sample metadata"
)

parser.add_argument(
    "--remove-positions",
    dest="remove_positions",
    required=False,
    action="store_true",
    help="remove chromosome, position on it and strand direction for each gene"
)

args = parser.parse_args()

outpath = Path(os.path.split(args.outfile)[0]).resolve()

if not outpath.exists():
    os.makedirs(outpath)

if os.path.isfile(args.inpath):
    adata = sc.read_loom(filename=args.inpath)
else:
    adata = sc.read_10x_mtx(path=args.inpath)

adata.raw = adata
adata.obs.index = pd.Index(map(lambda barcode: re.sub("[^ATCG]","",re.sub("^.*:","",barcode)), adata.obs.index))
adata.var["symbol"] = list(adata.var.index)
if "Accession" in adata.var.columns:
    adata.var.rename(columns={"Accession":"ensemblid"}, inplace=True)
if args.remove_positions:
    for column in ["Chromosome", "Start", "End", "Strand"]:
        if column in adata.var.columns:
            del adata.var[column]

for metadatum in args.sample_info:
    key, value = metadatum.split("=")
    adata.uns[key] = value

for alias_type in ["genename","geneid","ensemblid"]:
    adt.pp.set_ncbi_reference_name(adata, annotations="var", in_alias_type=alias_type, copy=False)
adata = adt.pp.var_names_merge_duplicates(adata, var_names_column="symbol")

adata.write_h5ad(filename=args.outfile, compression="gzip")