#!/usr/bin/env python

from pathlib import Path

import os
import argparse
import scanpy as sc

import anndatatools as adt

parser = argparse.ArgumentParser(
    prog="Loading 10X sparse matrix sc-RNAseq data",
    description="""Convert sc-RNAseq data from 10X sparse matrix format into hdf5 format. The structure of 10X sparse matrix format is a directory containing three files:
    - matrix.mtx.gz (sparse matrix in the Market Exchange MEX format) -- also named coordinate list format, which corresponds to compressed reordered sparse counting data)
    - barcodes.tsv.gz (information about each cell)
    - features.tsv.gz (information about each gene)""",
    usage="python load_10X.py [-h] <PATH> <FILE> [--sample-info <KEY=VALUE ...>]",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    "inpath",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="directory to 10X sparse matrix data"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="outfile in csv format"
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

args = parser.parse_args()

outpath = Path(os.path.split(args.outfile)[0]).resolve()

if not outpath.exists():
    os.makedirs(outpath)

adata = sc.read_10x_mtx(path=args.inpath)
adata.raw = adata
adata.var["symbol"] = list(adata.var.index)

for metadatum in args.sample_info:
    key, value = metadatum.split("=")
    adata.uns[key] = value

adt.pp.set_ncbi_reference_name(adata, annotations="var", copy=False)
adata = adt.pp.var_names_merge_duplicates(adata, var_names_column="symbol")

adata.write_h5ad(filename=args.outfile, compression="gzip")