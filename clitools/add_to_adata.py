#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Optional, Union

from pathlib import Path

import os
import argparse

import anndata as ad
import pandas as pd
import scanpy as sc
import anndatatools as adt

PathLike = Union[str,Path]

parser = argparse.ArgumentParser(
    prog="Add information to anndata object",
    description="""Add observation and/or variable-related information from csv or tsv files to anndata object.""",
    usage="python add_to_adata.py [-h] <FILE> <FILE> [--obs <FILE] [--variable <FILE]",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="input file (loom or h5ad format)"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="output file (loom or h5ad format)"
)

parser.add_argument(
    "--obs",
    dest="obs",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=None,
    metavar="PATH",
    help="file containing new data about observations (csv or tsv format)"
)

parser.add_argument(
    "--obs-type",
    dest="obs_type",
    type=type,
    required=False,
    choices=[str, int, float, bool],
    default=None,
    metavar="TYPE",
    help="observation type (str, int, float, bool)"
)

parser.add_argument(
    "--var-type",
    dest="var_type",
    type=type,
    required=False,
    choices=[str, int, float, bool],
    default=None,
    metavar="TYPE",
    help="observation type (str, int, float, bool)"
)

parser.add_argument(
    "--var",
    dest="var",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=None,
    metavar="PATH",
    help="file containing new data about variables (csv or tsv format)"
)

parser.add_argument(
    "--sep",
    dest="sep",
    type=str,
    required=False,
    default=",",
    metavar="CHAR",
    help="field delimiter for csv format (default: `,`)"
)

args = parser.parse_args()

if args.obs is None and args.var is None:
    argparse.ArgumentError("missing argument (--var or --obs is required")
else:
    obs_df = pd.read_csv(args.obs, index_col=0, sep=args.sep, dtype=args.obs_type) if args.obs else None
    var_df = pd.read_csv(args.var, index_col=0, sep=args.sep, dtype=args.var_type) if args.var else None

if str(args.infile).endswith("h5ad"):
    adata = sc.read_h5ad(filename=args.infile)
elif str(args.infile).endswith("loom"):
    adata = sc.read_loom(filename=args.infile)
else:
    raise argparse.ArgumentError("unable to synchronously open infile (required format: h5ad or loom)")

for column_name in obs_df.columns:
    if column_name in adata.obs.columns:
        del adata.obs[column_name]
adata.obs = adata.obs.merge(how='left',right=obs_df, left_index=True, right_index=True)

for column_name in var_df.columns:
    if column_name in adata.var.columns:
        del adata.var[column_name]
adata.var = adata.var.merge(how='left',right=var_df, left_index=True, right_index=True)

if str(args.outfile).endswith("h5ad"):
    adata.write_h5ad(filename=args.outpath, compression="gzip")
elif str(args.outfile).endswith("loom"):
    adata.write_loom(filename=args.outpath, write_obsm_varm=True)
else:
    raise argparse.ArgumentError("unable to synchronously create outfile (required format: h5ad or loom)")
