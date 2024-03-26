#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import os, argparse

import anndata as ad, pandas as pd

import numpy as np
from scipy.sparse import issparse

parser = argparse.ArgumentParser(
    prog="anndata to csv conversion",
    description="""Convert anndata object to dataframe object by taking counts
    and optionaly supplemental observation and save it in the specified location in csv format.""",
    usage="""python anndata2csv.py [-h] <path> <path> [<args>]"""
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="infile in h5ad format"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="outfile in csv format"
)

parser.add_argument(
    "-s", "--sep",
    dest="sep",
    type=str,
    required=False,
    default=",",
    metavar="CHAR",
    help="field delimiter for the outfile (default: `,`)"
)

parser.add_argument(
    "--layer",
    dest="layer",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="layer saved instead of attribute `X` in anndata instance"
)

parser.add_argument(
    "--observations",
    dest="obs",
    nargs="*",
    metavar="LITERAL",
    action="store",
    help="supplemental information to save from `obs` attribute in anndata instance"
)

parser.add_argument(
    "--exponential",
    dest="exp",
    required=False,
    action="store_true",
    help="perform exponential transformation on counts"
)

args = parser.parse_args()

if len(args.sep) != 1:
    raise argparse.ArgumentTypeError("`sep` argument must be a single character")

if not Path(os.path.dirname(args.outfile)).exists():
    os.makedirs(Path(os.path.dirname(args.outfile)))

print("Loading data...")
adata = ad.read_h5ad(args.infile)

print("Converting data...")
if args.layer:
    if issparse(adata.layers[args.layer]):
        df = pd.DataFrame(adata.layers[args.layer].toarray(), index=adata.obs.index, columns = adata.var.index)
    else:
        df = pd.DataFrame(adata.layers[args.layer], index=adata.obs.index, columns = adata.var.index)
elif issparse(adata.X):
    df = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns = adata.var.index)
else:
    df = pd.DataFrame(adata.X, index=adata.obs.index, columns = adata.var.index)
    
if args.exp:
    if "log1p" in adata.uns_keys() and adata.uns["log1p"].get('base') is not None:
        df = np.expm1(df * np.log(adata.uns['log1p']['base']))
    else:
        df = np.expm1(df)

if args.obs is not None:
    df.loc[:,args.obs] = adata.obs.loc[:,args.obs]

print("Saving data...")
print(df)
df.to_csv(args.outfile, sep=args.sep, index=True)
