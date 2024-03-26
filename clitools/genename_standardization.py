#!/usr/bin/env python

import os, argparse
from pathlib import Path

from utils.genesyn import GeneSynonyms

import pandas as pd

parser = argparse.ArgumentParser(
    prog="gene name standardization",
    description="""Update gene names by their NCBI reference name.""",
    usage="""python genename_standardization.py [-h] <path> <path> [<args>]"""
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="infile in csv format"
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
    default=",",
    required=False,
    metavar="CHAR",
    help="field delimiter for the outfile (default: `,`)"
)

parser.add_argument(
    "--axis",
    dest="axis",
    choices=["index", "columns"],
    default="index",
    required=False,
    metavar="[index | columns]",
    help="axis to change in dataframe instance (default: index)"
)

args = parser.parse_args()

if not Path(os.path.dirname(args.outfile)).exists:
    os.makedirs(os.path.dirname(args.outfile))

df = pd.read_csv(args.infile, index_col=0, sep=args.sep)
GeneSynonyms()(df, axis=args.axis, copy=False)
df.to_csv(args.outfile, sep=args.sep)
