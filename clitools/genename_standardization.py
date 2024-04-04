#!/usr/bin/env python

import warnings

import os, argparse
from pathlib import Path

from databases.genesyn import GeneSynonyms

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
    help="infile in txt, csv or tsv format"
)

parser.add_argument(
    "outfile",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="outfile in same format as infile"
)

parser.add_argument(
    "-s", "--sep",
    dest="sep",
    type=str,
    default=",",
    required=False,
    metavar="CHAR",
    help="field delimiter for the outfile if extension infile is csv or tsv (default: `,`)"
)

parser.add_argument(
    "--axis",
    dest="axis",
    choices=["index", "columns"],
    default="index",
    required=False,
    metavar="[index | columns]",
    help="axis to change in dataframe instance if extension infile is csv or tsv (default: index)"
)

parser.add_argument(
    "-q", "--quiet",
    dest="quiet",
    required=False,
    action="store_true",
    help="hidden information about running programm"
)

args = parser.parse_args()

if args.quiet:
    warnings.filterwarnings("ignore")

if not Path(os.path.dirname(args.outfile)).exists:
    os.makedirs(os.path.dirname(args.outfile))

file_extension = str(args.infile).split(".")[-1]

if file_extension == "txt":
    with open(args.infile, "r") as file:
        gene_list = file.readlines()
    gene_list = GeneSynonyms()(gene_list)
    with open(args.outfile, "w") as file:
        for gene in gene_list:
            file.write(f"{gene}\n")
elif file_extension == "csv" or file_extension == "tsv":
    output = pd.read_csv(args.infile, index_col=0, sep=args.sep)
    GeneSynonyms()(output, axis=args.axis, copy=False)
    output.to_csv(args.outfile, sep=args.sep)
else:
    raise OSError(f"extension not supported for `{args.infile}`: available extension are txt, csv and tsv")
