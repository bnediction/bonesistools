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
    "--organism",
    dest="organism",
    choices=["mouse","human","escherichia coli"],
    default="mouse",
    required=False,
    metavar="[mouse | human | escherichia coli]",
    help="gene-related organism (default: mouse)"
)

parser.add_argument(
    "--in-alias-type",
    dest="in_alias_type",
    default="genename",
    required=False,
    metavar="[genename | geneid | ensemblid | <database>]",
    help="input gene alias type (default: genename)"
)

parser.add_argument(
    "--out-alias-type",
    dest="out_alias_type",
    default="referencename",
    required=False,
    metavar="[referencename | geneid | ensemblid | <database>]",
    help="output gene alias type (default: referencename)"
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

genesynonyms = GeneSynonyms(organism=args.organism)

if file_extension == "txt":
    with open(args.infile, "r") as file:
        gene_list = [line.replace("\n","") for line in file]
    gene_list = genesynonyms(gene_list, in_alias_type=args.in_alias_type, out_alias_type=args.out_alias_type)
    with open(args.outfile, "w") as file:
        for gene in gene_list:
            file.write(f"{gene}\n")
elif file_extension == "csv" or file_extension == "tsv":
    output = pd.read_csv(args.infile, index_col=0, sep=args.sep)
    genesynonyms(
        output,
        in_alias_type=args.in_alias_type,
        out_alias_type=args.out_alias_type,
        axis=args.axis,
        copy=False
    )
    output.to_csv(args.outfile, sep=args.sep)
else:
    raise OSError(f"extension not supported for `{args.infile}`: available extension are txt, csv and tsv")
