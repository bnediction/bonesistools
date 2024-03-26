#!/usr/bin/env python

import argparse, os

from pathlib import Path

import pandas as pd, json

def file2signatures(file):
    signatures_d = dict()
    for sheetname, signature in file.items():
        if not sheetname=="Description":
            cell_type = sheetname.split(".txt", 1)[0]
            gene_symbols = [gene for gene in list(signature["Gene Symbol"]) if isinstance(gene, str)]
            signatures_d[cell_type] = gene_symbols
    return signatures_d

def df2signatures(df):
    signatures_d = dict()
    for cell_type, signatures in df.items():
        gene_symbols = [gene for gene in signatures if isinstance(gene, str)]
        signatures_d[cell_type] = gene_symbols
    return signatures_d

parser = argparse.ArgumentParser(
    prog="Load signatures data",
    description="""Load signatures data from two files,
    one in a table format and the other one in list format.""",
    usage="python load_signatures.py <args>"
)

parser.add_argument(
    "-t", "--table-infile",
    dest="table_infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    help="path to table signatures file"
)

parser.add_argument(
    "-l", "--list-infile",
    dest="list_infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    help="path to list signatures file"
)

parser.add_argument(
    "-o", "--outfile",
    dest="outfile",
    type=lambda x: Path(x).resolve(),
    required=True,
    help="output file"
)

args = parser.parse_args()

outpath = os.path.dirname(args.outfile)
if not outpath:
    os.makedirs(outpath)

table_signatures_df = pd.read_excel(io=args.table_infile, sheet_name=None)
table_signatures_d = file2signatures(table_signatures_df)

list_signatures_df = pd.read_excel(io=args.list_infile, sheet_name=0)
list_signatures_df.columns = list(list_signatures_df.iloc[0])
list_signatures_df.drop([0, 1], axis=0, inplace=True)
list_signatures_d = df2signatures(list_signatures_df)

signatures_d = {
    **table_signatures_d,
    **list_signatures_d
}

with open(f"{args.outfile}", "w") as file:
    json.dump(signatures_d, file, indent=1)
