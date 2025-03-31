#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Union

from pathlib import Path

import os
import argparse
import re

import anndata as ad
import pandas as pd
import scanpy as sc
import anndatatools as adt

PathLike = Union[str,Path]

parser = argparse.ArgumentParser(
    prog="Convert single-omics counting file format",
    description="""Convert single-omics counting data into the desired format.
    In case of 10x sparse matrix format, it is a directory containing three files:
    - matrix.mtx.gz (sparse matrix in the Market Exchange MEX format) -- also named coordinate list format, which corresponds to compressed reordered sparse counting data)
    - barcodes.tsv.gz (information about each cell)
    - features.tsv.gz (information about each gene)""",
    usage="python adata_conversion.py [-h] <PATH|FILE> <PATH|FILE> --from <h5ad|loom|10x> --to <h5ad|loom|csvs> [--metadata <KEY=VALUE ...> <args>]",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    "input",
    type=lambda x: Path(x).resolve(),
    metavar="PATH|FILE",
    help="input data directory or file"
)

parser.add_argument(
    "output",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="output data directory or file"
)

parser.add_argument(
    "--from",
    dest="from",
    type=str,
    choices=["h5ad", "loom", "10x"],
    metavar="[h5ad | loom | 10x]",
    required=True,
    help="matrix data input format"
)

parser.add_argument(
    "--to",
    dest="to",
    type=str,
    choices=["h5ad", "loom", "csv", "csvs"],
    metavar="[h5ad | loom | csv | csvs]",
    required=True,
    help="matrix data input format"
)

parser.add_argument(
    "--layer",
    dest="layer",
    type=str,
    required=False,
    default=None,
    metavar="LITERAL",
    help="layer saved if `--to=csv` (if not specified, save adata.X)"
)

parser.add_argument(
    "--remove-positions",
    dest="remove_positions",
    required=False,
    action="store_true",
    help="remove chromosome, position on it and strand directions for each gene"
)

parser.add_argument(
    "--metadata",
    dest="metadata",
    type=str,
    nargs="*",
    required=False,
    default=None,
    metavar="KEY=VALUE",
    help="dataset-related metadata"
)

parser.add_argument(
    "--genename-standardization",
    dest="genename_standardization",
    required=False,
    action="store_true",
    help="convert gene names by their NCBI reference names"
)

parser.add_argument(
    "--compression",
    dest="compression",
    required=False,
    action="store_true",
    help="output file compression (if output format is h5ad)"
)

def add_metadata(
    adata: ad.AnnData,
    **metadata
) -> None:
    for k, v in metadata.items():
        adata.uns[k] = v
    
args = parser.parse_args()

if args.__getattribute__("from") == args.__getattribute__("to"):
    raise argparse.ArgumentError("Argument --from and --to must be different")

if args.__getattribute__("to") == "csvs":
    os.makedirs(name=args.output, exist_ok=True)
else:
    os.makedirs(name=os.path.dirname(args.output), exist_ok=True)

if args.__getattribute__("from") == "h5ad":
    adata = sc.read_h5ad(filename=args.input)
elif args.__getattribute__("from") == "loom":
    adata = sc.read_loom(filename=args.input)
elif args.__getattribute__("from") == "10x":
    adata = sc.read_10x_mtx(path=args.input)

adata.obs.index = pd.Index(map(lambda barcode: re.sub("[^ATCG]","",re.sub("^.*:","",barcode)), adata.obs.index))

if args.remove_positions:
    for column in ["Chromosome", "Start", "End", "Strand"]:
        if column in adata.var.columns:
            del adata.var[column]

if args.metadata:
    split = lambda metadata: [metadatum.split("=") for metadatum in metadata]
    metadata_d = {info[0]: info[1] for info in split(args.metadata)}
    add_metadata(adata, **metadata_d)

if args.genename_standardization:
    adata.var["symbol"] = list(adata.var.index)
    for alias_type in ["genename","geneid","ensemblid"]:
        adt.pp.set_ncbi_reference_name(adata, annotations="var", in_alias_type=alias_type, copy=False)
    adata = adt.pp.var_names_merge_duplicates(adata, var_names_column="symbol")

if args.__getattribute__("to") == "h5ad":
    adata.write_h5ad(filename=args.output, compression="gzip" if args.compression else None)
elif args.__getattribute__("to") == "loom":
    adata.write_loom(filename=args.output, write_obsm_varm=True)
elif args.__getattribute__("to") == "zarr":
    adata.write_zarr(store=args.output)
elif args.__getattribute__("to") == "csv":
    adt.tl.anndata_to_dataframe(
        adata=adata,
        layer=args.layer
    ).to_csv(
        path_or_buf=args.output,
        sep=",",
        index=True
    )
elif args.__getattribute__("to") == "csvs":
    adata.write_csvs(dirname=args.output, sep=",")
    adt.tl.to_csv_or_mtx(
        adata=adata,
        filename=Path(f"{args.output}/matrix")
    )
    if adata.layers.keys():
        os.makedirs(name=Path(f"{args.output}/layers"), exist_ok=True)
        for layer in adata.layers:
            adt.tl.to_csv_or_mtx(
                adata=adata,
                filename=Path(f"{args.output}/layers/{layer}"),
                layer=layer
            )
