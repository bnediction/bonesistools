#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Optional, Union, Sequence
from collections import OrderedDict
import boolpy as bp

import sys
import argparse
import json
from pathlib import Path
from utils.argtype import Range
from utils.stdout import Section

import pandas as pd
import decoupler as dc
import numpy as np

import itertools
import networkx as nx
import graphtools as gtl

from utils.genesyn import GeneSynonyms

def collectri_to_grn(
    collectri: pd.DataFrame,
    sign_label: str = "weight",
    remove_pmid: bool = False
) -> nx.MultiDiGraph:
    if sign_label is not None:
        collectri = collectri.rename(columns = {sign_label:"sign"})
    if remove_pmid is True:
        remove_pmid = "PMID" in collectri.columns
    return nx.from_pandas_edgelist(
        df = collectri.drop("PMID", axis=1) if remove_pmid else collectri,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph
    )

def gene_removal(df: pd.DataFrame, graph: nx.Graph, copy: bool=True) -> Union[pd.DataFrame, None]:
    df = df.copy() if copy is True else df
    genes_to_remove = list()
    for gene in df.columns:
        if gene not in graph.nodes:
            genes_to_remove.append(gene)
    df.drop(labels=genes_to_remove, axis="columns", inplace=True)
    return df if copy is True else None

def sign_likelihood(
    interaction_scores: dict,
    gene_set: Optional[Sequence[str]] = None,
    minimum_path_number: int = 3,
    relative_threshold: float = 0.75,
    enable_loop: bool = False
):

    if not (0 < relative_threshold < 1):
        raise ValueError("`relative_threshold` must be between 0 and 1: `relative_threshold` = {relative_threshold}")
    if gene_set is None:
        gene_set = set(interaction_scores.keys())
    else:
        gene_set = set(gene_set).intersection(set(interaction_scores.keys()))

    interaction_signs = {gene: dict() for gene in gene_set}

    for u, v in itertools.combinations(gene_set, 2):
        
        from_u = interaction_scores[u][v]
        from_v = interaction_scores[v][u]
        _is_source = [False, False]
        _sign = [0, 0]
        
        if from_u.path_number >= minimum_path_number:
            if abs(from_u.score) / from_u.maxscore >= from_u.maxscore * relative_threshold:
                _is_source[0] = True
                _sign[0] = 1 if from_u.score > 0 else -1
        if from_v.path_number >= minimum_path_number:
            if abs(from_v.score) / from_v.maxscore >= from_v.maxscore * relative_threshold:
                _is_source[1] = True
                _sign[1] = 1 if from_v.score > 0 else -1
        
        if enable_loop is True:
            if _is_source[0] is True:
                interaction_signs[u][v] = _sign[0]
            if _is_source is True:
                interaction_signs[v][u] = _sign[1]
        else:
            if _is_source[0] ^ _is_source[1]:
                if _is_source[0] is True:
                    interaction_signs[u][v] = _sign[0]
                elif _is_source[1] is True:
                    interaction_signs[v][u] = _sign[1]
    
    return interaction_signs

parser = argparse.ArgumentParser(
    prog="computation of inter-cluster velocities",
    description="""compute velocity between cluster with respect to binarized meta-observations""",
    usage=""""python velocity.py [-h] -i <path> <path> [<args>]"""
)

parser.add_argument(
    "infile",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="infile in csv format"
)

parser.add_argument(
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="output path"
)

parser.add_argument(
    "--depth", "--radius",
    dest="radius",
    type=int,
    required=False,
    default=3,
    metavar="INT",
    help="maximum path length between a source and a target (default: 3)"
)

parser.add_argument(
    "--min-path-number",
    dest="min_path_number",
    type=int,
    required=False,
    default=1,
    metavar="INT",
    help="minimum number of paths for gene pairwise required for considering a gene as being a source (default: 3)"
)

parser.add_argument(
    "--base",
    dest="base",
    type=int,
    required=False,
    default=2,
    metavar="INT",
    help="base in the non-exponential weighting function (default: 2)"
)

parser.add_argument(
    "--relative-threshold",
    dest="threshold",
    type=float,
    action=Range,
    min=0.,
    max=1.,
    required=False,
    default=0.75,
    help="percentage of the maximum path score above which the path score must be for being consider as source-target gene pairwise (default: 0.75)"
)

parser.add_argument(
    "--enable-loop",
    dest="enable_loop",
    required=False,
    action="store_true",
    help="allow a gene pairwise to be mutually influenced by the other one"
)

parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    required=False,
    action="store_true",
    help="display information about running programm"
)

args = parser.parse_args()
if args.base <= 1:
    raise argparse.ArgumentError("incorrect value for `base` argument : {args.base}")

section = Section(verbose = args.verbose)
nexponential_fun = lambda base, radius: 1 / base**np.arange(0, radius)
bdc = bp.algebra.BooleanDifferentialCalculus()

print(f"Loading data...")

meta_bin = pd.read_csv(args.infile, index_col=0)

collectri_db = dc.get_collectri(organism="mouse", split_complexes=True)
grn = collectri_to_grn(collectri_db, sign_label="weight", remove_pmid=True)

if args.verbose:
    print(f"\tgrn: {len(grn.nodes)} genes; {len(grn.edges)} interactions", file=sys.stderr)

gene_synonyms = GeneSynonyms()
gene_synonyms(data=meta_bin, axis=1, copy=False)
gene_synonyms(data=grn, copy=False)
gene_set_before_cleaning = set(meta_bin.columns)
gene_removal(meta_bin, grn, copy=False)
gene_set = set(meta_bin.columns)

if args.verbose:
    print(f"\tdataframe: {len(gene_set_before_cleaning)} genes; {len(gene_set_before_cleaning)- len(gene_set)}/{len(gene_set_before_cleaning)} genes removed (no matching with grn genes)", file=sys.stderr)

print("Successors checking...")

section("Path extraction using depth-first extraction algorithm")
interaction_scores = gtl.grn.scoring(
    graph=grn,
    weights=nexponential_fun(base=args.base, radius=args.radius),
    radius=args.radius,
    gene_set=gene_set
)

section("Sign likelihood between gene pairwise")
interaction_signs = sign_likelihood(
    interaction_scores=interaction_scores,
    gene_set=gene_set,
    minimum_path_number=args.min_path_number,
    relative_threshold=args.threshold,
    enable_loop=args.enable_loop
)

with open(f"{args.outpath}/sign_likelihood.json", "w") as outfile:
    json.dump(interaction_signs, outfile)

section("Predecessor test using differential boolean calculus")

score_matrix = OrderedDict({condition: {} for condition in meta_bin.index})
for c1, c2 in itertools.product(meta_bin.index, repeat=2):
    score_matrix[c1][c2] = 0

for source, targets in interaction_signs.items():
    for target, sign in targets.items():
        pair_df = meta_bin.loc[:, [source, target]]
        for c1, c2 in itertools.product(meta_bin.index, repeat=2):
            _predecessor = bdc.pairwise_predecessor_test(
                source_v1 = pair_df.loc[c1, source],
                source_v2 = pair_df.loc[c2, source],
                target_v1 = pair_df.loc[c1, target],
                target_v2 = pair_df.loc[c2, target],
                sign=sign
            )
            if _predecessor is True:
                score_matrix[c1][c2] += 1
            elif _predecessor is False:
                score_matrix[c1][c2] -= 1
            else:
                pass

score_df = pd.DataFrame.from_dict(score_matrix, orient="index")

print("Saving data...")

score_df.to_csv(f"{args.outpath}/pairwise_predecessor_scores.csv", sep=",", index=True)

if args.verbose:
    print(f"\n{score_df}\n")