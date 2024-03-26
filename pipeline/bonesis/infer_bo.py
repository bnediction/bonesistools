#!/usr/bin/env python

import sys
import json
from pathlib import Path
from argparse import ArgumentParser

from tqdm import tqdm

import pandas as pd
import decoupler as dc

import networkx as nx
import bonesis
from bonesis.asp_encoding import clingo_encode

from plzf_rara_model import (
    important_nodes,
    bomodel,
    load_bin,
    collectri_to_grn
)

from utils.genesyn import GeneSynonyms

parser = ArgumentParser(
    prog="Boolean network inference",
    description="""From binarized meta-observations and specified trajectories,
    infer a Most Permissive Boolean Network""",
    usage="""python infer_bo.py [-h] <action> <path> --bin-metastate <path> [<args>]"""
)

parser.add_argument(
    "action",
    metavar="[filter_stage1 | filter_stage2 | one | one-min | one-sub]",
    choices=["filter_stage1", "filter_stage2", "one", "one-min", "one-sub"]
)

parser.add_argument(
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    metavar="PATH",
    help="output path"
)

parser.add_argument(
    "--bin-metastates",
    dest="bin_metastates",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="file with binarized clusters"
)

parser.add_argument(
    "--model-specification",
    dest="model_specification",
    type=lambda x: Path(x).resolve(),
    required=False,
    default="plzf_rara_model.txt",
    metavar="PATH",
    help="file with binarized clusters (default: plzf_rara_model.txt)"
)

parser.add_argument(
    "--filter-grn",
    dest="filter_grn",
    type=lambda x: Path(x).resolve(),
    required=False,
    metavar="PATH",
    help="file with one node per line"
)

parser.add_argument(
    "--force-nodes",
    dest="force_nodes",
    type=lambda x: Path(x).resolve(),
    required=False,
    metavar="PATH",
    help="json file with node list"
)

parser.add_argument(
    "--maximize-important-nodes",
    dest="maximize_important_nodes",
    required=False,
    action="store_true"
)

parser.add_argument(
    "--force-important-nodes",
    dest="force_important_nodes",
    required=False,
    action="store_true"
)

parser.add_argument(
    "--clingo-opt-strategy",
    dest="clingo_opt_strategy",
    type=str
)

parser.add_argument(
    "--quiet",
    dest="quiet",
    required=False,
    action="store_true"
)

args = parser.parse_args()

bonesis.settings["quiet"] = args.quiet

genename = GeneSynonyms()

collectri_db = dc.get_collectri(organism="mouse", split_complexes=True)
grn = collectri_to_grn(collectri_db, sign_label="weight", remove_pmid=True)
genename.graph_standardization(grn, copy=False)
if args.filter_grn:
    with open(args.filter_grn) as fp:
        nodes = [line.strip() for line in fp.readlines()]
    grn = grn.subgraph(nodes)

print(f"GRN has {len(grn.nodes)} nodes and {len(grn.edges)} edges", file=sys.stderr)

pkn_options = {
    "canonic": True,
    "maxclause": 8,
}
if args.action.startswith("filter"):
    pkn_options["canonic"] = False
if args.action == "filter_stage1":
    pkn_options["allow_skipping_nodes"] = True

meta_bin = load_bin(args.bin_metastates, gene_synonyms = genename)

pkn = bonesis.domains.InfluenceGraph(grn, **pkn_options)
bo = bonesis.BoNesis(pkn, meta_bin)
bomodel(bo, args.model_specification)

if args.action == "filter_stage1":
    bo.maximize_nodes()
    if args.force_nodes:
        with open(args.force_nodes) as fp:
            nodes = json.load(fp)
        for node in nodes:
            bo.custom(f"node({clingo_encode(node)}).")
    if args.force_important_nodes:
        for node in important_nodes:
            bo.custom(f"node({clingo_encode(node)}).")
    elif args.maximize_important_nodes:
        for node in important_nodes:
            bo.custom(f"important_node({clingo_encode(node)}).")
        bo.custom("#maximize { 1@100,N: important_node(N),node(N) }.")

    def interm_solution(nodes):
        with open(f"{args.outpath}/filter_stage1-last-model.json", "w") as fp:
            json.dump(list(sorted(nodes)), fp, indent=2)

    clingo_opt_strategy = args.clingo_opt_strategy or "bb,dec"
    view = bonesis.NodesView(bo, mode="optN", progress=tqdm,
                                intermediate_model_cb=interm_solution,
                                clingo_opt_strategy=clingo_opt_strategy)
    view.standalone(output_filename=f"{args.outpath}/filter_stage1.sh")
    solution = next(iter(view))
    for node in solution:
        print(node)

if args.action == "filter_stage2":
    bo.maximize_strong_constants()
    view = bonesis.NonStrongConstantNodesView(bo, mode="optN",
                                  clingo_opt_strategy="usc",
                                  clingo_options=["--opt-usc-shrink=inv"])
    view.standalone(output_filename=f"{args.outpath}/filter_stage2.sh")
    solution = next(iter(view))
    for node in solution:
        print(node)

def write_solution(solution, name):
    f = solution[1]
    f.save(f"{name}.bn")
    df = pd.DataFrame(solution[2])
    df.to_csv(f"{name}.csv")
    noi = set(f) - set(f.constants())
    with open(f"{name}.noi.txt", "w") as fp:
        fp.write("".join([f"{n}\n" for n in noi]))
    ig = f.influence_graph()
    nx.drawing.nx_pydot.write_dot(ig, f"{name}.dot")

if args.action == "one":
    view = bonesis.InfluenceGraphView(bo, extra=("boolean-network", "configurations"))
    solution = next(iter(view))
    write_solution(solution, f"{args.outpath}/bn-1")

if args.action == "one-min":
    bo.custom("edge(A,B) :- clause(B,_,A,_). #minimize { 1@1,A,B: edge(A,B) }.")
    bo.custom("#maximize { 1@10,N: constant(N) }.")
    view = bonesis.InfluenceGraphView(bo, mode="optN", clingo_opt_strategy="usc",
                                      extra=("boolean-network",
                                             "configurations"),
                                      progress=tqdm)
    view.standalone(output_filename=f"{args.outpath}/one-min.sh")
    solution = next(iter(view))
    write_solution(solution, "min-1")

if args.action == "one-sub":
    view = bonesis.InfluenceGraphView(bo, solutions="subset-minimal", extra=("boolean-network", "configurations"))
    view.standalone(output_filename=f"{args.outpath}/one-sub.sh")
    solution = next(iter(view))
    write_solution(solution, f"{args.outpath}/sub-1")
