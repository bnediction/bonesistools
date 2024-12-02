#!/usr/bin/env python

from pathlib import Path

import argparse

from colomoto import minibn

parser = argparse.ArgumentParser(
    prog="Boolean Network nodes",
    description="""Retrieve nodes from a Boolean Network (bn)""",
    usage="python boolean_network_nodes.py [-h] <FILE ...> [-o <FILE>]"
)

parser.add_argument(
    dest="infile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="bnet file"
)

parser.add_argument(
    "-o",
    dest="outfile",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=None,
    metavar="FILE",
    help="txt file (if not specified, send it to stdout)"
)

args = parser.parse_args()

f = minibn.BooleanNetwork.load(args.infile)
nodes = list(f.keys())

if args.outfile:
    with open(args.outfile, "w") as file:
        for node in nodes:
            file.write(f"{node}\n")
else:
    for node in nodes:
        print(node)
