#!/usr/bin/env python

from pathlib import Path

import argparse

import ginsim
import biolqm

parser = argparse.ArgumentParser(
    prog="Conversion to a pure Boolean network",
    description="""Convert regulatory network (format zginml) into pure Boolean network (format bnet)""",
    usage="python zginml_to_bnet.py [-h] <FILE> [-o <FILE>]"
)

parser.add_argument(
    dest="infile",
    type=lambda x: Path(x).resolve(),
    metavar="FILE",
    help="zginml file"
)

parser.add_argument(
    "-o",
    dest="outfile",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=None,
    metavar="FILE",
    help="bnet file (if not specified, send it to stdout)"
)

args = parser.parse_args()

lrg = ginsim.load(str(args.infile))
bn = biolqm.to_minibn(lrg.getModel())

if args.outfile:
    with open(args.outfile, "w") as file:
        file.write(str(bn))
else:
    print(bn)
