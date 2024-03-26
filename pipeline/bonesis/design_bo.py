#!/usr/bin/env python

from pathlib import Path

import argparse

from itertools import combinations

parser = argparse.ArgumentParser(
    prog="Convert a text file describing trajectories into comprehensible text for bonesis package",
    description="""conversion of a text file describing the trajectories where each line is in the form: `node_1 -> ... -> node_k`. \
    The output stream provides specifications comprehensible for bonesis.""",
    usage="python design_bo.py [-h] <path>"
)

parser.add_argument(
    dest="infile",
    type=lambda x: Path(x).resolve(),
    metavar="<path>",
    help="path to text infile"
)

args = parser.parse_args()

def read_trajectories(
    file: Path
) -> list:
    trajectories = list()
    with open(file, "r") as file:
        for line in file:
            trajectory = line.replace("\n", "").split(" -> ")
            trajectories.append(trajectory)
    return trajectories

def write_bonesis_model(
    trajectories: list
) -> None:
    stable_states = list()
    for trajectory in trajectories:
        if len(trajectory) == 1:
            continue
        bo_trajectory = str()
        fp = trajectory[-1]
        for config in trajectory:
            if config != fp:
                bo_trajectory += f"~bo.obs('{config}') >= "
            else:
                bo_trajectory += f"bo.fixed(~bo.obs('{config}'))"
                stable_states.append(fp)
        print(bo_trajectory)
    for a, b in combinations(stable_states, 2):
        print(f"~bo.obs('{a}') != ~bo.obs('{b}')")
    return None

trajectories = read_trajectories(args.infile)
write_bonesis_model(trajectories)
