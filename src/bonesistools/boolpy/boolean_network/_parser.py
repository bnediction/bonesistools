#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from boolean import BooleanAlgebra

if TYPE_CHECKING:
    from ._network import BooleanNetwork, BooleanNetworkEnsemble


def read_bnet(
    file: Union[str, Path],
    ba: Optional[BooleanAlgebra] = None,
    check: bool = True,
) -> BooleanNetwork:
    """
    Read a Boolean network from a `.bnet` file.

    Parameters
    ----------
    file: str or Path
        Path to the `.bnet` file.
    ba: BooleanAlgebra (optional, default: None)
        Boolean algebra used to parse and store Boolean expressions. If `None`,
        a new BooleanAlgebra instance is created.
    check: bool (default: True)
        If `True`, validate that all symbols referenced by rules are defined as
        network components.

    Returns
    -------
    BooleanNetwork
        Parsed Boolean network.
    """

    from ._network import BooleanNetwork

    file = Path(file)

    rules = {}

    for line in file.read_text().splitlines():
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        component, rule = line.split(",", maxsplit=1)
        component = component.strip()
        rule = rule.strip().replace("!", "~")

        rules[component] = rule

    bn = BooleanNetwork(rules, ba=ba, check=False)

    if check:
        bn.validate()

    return bn


def read_bnet_directory(
    directory: Union[str, Path],
    ba: Optional[BooleanAlgebra] = None,
    recursive: bool = False,
    check: bool = True,
) -> BooleanNetworkEnsemble:
    """
    Read all `.bnet` files from a directory into a Boolean network ensemble.

    Parameters
    ----------
    directory: str or Path
        Directory containing `.bnet` files.
    ba: BooleanAlgebra (optional, default: None)
        Boolean algebra used to parse and store Boolean expressions. If `None`,
        a new BooleanAlgebra instance is created.
    recursive: bool (default: False)
        If `True`, recursively search subdirectories for `.bnet` files.
    check: bool (default: True)
        If `True`, validate that all symbols referenced by rules are defined as
        network components.

    Returns
    -------
    BooleanNetworkEnsemble
        Ensemble containing all loaded Boolean networks.

    Raises
    ------
    FileNotFoundError
        If `directory` does not exist.
    NotADirectoryError
        If `directory` is not a directory.
    ValueError
        If no `.bnet` file is found.
    """

    from ._network import BooleanNetworkEnsemble

    directory = Path(directory)

    if not directory.exists():
        raise FileNotFoundError(f"directory does not exist: {directory}")

    if not directory.is_dir():
        raise NotADirectoryError(f"path is not a directory: {directory}")

    pattern = "**/*.bnet" if recursive else "*.bnet"

    files = sorted(directory.glob(pattern))

    if len(files) == 0:
        raise ValueError(f"no '.bnet' file found in directory: {directory}")

    bns = [read_bnet(file, ba=ba, check=check) for file in files]

    return BooleanNetworkEnsemble(bns=bns)
