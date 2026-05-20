#!/usr/bin/env python

"""
Parsers for Boolean algebra objects.
"""

import json
from pathlib import Path
from typing import Union

from ._hypercube import Hypercube


def read_hypercube(file: Union[str, Path]) -> Hypercube:
    """
    Read a hypercube from a JSON file.

    Parameters
    ----------
    file: str or Path
        JSON file storing a mapping from component names to partial Boolean
        values.

    Returns
    -------
    Hypercube
        Parsed hypercube.

    Raises
    ------
    ValueError
        If one of the stored values cannot be interpreted as a
        PartialBoolean value.

    Examples
    --------
    A valid JSON file has the following form:

    {
        "Gata1": 1,
        "Spi1": 0,
        "Tal1": "*"
    }
    """

    file = Path(file)

    with open(file) as fp:
        mapping = json.load(fp)

    return Hypercube(mapping)
