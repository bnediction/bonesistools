#!/usr/bin/env python

"""
Parsers for Boolean algebra objects.
"""

import json
from pathlib import Path
from typing import Dict, Optional, Union

from ..._typing import FileOrientation
from ..._validation import _as_orientation
from ..._warnings import _warn_deprecated_argument
from ._hypercube import Hypercube


def read_hypercube(file: Union[str, Path]) -> Hypercube:
    """
    Read a hypercube from a JSON file.

    Examples
    --------
    A valid JSON file has the following form:

    {
        "Gata1": 1,
        "Spi1": 0,
        "Tal1": "*"
    }

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
    """

    file = Path(file)

    with open(file) as fp:
        mapping = json.load(fp)

    return Hypercube(mapping)


def read_hypercubes(
    file: Union[str, Path],
    orientation: FileOrientation = "columns",
    sep: str = ",",
    *,
    axis: Optional[FileOrientation] = None,
) -> Dict[str, Hypercube]:
    """
    Read named hypercubes from a CSV, TSV or JSON file.

    For tabular files, missing values are interpreted as free values `"*"`.

    JSON files must store a mapping from hypercube names to component mappings.
    Missing JSON values (`null`) are interpreted as free values `"*"`.

    Parameters
    ----------
    file: str or Path
        Input file.
    orientation: {"columns", "rows"} (default: "columns")
        File orientation. If `"columns"`, each column stores one hypercube. If
        `"rows"`, each row stores one hypercube. Only used for tabular files.
    sep: str (default: ",")
        Field delimiter used for tabular files.
    axis: {"columns", "rows"}, optional
        Deprecated alias for `orientation`.

    Returns
    -------
    Dict[str, Hypercube]
        Mapping from hypercube names to parsed hypercubes.

    Raises
    ------
    ValueError
        If the file extension or orientation value is not supported, or if one
        of the stored values cannot be interpreted as a PartialBoolean value.
    """

    if axis is not None:
        _warn_deprecated_argument("axis", "orientation", stacklevel=2)
        orientation = axis

    orientation = _as_orientation(orientation)
    file = Path(file)
    suffix = file.suffix.lower()

    if suffix == ".json":

        with open(file) as fp:
            mapping = json.load(fp)

        return {
            name: Hypercube(
                {
                    component: "*" if value is None else value
                    for component, value in hypercube.items()
                }
            )
            for name, hypercube in mapping.items()
        }

    if suffix == ".tsv":
        sep = "\t"

    elif suffix != ".csv":
        raise ValueError(
            f"unsupported input format: '{file.suffix}' "
            "(supported formats: .csv, .tsv, .json)"
        )

    import pandas as pd

    data = pd.read_csv(file, index_col=0, sep=sep)
    data = data.where(data.notna(), "*")

    if orientation == "rows":
        data = data.T

    hypercubes = {}

    for name in data.columns:
        mapping = {str(gene): value for gene, value in data[name].to_dict().items()}
        hypercubes[str(name)] = Hypercube(mapping)

    return hypercubes
