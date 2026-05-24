#!/usr/bin/env python

"""
Parsers for Boolean algebra objects.
"""

import json
from pathlib import Path
from typing import Union, Dict

from ..._compat import Literal

from ._hypercube import Hypercube

Axis = Literal["columns", "rows"]


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
    axis: Axis = "columns",
    sep: str = ",",
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
    axis: str (default: "columns")
        Axis along which hypercubes are stored. Supported values are
        `"columns"` and `"rows"`. Only used for tabular files.
    sep: str (default: ",")
        Field delimiter used for tabular files.

    Returns
    -------
    Dict[str, Hypercube]
        Mapping from hypercube names to parsed hypercubes.

    Raises
    ------
    ValueError
        If the file extension or axis value is not supported, or if one of the
        stored values cannot be interpreted as a PartialBoolean value.
    """

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

    if axis not in {"columns", "rows"}:
        raise ValueError(
            f"invalid argument value for 'axis': "
            f"expected 'columns' or 'rows' but received {axis!r}"
        )

    import pandas as pd

    data = pd.read_csv(file, index_col=0, sep=sep)
    data = data.where(data.notna(), "*")

    if axis == "rows":
        data = data.T

    return {name: Hypercube(data[name].to_dict()) for name in data.columns}
