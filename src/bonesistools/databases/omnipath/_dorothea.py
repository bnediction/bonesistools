#!/usr/bin/env python

from pathlib import Path
from typing import Any, List, Optional, Union

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

import networkx as nx
import pandas as pd

from ..ncbi import OutputIdentifierType, GeneSynonyms

DorotheaWrapper = Literal["op", "get"]


def load_dorothea_grn(
    organism: Union[str, int] = "mouse",
    levels: Optional[List[str]] = None,
    genesyn: Optional[GeneSynonyms] = None,
    gene_identifier_type: OutputIdentifierType = "official_name",
    wrapper: DorotheaWrapper = "op",
    reload: bool = False,
    **kwargs: Any,
) -> nx.MultiDiGraph:
    """
    Load a transcriptional regulatory network derived from DoRothEA.

    The DoRothEA resource is retrieved through `decoupler` and converted into a
    signed NetworkX MultiDiGraph. Two decoupler wrappers can be used because
    current and historical decoupler APIs may expose different DoRothEA
    retrieval functions.

    By default, this function uses `decoupler.op.dorothea`, which corresponds
    to the current OmniPath-oriented decoupler interface. The alternative
    `decoupler.get_dorothea` wrapper can be selected explicitly for
    compatibility with historical analyses or environments.

    Retrieved resources are cached locally in an `.cache` directory next to
    this module. The cache is separated by wrapper, organism and confidence
    levels, so that results obtained through different decoupler wrappers are
    never mixed.

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Accepted values depend on the selected decoupler
        wrapper.
    levels: list of str, optional
        DoRothEA confidence levels to retrieve. Valid levels are typically
        `"A"`, `"B"`, `"C"` and `"D"` for `decoupler.op.dorothea`; accepted
        values may differ for `decoupler.get_dorothea` and older decoupler
        versions. If None, use `["A", "B", "C"]`.
    genesyn: GeneSynonyms, optional
        GeneSynonyms object used to convert graph node identifiers.
    gene_identifier_type: OutputIdentifierType (default: "official_name")
        Output gene identifier type used when `genesyn` is provided.
    wrapper: {"op", "get"} (default: "op")
        Decoupler wrapper used to retrieve DoRothEA.

        If `"op"`, use `decoupler.op.dorothea`.

        If `"get"`, use `decoupler.get_dorothea`.

        No automatic fallback is performed between wrappers, because the two
        wrappers may retrieve different resources or apply different
        preprocessing pipelines.
    reload: bool (default: False)
        If True, ignore the local cache and retrieve the resource again through
        decoupler.
    **kwargs: Any
        Keyword arguments passed to the selected decoupler DoRothEA wrapper.

    Returns
    -------
    nx.MultiDiGraph
        Signed regulatory network. Edges contain the attributes returned by
        decoupler, with `weight` renamed to `sign` and converted to -1 or 1.

    Raises
    ------
    TypeError
        If `organism` or `genesyn` has an unsupported type.
    ValueError
        If `wrapper` is not `"op"` or `"get"`.
    AttributeError
        If the selected decoupler wrapper is not available in the installed
        decoupler version.
    """

    if levels is None:
        levels = ["A", "B", "C"]

    if not isinstance(organism, (str, int)):
        raise TypeError(
            f"unsupported argument type for 'organism': "
            f"expected {str} or {int} but received {type(organism)}"
        )

    if wrapper not in ["op", "get"]:
        raise ValueError(
            f"invalid argument value for 'wrapper': "
            f"expected 'op' or 'get' but received {wrapper!r}"
        )

    cache_dir = Path(__file__).resolve().parent / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    levels_key = "".join(levels)
    cache_file = cache_dir / f"dorothea_{wrapper}_{organism}_{levels_key}.csv"

    if cache_file.exists() and not reload:
        dorothea_db = pd.read_csv(cache_file)

    else:
        import decoupler as dc  # type: ignore

        if wrapper == "op":
            dorothea_db = dc.op.dorothea(
                organism=organism,
                levels=levels,
                **kwargs,
            )

        else:
            dorothea_db = dc.get_dorothea(
                organism=organism,
                levels=levels,
                **kwargs,
            )

        dorothea_db.to_csv(cache_file, index=False)

    dorothea_db = dorothea_db.rename(columns={"weight": "sign"})
    dorothea_db["sign"] = dorothea_db["sign"].apply(lambda x: -1 if x < 0 else 1)

    grn = nx.from_pandas_edgelist(
        df=dorothea_db,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph,
    )

    if genesyn is None:
        return grn

    if isinstance(genesyn, GeneSynonyms):
        genesyn(
            grn,
            input_identifier_type="name",
            output_identifier_type=gene_identifier_type,
            copy=False,
        )
        return grn

    raise TypeError(
        f"unsupported argument type for 'genesyn': "
        f"expected {GeneSynonyms} but received {type(genesyn)}"
    )
