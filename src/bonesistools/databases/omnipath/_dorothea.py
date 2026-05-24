#!/usr/bin/env python

from pathlib import Path
import warnings
from typing import (
    Any,
    List,
    Optional,
    Union,
    cast,
)

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

import networkx as nx
import pandas as pd

from ..ncbi import OutputIdentifierType, GeneSynonyms

DorotheaWrapper = Literal["op", "get"]

GET_DOROTHEA_DECOUPLER_REQUIREMENT = "decoupler<2.0.0"
OP_DOROTHEA_DECOUPLER_REQUIREMENT = "decoupler>=2.0.0"


def _get_decoupler_version(decoupler_module: Any) -> str:
    return getattr(decoupler_module, "__version__", "unknown")


def _raise_missing_dorothea_wrapper(
    wrapper: DorotheaWrapper,
    decoupler_module: Any,
) -> None:
    decoupler_version = _get_decoupler_version(decoupler_module)

    if wrapper == "get":
        raise AttributeError(
            "`decoupler.get_dorothea` is not available in the installed "
            f"decoupler version ({decoupler_version}). `wrapper='get'` "
            f"requires {GET_DOROTHEA_DECOUPLER_REQUIREMENT}, or an existing "
            "local cache file generated with a compatible decoupler version."
        )

    raise AttributeError(
        "`decoupler.op.dorothea` is not available in the installed decoupler "
        f"version ({decoupler_version}). `wrapper='op'` requires "
        f"{OP_DOROTHEA_DECOUPLER_REQUIREMENT}, or an existing local cache "
        "file generated with a compatible decoupler version."
    )


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
    this module. The complete selected decoupler resource is cached per wrapper
    and organism, then filtered locally according to `levels`.

    Parameters
    ----------
    organism: str or int (default: "mouse")
        Organism of interest. Accepted values depend on the selected decoupler
        wrapper.
    levels: list of str, optional
        DoRothEA confidence levels to retrieve. Valid levels are typically
        `"A"`, `"B"`, `"C"` and `"D"` for `decoupler.op.dorothea`; accepted
        values may differ for `decoupler.get_dorothea` and older decoupler
        versions. If None, use `["A", "B", "C"]`. The cache itself stores the
        complete selected decoupler resource; `levels` only filters the returned
        graph. With `wrapper="get"`, confidence levels are not passed to
        decoupler because the legacy wrapper is not level-parameterized; in
        practice this resource may contain only confidence level `"A"`.
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
        decoupler version. `wrapper="op"` requires decoupler>=2.0.0 when no
        cache is available, while `wrapper="get"` requires decoupler<2.0.0
        when no cache is available.
    Warns
    -----
    UserWarning
        If `wrapper="get"` is used with requested confidence levels other than
        `"A"`.
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

    if wrapper == "get" and any(level != "A" for level in levels):
        warnings.warn(
            "`wrapper='get'` uses the legacy `decoupler.get_dorothea` resource, "
            "which is not downloaded by confidence level and may only contain "
            "level 'A'; requested levels are applied only as a local filter.",
            UserWarning,
            stacklevel=2,
        )

    cache_dir = Path(__file__).resolve().parent / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_file = cache_dir / f"dorothea_{wrapper}_{organism}_complete.csv"

    if cache_file.exists() and not reload:
        dorothea_db = pd.read_csv(cache_file)

    else:
        import decoupler as _dc  # type: ignore

        dc = cast(Any, _dc)

        if wrapper == "op":
            if not hasattr(dc, "op") or not hasattr(dc.op, "dorothea"):
                _raise_missing_dorothea_wrapper(wrapper, dc)

            dorothea_db = dc.op.dorothea(
                organism=cast(str, organism),
                levels=["A", "B", "C", "D"],
                **kwargs,
            )

        else:
            if not hasattr(dc, "get_dorothea"):
                _raise_missing_dorothea_wrapper(wrapper, dc)

            dorothea_db = dc.get_dorothea(
                organism=cast(str, organism),
                **kwargs,
            )

        dorothea_db.to_csv(cache_file, index=False)

    dorothea_columns = list(dorothea_db.columns)
    dorothea_db = pd.DataFrame.from_records(
        (row for row in dorothea_db.to_dict("records") if row["confidence"] in levels),
        columns=dorothea_columns,
    )
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
