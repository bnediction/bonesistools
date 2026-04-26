#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Mapping,
    List,
    Any
)

import networkx as nx

from ..ncbi import (
    OutputIdentifierType,
    GeneSynonyms
)

def load_dorothea_grn(
    organism: Union[str, int] = "mouse",
    levels: List[str] = ["A", "B", "C"],
    genesyn: Optional[GeneSynonyms] = None,
    gene_identifier_type: OutputIdentifierType = "official_name",
    **kwargs: Mapping[str, Any]
)-> nx.MultiDiGraph:
    """
    Provide a Graph Regulatory Network (GRN) derived from DoRothEA database [1].

    Parameters
    ----------
    organism: str | int (default: 'mouse')
        Common name or identifier of the organism of interest.
        Identifier can be NCBI ID, EnsemblID or latin name.
    levels: List[str] (default: ['A', 'B', 'C'])
        List of confidence levels in regulation.
        Confidence score range are between A and D, A being the most confident and D being the less confident.
    gene_synonyms: GeneSynonyms (optional, default: None)
        If GeneSynonyms object is passed, then gene identifiers are converted into the desired identifier format.
    input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
        Gene identifier input format.
    output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
        Gene identifier output format.
    **kwargs: Mapping[str, Any]
        Keyword-arguments passed to function 'omnipath.interactions.Dorothea.get'.
    
    Returns
    -------
    Return graph from DoRothEA database.

    References
    ----------
    [1] Garcia-Alonso et al. (2019). Benchmark and integration of resources
    for the estimation of human transcription factor activities.
    Genome research 29(8), 1363-1375 (https://genome.cshlp.org/content/29/8/1363)
    """

    if not isinstance(organism, (str, int)):
        raise TypeError(f"unsupported argument type for 'organism': expected {str} or {int} but received {type(organism)}")
    
    import decoupler as dc # type: ignore
    
    try:
        dorothea_db = dc.get_dorothea(organism=organism, levels=levels, **kwargs)
    except:
        dorothea_db = dc.op.dorothea(organism=organism, levels=levels, **kwargs)
    dorothea_db = dorothea_db.rename(columns = {"weight":"sign"})
    dorothea_db["sign"] = dorothea_db["sign"].apply(lambda x: -1 if x < 0 else 1)

    grn = nx.from_pandas_edgelist(
        df = dorothea_db,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph
    )
    if genesyn is None:
        return grn
    elif isinstance(genesyn, GeneSynonyms):
        genesyn(
            grn,
            input_identifier_type="name",
            output_identifier_type=gene_identifier_type,
            copy=False
        )
        return grn
    else:
        raise TypeError(f"unsupported argument type for 'gene_synonyms': expected {GeneSynonyms} but received {type(genesyn)}")
