#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Mapping,
    List,
    Any
)

import networkx as nx

from ..ncbi import GeneSynonyms

def load_dorothea_grn(
    organism: Union[str, int] = "mouse",
    levels: List[str] = ["A", "B", "C"],
    genesyn: Optional[GeneSynonyms] = None,
    gene_type: str = "genename",
    alias_type: str = "referencename",
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
    gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
        Gene identifier input format.
    alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
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
    
    dorothea_db = dc.get_dorothea(organism=organism, levels=levels, **kwargs)
    dorothea_db = dorothea_db.rename(columns = {"weight":"sign"})

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
            gene_type=gene_type,
            alias_type=alias_type,
            copy=False
        )
        return grn
    else:
        raise TypeError(f"unsupported argument type for 'gene_synonyms': expected {GeneSynonyms} but received {type(genesyn)}")
