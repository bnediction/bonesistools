#!/usr/bin/env python

from typing import (
    Optional,
    Union,
    Mapping,
    Any
)

import networkx as nx

from ..ncbi import GeneSynonyms

def load_collectri_grn(
    organism: Union[str, int] = "mouse",
    split_complexes: bool = False,
    remove_pmid: bool = False,
    genesyn: Optional[GeneSynonyms] = None,
    gene_type: str = "genename",
    alias_type: str = "referencename",
    **kwargs: Mapping[str, Any]
)-> nx.MultiDiGraph:
    """
    Provide a Graph Regulatory Network (GRN) derived from Collectri database [1].

    Parameters
    ----------
    organism: str | int (default: 'mouse')
        Common name or identifier of the organism of interest.
        Identifier can be NCBI ID, EnsemblID or latin name.
    split_complexes: bool (default: False)
        Specify whether to split complexes into subunits.
    remove_pmid: bool (default: False)
        Specify whether to remove PMIDs in node labels.
    gene_synonyms: GeneSynonyms (optional, default: None)
        If GeneSynonyms object is passed, then gene identifiers are converted into the desired identifier format.
    gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
        Gene identifier input format.
    alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
        Gene identifier output format.
    **kwargs: Mapping[str, Any]
        Keyword-arguments passed to function 'omnipath.interactions.CollecTRI.get'.
    
    Returns
    -------
    Return graph from Collectri database.

    References
    ----------
    [1] Müller-Dott et al. (2023). Expanding the coverage of regulons from high-confidence
    prior knowledge for accurate estimation of transcription factor activities.
    Nucleic Acids Research, 51(20), 10934-10949 (https://doi.org/10.1093/nar/gkad841)
    """

    if not isinstance(organism, (str, int)):
        raise TypeError(f"unsupported argument type for 'organism': expected {str} or {int} but received {type(organism)}")
    if not isinstance(split_complexes, bool):
        raise TypeError(f"unsupported argument type for 'split_complexes': expected {bool} but received {type(split_complexes)}")
    
    import decoupler as dc # type: ignore
    
    collectri_db = dc.get_collectri(organism=organism, split_complexes=split_complexes, **kwargs)
    collectri_db = collectri_db.rename(columns = {"weight":"sign"})
    if isinstance(remove_pmid, bool):
        if remove_pmid:
            collectri_db = collectri_db.drop("PMID", axis=1)
        else:
            pass
    else:
        raise TypeError(f"unsupported argument type for 'remove_pmid': expected {bool} but received {type(remove_pmid)}")
    grn = nx.from_pandas_edgelist(
        df = collectri_db,
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

