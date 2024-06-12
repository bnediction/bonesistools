#!/usr/bin/env python

from typing import Optional, Union

import decoupler as dc

import networkx as nx

from .genesyn import GeneSynonyms

def load_grn(
    organism: Union[str,int]="human",
    split_complexes=False,
    remove_pmid: bool=False,
    gene_synonyms: Optional[GeneSynonyms] = None,
    in_alias_type: str="genename",
    out_alias_type: str="referencename",
    **kwargs
)-> nx.MultiDiGraph:
    """
    Provide a Graph Regulatory Network (GRN) derived from Collectri database
    (<Müller-Dott et al. (2023), Expanding the coverage of regulons from high-confidence
    prior knowledge for accurate estimation of transcription factor activities>)

    Parameters
    ----------
    organism
        common name or identifier of the organism of interest (default: human).
        Identifier can be NCBI ID, EnsemblID or latin name.
    split_complexes
        specify whether to split complexes into subunits
    remove_pmid
        specify whether to remove PMIDs in node labels
    kwargs
        keyword-arguments passed to`omnipath.interactions.CollecTRI.get    
        
    Returns
    -------
    Return graph from Collectri database.
    """

    if not isinstance(organism, (str, int)):
        raise TypeError(f"parameter `organism` is not a string or integer instance")
    if not isinstance(split_complexes, bool):
        raise TypeError(f"parameter `split_complexes` is not a boolean instance")
    collectri_db = dc.get_collectri(organism=organism, split_complexes=split_complexes, **kwargs)
    collectri_db = collectri_db.rename(columns = {"weight":"sign"})
    if isinstance(remove_pmid, bool):
        if remove_pmid:
            collectri_db = collectri_db.drop("PMID", axis=1)
        else:
            pass
    else:
        raise ValueError(f"parameter `remove_pmid` is not a boolean instance")
    grn = nx.from_pandas_edgelist(
        df = collectri_db,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph
    )
    if gene_synonyms is None:
        return grn
    elif isinstance(gene_synonyms, GeneSynonyms):
        gene_synonyms.graph_standardization(
            graph=grn,
            in_alias_type=in_alias_type,
            out_alias_type=out_alias_type,
            copy=False
        )
        return grn
    else:
        raise TypeError(f"parameter `gene_synonyms` is not a {GeneSynonyms} instance")

