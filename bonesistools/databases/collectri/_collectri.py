#!/usr/bin/env python

from typing import Optional, Union

import networkx as nx

from ..ncbi import GeneSynonyms

def load_grn(
    organism: Union[str, int] = "human",
    split_complexes: bool = False,
    remove_pmid: bool = False,
    gene_synonyms: Optional[GeneSynonyms] = None,
    input_type: str = "genename",
    output_type: str = "referencename",
    **kwargs
)-> nx.MultiDiGraph:
    """
    Provide a Graph Regulatory Network (GRN) derived from Collectri database.

    Parameters
    ----------
    organism
        Common name or identifier of the organism of interest (default: human).
        Identifier can be NCBI ID, EnsemblID or latin name.
    split_complexes
        Specify whether to split complexes into subunits.
    remove_pmid
        Specify whether to remove PMIDs in node labels.
    kwargs
        Keyword-arguments passed to omnipath.interactions.CollecTRI.get().
    
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
    
    import decoupler as dc
    
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
    if gene_synonyms is None:
        return grn
    elif isinstance(gene_synonyms, GeneSynonyms):
        gene_synonyms.graph_standardization(
            graph=grn,
            input_type=input_type,
            output_type=output_type,
            copy=False
        )
        return grn
    else:
        raise TypeError(f"unsupported argument type for 'gene_synonyms': expected {GeneSynonyms} but received {type(gene_synonyms)}")

