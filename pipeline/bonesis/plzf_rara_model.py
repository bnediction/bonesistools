#!/usr/bin/env python

from pathlib import Path
from utils.genesyn import GeneSynonyms

import pandas as pd

import networkx as nx
from bonesis import BoNesis

important_nodes = []

def bomodel(bo: BoNesis, file: Path = Path("plzf_rara_model.txt")) -> None:
    with open(file, "r") as file:
        for line in file:
            eval(line)
    return None

def load_bin(file: Path, gene_synonyms: GeneSynonyms = None) -> dict:
    meta_bin = pd.read_csv(file, index_col=0)
    if gene_synonyms is not None and isinstance(gene_synonyms, GeneSynonyms):
        gene_synonyms.df_standardization(meta_bin, axis=1, copy=False)
    return {config: genes.dropna().to_dict() for config, genes in meta_bin.iterrows()}

def collectri_to_grn(
    collectri: pd.DataFrame,
    sign_label: str = "weight",
    remove_pmid: bool = False
    ) -> nx.MultiDiGraph:
    if sign_label is not None:
        collectri = collectri.rename(columns = {sign_label:"sign"})
    if remove_pmid is True:
        remove_pmid = "PMID" in collectri.columns
    return nx.from_pandas_edgelist(
        df = collectri.drop("PMID", axis=1) if remove_pmid else collectri,
        source="source",
        target="target",
        edge_attr=True,
        create_using=nx.MultiDiGraph
    )
