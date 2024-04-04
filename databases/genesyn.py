#!/usr/bin/env python

import ctypes

import warnings

from typing import Union, Any, Sequence, Dict, List, Tuple
from pandas._typing import Axis
try:
    from collections import Sequence as SequenceInstance
except:
    from collections.abc import Sequence as SequenceInstance

import os
from pathlib import Path

from pandas import DataFrame

import networkx as nx
from networkx import Graph

class GeneSynonyms(object):

    def __init__(self, ncbi_file: Path = None) -> None:
        if isinstance(ncbi_file, Path):
            self.ncbi_file = ncbi_file
        elif ncbi_file is None:
            self.ncbi_file = Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/.mus_musculus_gene_info.tsv")
        else:
            raise TypeError(f"`ncbi_file` is not a Path")
        self.gene_synonyms = self.__synonyms_from_NCBI(self.ncbi_file)
        self.__upper_gene_synonyms = {gene.upper(): self.gene_synonyms[gene] for gene in self.gene_synonyms.keys()}
        return None
    
    def __get__(self, attribute: str = None) -> Any:
        if attribute is None or attribute == "gene_synonyms":
            return {gene: ncbi_reference_name.value.decode() for gene, ncbi_reference_name in self.gene_synonyms.items()}
        elif attribute == "__upper_gene_synonyms":
            raise AttributeError(f"`__upper_gene_synonyms` attribute is private")
        else:
            return getattr(self, attribute)
    
    def __set__(self, ncbi_file: Path) -> None:
        if isinstance(ncbi_file, Path):
            self.ncbi_file = ncbi_file
        else:
            raise TypeError(f"`ncbi_file` is not a Path")
        self.gene_synonyms = self.__synonyms_from_NCBI(self.gi_file)
        self.__upper_gene_synonyms = {gene.upper(): self.gene_synonyms[gene] for gene in self.gene_synonyms.keys()}
        return None
    
    def __call__(self, data: Union[Sequence[Tuple[str, str, Dict[str, int]]], DataFrame, Graph], *args, **kwargs):
        if (isinstance(data, SequenceInstance) and not isinstance(data, str)) or isinstance(data, set):
            return self.sequence_standardization(data, *args, **kwargs)
        elif isinstance(data, DataFrame):
            return self.df_standardization(data, *args, **kwargs)
        elif isinstance(data, Graph):
            return self.graph_standardization(data, *args, **kwargs)
        else:
            raise TypeError(f"fail to convert gene name: `data` has incorrect type")
   
    def __synonyms_from_NCBI(self, gi_file: Path) -> dict:
        """
        Create a dictionary matching each gene name to its NCBI reference gene name.
        For speeding up the task facing a large matrix from NCBI, the parsing of the NCBI gene data is run with awk.

        Parameters
        ----------
        gi_file
            path to the NCBI gene info data

        Returns
        -------
        Return a dictionary where keys correspond to gene name and values correspond to reference gene name
        """

        gi_file_cut = Path(f"{gi_file}_cut")
        command_parsing = "awk -F'\t' '{print $3 \"\t\" $5 \"\t\" $11}' " + str(gi_file) + " | tr \| '\t' > " + str(gi_file_cut) + " ; sed -i 1d " + str(gi_file_cut)
        os.system(command_parsing)

        gene_synonyms_dict = dict()
        reference_names = set()

        with open (gi_file_cut, "r") as file_synonyms:
            for gene in file_synonyms:
                gene = gene.strip()
                gene_synonyms_list = gene.split("\t")
                ncbi_reference_name = gene_synonyms_list.pop(0)
                res = [_synonym for _synonym in gene_synonyms_list if (_synonym != "-" and _synonym != ncbi_reference_name)]

                gene_synonyms_dict[ncbi_reference_name] = ctypes.create_string_buffer(ncbi_reference_name.encode())
                reference_names.add(ncbi_reference_name)

                for gene in res:
                    if gene not in reference_names and gene not in gene_synonyms_dict:
                        # warning with NCBI list of synonyms: a noun can be the synonym of several reference names. Arbitrary, the choosen one is the first.
                        gene_synonyms_dict[gene] = ctypes.create_string_buffer(ncbi_reference_name.encode())

        os.system(f"rm {str(gi_file_cut)}")
        return gene_synonyms_dict
    
    def get_reference_gene_name(self, gene_name: str) -> str:
        """
        Provide the reference name with respect to a gene name.

        Parameters
        ----------
        gene_name
            name of a gene
        gene_synonyms_dict
            dictionary where keys correspond to gene name and values correspond to reference gene name

        Returns
        -------
        Given a gene name, return its reference name.
        """

        _gene_name = gene_name.upper()
        if _gene_name in self.__upper_gene_synonyms:
            return self.__upper_gene_synonyms[_gene_name].value.decode()
        else:
            warnings.warn(f"NCBI does not find a correspondance for {gene_name}.", stacklevel=10)
            return gene_name
    
    def sequence_standardization(self, gene_sequence: Sequence[str]) -> Sequence[str]:
        """
        Create a copy of the input Sequence, with each gene name replaced by its reference name.

        Parameters
        ----------
        interaction_list
            list of tuples containing string (source) + string (target) + dict (sign = -1 or 1)
        
        Returns
        -------
        return an interaction list where each gene name is converted into its reference value.
        """
        
        standardized_gene_sequence = list()
        for gene in gene_sequence:
            standardized_gene_sequence.append(self.get_reference_gene_name(gene))
        
        standardized_gene_sequence = type(gene_sequence)(standardized_gene_sequence)
        
        return standardized_gene_sequence

    def interaction_list_standardization(self, interactions_list: Sequence[Tuple[str, str, Dict[str, int]]]) -> List[Tuple[str, str, Dict[str, int]]]:
        """
        Create a copy of the input list of pairwise interactions, with each gene name replaced by its reference name.

        Parameters
        ----------
        interaction_list
            list of tuples containing string (source) + string (target) + dict (sign = -1 or 1)

        Returns
        -------
        return an interaction list where each gene name is converted into its reference value.
        """

        standardized_interactions_list = list()
        for interaction in interactions_list:
            source = self.get_reference_gene_name(interaction[0])
            target = self.get_reference_gene_name(interaction[1])
            standardized_interactions_list.append((source, target, interaction[2]))

        return standardized_interactions_list

    def df_standardization(
        self,
        df: DataFrame,
        axis: Axis = 0,
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Replace gene name with its reference gene name into `df`.

        Parameters
        ----------
        df
            dataframe where names must be standardized
        axis
            whether to rename labels from the index (0 or `index`) or columns (1 or `columns`)
        copy
            return a copy instead of updating `df`
        
        Returns
        -------
        Depending on `inplace`, update or return dataframe with standardized gene name.
        """

        df = df.copy() if copy is True else df

        synonyms = list()

        if axis == 0 or axis == "index":
            gene_iterator = iter(df.index)
        elif axis == 1 or axis == "columns":
            gene_iterator = iter(df.columns)
        else:
            raise ValueError(f"No axis named {axis} for object type DataFrame")
        for gene in gene_iterator:
            ncbi_reference_name = self.get_reference_gene_name(gene)
            synonyms.append(ncbi_reference_name)
        if axis == 0 or axis == "index":
            df.index = synonyms
        elif axis == 1 or axis == "columns":
            df.columns = synonyms

        if copy is True:
            return df
    
    def graph_standardization(
        self,
        graph: Graph,
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene name with its reference gene name into `graph`.

        Parameters
        ----------
        graph
            graph where nodes must be standardized
        copy
            return a copy instead of updating `graph`
        
        Returns
        -------
        Depending on `inplace`, update or return graph with standardized gene name.
        """

        synonym_mapping = dict()
        for gene in graph.nodes:
            synonym_mapping[gene] = self.get_reference_gene_name(gene)
        if copy is True:
            return nx.relabel_nodes(graph, mapping=synonym_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=synonym_mapping, copy=False)
            return None
