#!/usr/bin/env python

import ctypes
from collections import namedtuple
try:
    from sortedcontainers import SortedSet
except:
    pass

import warnings

from typing import Union, Any, Sequence, Dict, List, Tuple
from pandas._typing import Axis
try:
    from collections import Sequence as SequenceInstance
except:
    from collections.abc import Sequence as SequenceInstance

import os
from pathlib import Path

import re

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
        self.gene_alias_mapping = self.__alias_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_alias_mapping["genename"].items()}
        return None
    
    def __get__(self, attribute: str = "gene_alias_mapping") -> Any:
        if attribute == "__upper_gene_names_mapping":
            raise AttributeError(f"private attribute")
        else:
            return getattr(self, attribute)
    
    def __set__(self, ncbi_file: Path) -> None:
        if isinstance(ncbi_file, Path):
            self.ncbi_file = ncbi_file
        else:
            raise TypeError(f"`ncbi_file` is not a Path")
        self.gene_alias_mapping = self.__alias_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_alias_mapping["genename"].items()}
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

    def __alias_from_NCBI(self, gi_file: Path) -> dict:
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
        command_parsing = "awk -F'\t' 'NR>1 {print $2 \"\t\" $3 \"\t\" $5 \"\t\" $11 \"\t\" $6}' " + str(gi_file) + " > " + str(gi_file_cut)
        os.system(command_parsing)

        gene_alias_mapping = {
            "geneid": dict(),
            "genename": dict(),
            "ensemblid": dict(),
            "mgi": dict()
        }

        try:
            gene_names = SortedSet()
        except:
            gene_names = set()

        with open(gi_file_cut, "r") as file:
            for line in file:
                geneinfo = line.strip().split("\t")
                _geneid = geneinfo.pop(0)
                _reference_name = geneinfo.pop(0)
                _synonyms = [_synonym for _synonym in geneinfo.pop(0).split("|") + geneinfo.pop(0).split("|") \
                    if (_synonym != "-" and _synonym != _reference_name)]
                geneinfo = geneinfo[0]
                _mgi = re.findall("MGI:MGI:[0-9]*", geneinfo)
                _mgi = re.sub("^MGI:MGI:", "", _mgi[0]) if _mgi else None
                _ensemblid = re.findall("Ensembl:[A-Z]{7}[0-9]{11}", geneinfo)
                _ensemblid = re.sub("^Ensembl:", "", _ensemblid[0]) if _ensemblid else None

                gene_names.add(_reference_name)
                pointer_to_geneid = ctypes.create_string_buffer(_geneid.encode())
                if _geneid:
                    gene_alias_mapping["geneid"][_geneid] = namedtuple("GeneAlias", ["reference_genename", "ensemblid", "mgi"])(_reference_name, _ensemblid, _mgi)
                gene_alias_mapping["genename"][_reference_name] = pointer_to_geneid
                for _synonym in _synonyms:
                    if _synonym not in gene_names:
                        gene_names.add(_synonym)
                        gene_alias_mapping["genename"][_synonym] = pointer_to_geneid
                    else:
                        print(f"synonym {_synonym} multiple times with ref name: {_reference_name}")
                if _ensemblid:
                    gene_alias_mapping["ensemblid"][_ensemblid] = pointer_to_geneid
                if _mgi:
                    gene_alias_mapping["mgi"][_mgi] = pointer_to_geneid

        os.system(f"rm {str(gi_file_cut)}")
        return gene_alias_mapping
    
    def get_geneid(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the geneid with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|ensemblid|mgi

        Returns
        -------
        Given an alias, return its geneid
        """

        if alias_type == "genename":
            gene_name = alias.upper()
            if gene_name in self.__upper_gene_names_mapping:
                return self.__upper_gene_names_mapping[gene_name].value.decode()
            else:
                warnings.warn(f"no correspondance for gene name '{alias}'", stacklevel=10)
                return None
        elif alias_type in ["ensemblid", "mgi"]:
            if alias in self.gene_alias_mapping[alias_type]:
                return self.gene_alias_mapping[alias_type][alias].value.decode()
            else:
                warnings.warn(f"no geneid correspondance for {alias_type} '{alias}'", stacklevel=10)
                return None
        else:
            raise ValueError("Invalid value for argument 'alias_type'")
    
    def get_reference_name(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the reference name with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|geneid|ensemblid|mgi

        Returns
        -------
        Given an alias, return its reference name
        """
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_alias_mapping["geneid"]:
            return self.gene_alias_mapping["geneid"][geneid].reference_genename
        else:
            warnings.warn(f"no reference genename correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None
    
    def get_mgi(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the mgi with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|geneid|ensemblid

        Returns
        -------
        Given an alias, return its mgi
        """
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_alias_mapping["geneid"]:
            return self.gene_alias_mapping["geneid"][geneid].mgi
        else:
            warnings.warn(f"no mgi correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None

    def get_ensemblid(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the ensemblid with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|geneid|mgi

        Returns
        -------
        Given an alias, return its ensemblid
        """
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_alias_mapping["geneid"]:
            return self.gene_alias_mapping["geneid"][geneid].ensemblid
        else:
            warnings.warn(f"no ensemblid correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None

    def __convert(out_alias_type: str="referencename") -> str:
        """
        Convert gene alias.

        Parameters
        ----------
        out_alias_type
            geneid|referencename|ensemblid|mgi

        Returns
        -------
        Return a function converting gene labels.
        """

        if out_alias_type == "referencename":
            out_alias_type = "reference_gene_name"
        return eval(f"self.get_{out_alias_type}")

    def sequence_standardization(
            self,
            gene_sequence: Sequence[str],
            in_alias_type: str="genename",
            out_alias_type: str="referencename"
        ) -> Sequence[str]:
        """
        Create a copy of the input Sequence, with corresponding alias.

        Parameters
        ----------
        gene_sequence
            list of tuples containing string (source) + string (target) + dict (sign = -1 or 1)
        in_alias_type
            genename|geneid|ensemblid|mgi
        out_alias_type
            referencename|geneid|ensemblid|mgi
        
        Returns
        -------
        return a gene sequence where each gene alias is converted into the user-defined alias type.
        """
        
        convert = self.__convert(out_alias_type)
        standardized_gene_sequence = list()
        for gene in gene_sequence:
            standardized_gene_sequence.append(convert(gene,in_alias_type))
        
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
