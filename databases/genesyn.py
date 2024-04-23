#!/usr/bin/env python

import ctypes
from functools import partial
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

ORGANISMS = {"mouse", "human", "escherichia coli"}

FTP_GENE_INFO = {
    "mouse": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    "human": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "escherichia coli": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz"
}

NCBI_FILES = {
    "mouse": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/.ncbi_gi/.mus_musculus_gene_info.tsv"),
    "human": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/.ncbi_gi/.homo_sapiens_gene_info.tsv"),
    "escherichia coli": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/.ncbi_gi/.escherichia_coli_gene_info.tsv")
}

class GeneSynonyms(object):

    def __init__(self, organism: str="mouse", force_download: bool=False) -> None:
        organism = organism.lower().replace("-"," ")
        if organism in ORGANISMS:
            self.organism = organism
        else:
            raise ValueError("invalid value for argument 'organism'")
        self.ncbi_file = NCBI_FILES[organism]
        if force_download is True:
            command_parsing = f"wget --quiet --show-progress -cO {self.ncbi_file}.gz {FTP_GENE_INFO[self.organism]} && gunzip --quiet {self.ncbi_file}.gz"
            os.system(command_parsing)
        elif force_download is False:
            pass
        else:
            raise ValueError(f"invalid value for argument 'force_download' (expected {type(bool)}, got {type(force_download)})")
        self.force_download = force_download
        self.gene_aliases_mapping = self.__aliases_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_aliases_mapping["genename"].items()}
        self.databases = {_db for _db in self.gene_aliases_mapping["databases"].keys()}
        return None
    
    def __get__(self, attribute: str = "gene_aliases_mapping") -> Any:
        if attribute == "__upper_gene_names_mapping":
            raise AttributeError(f"private attribute")
        else:
            return getattr(self, attribute)
    
    def __set__(self, organism: str=None, force_download: bool=False) -> None:
        if organism is None and self.organism in ORGANISMS:
            pass
        elif organism in ORGANISMS:
            self.organism = organism
        else:
            raise ValueError("invalid value for argument 'organism'")
        self.ncbi_file = NCBI_FILES[self.organism]
        if force_download is True:
            command_parsing = f"wget --quiet --show-progress -cO {self.ncbi_file}.gz {FTP_GENE_INFO[self.organism]} && gunzip --quiet {self.ncbi_file}.gz"
            os.system(command_parsing)
        elif force_download is False:
            pass
        else:
            raise ValueError(f"invalid value for argument 'force_download' (expected {type(bool)}, got {type(force_download)})")
        self.force_download = force_download
        self.gene_aliases_mapping = self.__aliases_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_aliases_mapping["genename"].items()}
        self.databases = {_db for _db in self.gene_aliases_mapping["databases"].keys()}
        return None
    
    def __call__(self, data: Union[Sequence[Tuple[str, str, Dict[str, int]]], DataFrame, Graph], *args, **kwargs):
        if (isinstance(data, SequenceInstance) and not isinstance(data, str)) or isinstance(data, set):
            return self.sequence_standardization(data, *args, **kwargs)
        elif isinstance(data, DataFrame):
            return self.df_standardization(data, *args, **kwargs)
        elif isinstance(data, Graph):
            return self.graph_standardization(data, *args, **kwargs)
        else:
            raise TypeError(f"fail to convert gene name: 'data' has incorrect type")

    def __aliases_from_NCBI(self, gi_file: Path) -> dict:
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

        gene_aliases_mapping = {
            "geneid": dict(),
            "genename": dict(),
            "ensemblid": dict(),
            "databases": dict()
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
                _ensemblid = None
                _db_name_list = geneinfo.pop(0).split("|")
                _db_name_dict = dict()
                for _db_name in _db_name_list:
                    if _db_name == "-":
                        continue
                    _db, _name = _db_name.split(":", maxsplit=1)
                    if _db == "Ensembl":
                        _ensemblid = re.findall("[A-Z]{7}[0-9]{11}", _name)
                        _ensemblid = _ensemblid[0] if _ensemblid else None
                    else:
                        _db_name_dict[_db] = _name
                
                gene_names.add(_reference_name)
                pointer_to_geneid = ctypes.create_string_buffer(_geneid.encode())
                if _geneid:
                    gene_aliases_mapping["geneid"][_geneid] = namedtuple("GeneAlias", ["reference_genename", "ensemblid", "databases"])(_reference_name, _ensemblid, _db_name_dict)
                gene_aliases_mapping["genename"][_reference_name] = pointer_to_geneid
                for _synonym in _synonyms:
                    if _synonym not in gene_names:
                        gene_names.add(_synonym)
                        gene_aliases_mapping["genename"][_synonym] = pointer_to_geneid
                    else:
                        print(f"synonym {_synonym} multiple times with ref name: {_reference_name}")
                if _ensemblid:
                    gene_aliases_mapping["ensemblid"][_ensemblid] = pointer_to_geneid
                if _db_name_dict:
                    for _db, _name in _db_name_dict.items():
                        if _db not in gene_aliases_mapping["databases"]:
                            gene_aliases_mapping["databases"][_db] = dict()
                        gene_aliases_mapping["databases"][_db][_name] = pointer_to_geneid

        os.system(f"rm {str(gi_file_cut)}")
        return gene_aliases_mapping
    
    def get_databases(self):
        """
        Provide database names using specific nomenclature.
        """
        return self.databases

    def get_geneid(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the geneid with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|ensemblid|<database>
            see self.get_database() for enumerating database names

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
        elif alias_type == "ensemblid":
            if alias in self.gene_aliases_mapping[alias_type]:
                return self.gene_aliases_mapping[alias_type][alias].value.decode()
            else:
                warnings.warn(f"no geneid correspondance for {alias_type} '{alias}'", stacklevel=10)
                return None
        elif alias_type in self.databases:
            if alias in self.gene_aliases_mapping["databases"][alias_type]:
                return self.gene_aliases_mapping["databases"][alias_type][alias].value.decode()
            else:
                warnings.warn(f"no geneid correspondance for {alias_type} '{alias}'", stacklevel=10)
                return None
        else:
            raise ValueError("invalid value for argument 'alias_type'")
    
    def get_reference_name(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the reference name with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names

        Returns
        -------
        Given an alias, return its reference name
        """
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_aliases_mapping["geneid"]:
            return self.gene_aliases_mapping["geneid"][geneid].reference_genename
        else:
            warnings.warn(f"no reference genename correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None

    def get_ensemblid(self, alias: str, alias_type: str="genename") -> str:
        """
        Provide the ensemblid with respect to a gene alias.

        Parameters
        ----------
        alias
            gene alias
        alias_type
            genename|geneid|<database>
            see self.get_database() for enumerating database names

        Returns
        -------
        Given an alias, return its ensemblid
        """
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_aliases_mapping["geneid"]:
            return self.gene_aliases_mapping["geneid"][geneid].ensemblid
        else:
            warnings.warn(f"no ensemblid correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None
    
    def get_alias_from_database(self, database: str, alias: str, alias_type: str="genename") -> str:
        """
        Provide the database-defined gene name.

        Parameters
        ----------
        database
            see self.get_database() for enumerating database names
        alias
            gene alias
        alias_type
            genename|geneid|ensemblid

        Returns
        -------
        Given a gene identifier derived from a database, return its alias
        """

        if database not in self.get_databases():
            raise ValueError(f"invalid value for argument 'database' (got {database}, expected a value in {self.get_databases})")
    
        geneid = self.get_geneid(alias, alias_type) if alias_type != "geneid" else alias
        if geneid in self.gene_aliases_mapping["geneid"]:
            if database in self.gene_aliases_mapping["geneid"][geneid].databases:
                return self.gene_aliases_mapping["geneid"][geneid].databases[database]
            else:
                warnings.warn(f"no {database} correspondance for {alias_type} '{alias}'", stacklevel=10)
                return None
        else:
            warnings.warn(f"no {database} correspondance for {alias_type} '{alias}'", stacklevel=10)
            return None

    def __convert(self, out_alias_type: str="referencename") -> str:
        """
        Convert gene alias.

        Parameters
        ----------
        out_alias_type
            geneid|referencename|ensemblid|<database>
            see self.get_database() for enumerating database names

        Returns
        -------
        Return a function converting gene labels.
        """

        if out_alias_type == "referencename":
            out_alias_type = "reference_name"
        if out_alias_type in ["geneid", "reference_name", "ensemblid"]:
            return eval(f"self.get_{out_alias_type}")
        elif out_alias_type in self.get_databases():
            return partial(self.get_alias_from_database, database=out_alias_type)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute 'get_{out_alias_type}'")

    def sequence_standardization(
        self,
        gene_sequence: Sequence[str],
        in_alias_type: str="genename",
        out_alias_type: str="referencename",
        keep_if_missing: bool=True
    ) -> Sequence[str]:
        """
        Create a copy of the input Sequence, with corresponding aliases.

        Parameters
        ----------
        gene_sequence
            list of genes
        in_alias_type
            genename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        out_alias_type
            referencename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        keep_if_missing
            if conversion is performed from gene names towards theirs reference gene names,
            keep gene name instead of 'None' value if gene name is missing from NCBI database.
        
        Returns
        -------
        return a gene sequence where each gene alias is converted into the user-defined alias type.
        """
        
        keep_if_missing = keep_if_missing if (in_alias_type=="genename" and out_alias_type=="referencename") else False
        standardized_gene_sequence = list()
        alias_conversion = self.__convert(out_alias_type)
        for gene in gene_sequence:
            alias = alias_conversion(alias=gene, alias_type=in_alias_type)
            alias = gene if (keep_if_missing and alias is None) else alias
            standardized_gene_sequence.append(alias)
        
        standardized_gene_sequence = type(gene_sequence)(standardized_gene_sequence)
        
        return standardized_gene_sequence

    def interaction_list_standardization(
        self,
        interactions_list: Sequence[Tuple[str, str, Dict[str, int]]],
        in_alias_type: str="genename",
        out_alias_type: str="referencename",
        keep_if_missing: bool=True
    ) -> List[Tuple[str, str, Dict[str, int]]]:
        """
        Create a copy of the input list of pairwise interactions, with corresponding aliases.

        Parameters
        ----------
        interaction_list
            list of tuples containing string (source) + string (target) + dict (sign = -1 or 1)
        in_alias_type
            genename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        out_alias_type
            referencename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        keep_if_missing
            if conversion is performed from gene names towards theirs reference gene names,
            keep gene name instead of 'None' value if gene name is missing from NCBI database.

        Returns
        -------
        return an interaction list where each gene name is converted into the user-defined alias type.
        """

        keep_if_missing = keep_if_missing if (in_alias_type=="genename" and out_alias_type=="referencename") else False
        standardized_interactions_list = list()
        alias_conversion = self.__convert(out_alias_type)
        for interaction in interactions_list:
            source = alias_conversion(alias=interaction[0], alias_type=in_alias_type)
            source = interaction[0] if (keep_if_missing and source is None) else source
            target = alias_conversion(alias=interaction[0], alias_type=in_alias_type)
            target = interaction[1] if (keep_if_missing and target is None) else target
            standardized_interactions_list.append((source, target, interaction[2]))

        return standardized_interactions_list

    def df_standardization(
        self,
        df: DataFrame,
        axis: Axis=0,
        in_alias_type: str="genename",
        out_alias_type: str="referencename",
        keep_if_missing: bool=True,
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Replace gene name with its reference gene name into 'df'.

        Parameters
        ----------
        df
            DataFrame object where names are expected being standardized
        axis
            whether to rename labels from the index (0 or 'index') or columns (1 or 'columns')
        in_alias_type
            genename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        out_alias_type
            referencename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        keep_if_missing
            if conversion is performed from gene names towards theirs reference gene names,
            keep gene name instead of 'None' value if gene name is missing from NCBI database.
        copy
            return a copy instead of updating 'df'
        
        Returns
        -------
        Depending on 'copy', update or return DataFrame object with standardized gene name.
        """

        keep_if_missing = keep_if_missing if (in_alias_type=="genename" and out_alias_type=="referencename") else False
        df = df.copy() if copy is True else df
        alias_conversion = self.__convert(out_alias_type)

        aliases = list()

        if axis == 0 or axis == "index":
            gene_iterator = iter(df.index)
        elif axis == 1 or axis == "columns":
            gene_iterator = iter(df.columns)
        else:
            raise ValueError(f"No axis named {axis} for object type DataFrame")

        for gene in gene_iterator:
            alias = alias_conversion(alias=gene, alias_type=in_alias_type)
            alias = gene if (keep_if_missing and alias is None) else alias
            aliases.append(alias)

        if axis == 0 or axis == "index":
            df.index = aliases
        elif axis == 1 or axis == "columns":
            df.columns = aliases

        if copy is True:
            return df
    
    def graph_standardization(
        self,
        graph: Graph,
        in_alias_type: str="genename",
        out_alias_type: str="referencename",
        keep_if_missing: bool=True,
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene name with its reference gene name into 'graph'.

        Parameters
        ----------
        graph
            graph where nodes must be standardized
        in_alias_type
            genename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        out_alias_type
            referencename|geneid|ensemblid|<database>
            see self.get_database() for enumerating database names
        keep_if_missing
            if conversion is performed from gene names towards theirs reference gene names,
            keep gene name instead of 'None' value if gene name is missing from NCBI database.
        copy
            return a copy instead of updating 'graph'
        
        Returns
        -------
        Depending on 'copy', update or return graph with standardized gene name.
        """

        keep_if_missing = keep_if_missing if (in_alias_type=="genename" and out_alias_type=="referencename") else False
        aliases_mapping = dict()
        alias_conversion = self.__convert(out_alias_type)
        for gene in graph.nodes:
            alias = alias_conversion(alias=gene, alias_type=in_alias_type)
            alias = gene if (keep_if_missing and alias is None) else alias
            aliases_mapping[gene] = alias
        if copy is True:
            return nx.relabel_nodes(graph, mapping=aliases_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=aliases_mapping, copy=False)
            return None
