#!/usr/bin/env python

from typing import (
    Union,
    Mapping,
    Sequence,
    Dict,
    Tuple,
    Literal,
    Callable,
    Any,
    get_args
)
from pandas import DataFrame
from pandas._typing import Axis
try:
    from collections import Sequence as SequenceInstance
except:
    from collections.abc import Sequence as SequenceInstance
from networkx import Graph
from .._typing import MPBooleanNetwork

import ctypes
from collections import namedtuple
from functools import partial
try:
    from sortedcontainers import SortedSet
except:
    pass

import os, warnings
from pathlib import Path

import re

import networkx as nx

GeneType = Literal["genename", "geneid", "ensemblid"]
AliasType = Literal["referencename", "geneid", "ensemblid"]
InteractionList = Sequence[Tuple[str, str, Dict[str, int]]]

ORGANISMS = Literal[
    "mouse",
    "human",
    "escherichia coli"
]

FTP_GENE_INFO = {
    "mouse": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    "human": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "escherichia coli": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz"
}

NCBI_FILES = {
    "mouse": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/data/gi/mus_musculus_gene_info.tsv"),
    "human": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/data/gi/homo_sapiens_gene_info.tsv"),
    "escherichia coli": Path(f"{str(os.path.abspath(os.path.dirname(__file__)))}/data/gi/escherichia_coli_gene_info.tsv")
}

class GeneSynonyms(object):
    """
    Mapping between gene aliases.

    Parameters
    ----------
    organism: str (default: mouse)
        Common name of the organism of interest.
    force_download: bool (default: False)
        Request to the ncbi ftp protocol for downloading gene info data.
    show_warnings: bool (default: False)
        Print warning messages.
    """

    def __init__(
        self,
        organism: str = "mouse",
        force_download: bool = False,
        show_warnings: bool = False
    ) -> None:
        
        def valid_alias_type(*args):
            return Literal[args]

        organism = organism.lower().replace("-"," ")

        if organism in get_args(ORGANISMS):
            self.organism = organism
        else:
            raise ValueError(f"invalid argument value for 'organism': {organism}")
        self.ncbi_file = NCBI_FILES[organism]

        if force_download is True:
            command_parsing = f"wget --quiet --show-progress -cO {self.ncbi_file}.gz {FTP_GENE_INFO[self.organism]} && gunzip --quiet {self.ncbi_file}.gz"
            os.system(command_parsing)
        elif force_download is False:
            pass
        else:
            raise TypeError(f"unsupported argument type for '{force_download}': expected {bool} but received {type(force_download)}")

        self.force_download = force_download
        self.show_warnings = show_warnings
        self.gene_aliases_mapping = self.__aliases_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_aliases_mapping["genename"].items()}
        self.databases = Literal[{db for db in self.gene_aliases_mapping["databases"].keys()}]
        self.valid_alias_types = valid_alias_type("genename", "geneid", "ensemblid", *self.gene_aliases_mapping["databases"].keys())

        return None
    
    def __get__(
        self,
        attribute: str = "gene_aliases_mapping"
    ) -> Any:
    
        if attribute == "__upper_gene_names_mapping":
            raise AttributeError(f"attribute '{attribute}' not accessed: private")
        else:
            return getattr(self, attribute)
    
    def __set__(
        self,
        organism: str = None,
        force_download: bool = False,
        show_warnings: bool = False
    ) -> None:
        
        def valid_alias_type(*args):
            return Literal[args]
    
        if organism is None and self.organism in get_args(ORGANISMS):
            pass
        elif organism in get_args(ORGANISMS):
            self.organism = organism
        else:
            raise ValueError(f"invalid argument value for 'organism': {organism}")
        self.ncbi_file = NCBI_FILES[self.organism]
        
        if force_download is True:
            command_parsing = f"wget --quiet --show-progress -cO {self.ncbi_file}.gz {FTP_GENE_INFO[self.organism]} && gunzip --quiet {self.ncbi_file}.gz"
            os.system(command_parsing)
        elif force_download is False:
            pass
        else:
            raise TypeError(f"unsupported argument type for '{force_download}': expected {bool} but received {type(force_download)}")
        
        self.force_download = force_download
        self.show_warnings = show_warnings
        self.gene_aliases_mapping = self.__aliases_from_NCBI(self.ncbi_file)
        self.__upper_gene_names_mapping = {key.upper(): value for key, value in self.gene_aliases_mapping["genename"].items()}
        self.databases = {_db for _db in self.gene_aliases_mapping["databases"].keys()}
        self.valid_alias_types = valid_alias_type("genename", "geneid", "ensemblid", *self.gene_aliases_mapping["databases"].keys())

        return None
    
    def __call__(
        self,
        data: Union[InteractionList, DataFrame, Graph, MPBooleanNetwork], # type: ignore
        *args: Sequence[Any],
        **kwargs: Mapping[str, Any]
    ):
        if (isinstance(data, SequenceInstance) and not isinstance(data, str)) or isinstance(data, set):
            return self.convert_sequence(data, *args, **kwargs)
        elif isinstance(data, DataFrame):
            return self.convert_df(data, *args, **kwargs)
        elif isinstance(data, Graph):
            return self.convert_graph(data, *args, **kwargs)
        elif isinstance(data, MPBooleanNetwork):
            return self.convert_bn(data, *args, **kwargs)
        else:
            raise TypeError(f"unsupported argument type for 'data': {data}")

    def __aliases_from_NCBI(self, gi_file: Path) -> Dict:
        """
        Create a dictionary matching each gene identifier to its NCBI reference gene name.
        For speeding up the task facing a large matrix from NCBI, the parsing of the NCBI gene data is run with awk.

        Parameters
        ----------
        gi_file: Path
            Path to the NCBI gene info data.

        Returns
        -------
        Return a mapping between gene identifiers and NCBI reference gene names.
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
                    gene_aliases_mapping["geneid"][_geneid] = namedtuple("Identifiers", ["reference_genename", "ensemblid", "databases"])(_reference_name, _ensemblid, _db_name_dict)
                gene_aliases_mapping["genename"][_reference_name] = pointer_to_geneid
                for _synonym in _synonyms:
                    if _synonym not in gene_names:
                        gene_names.add(_synonym)
                        gene_aliases_mapping["genename"][_synonym] = pointer_to_geneid
                    elif self.show_warnings:
                        warnings.warn(f"synonym {_synonym} multiple times with reference name: {_reference_name}", stacklevel=10)
                if _ensemblid:
                    gene_aliases_mapping["ensemblid"][_ensemblid] = pointer_to_geneid
                if _db_name_dict:
                    for _db, _name in _db_name_dict.items():
                        if _db not in gene_aliases_mapping["databases"]:
                            gene_aliases_mapping["databases"][_db] = dict()
                        gene_aliases_mapping["databases"][_db][_name] = pointer_to_geneid

        os.remove(gi_file_cut)

        return gene_aliases_mapping
    
    def get_geneid(
        self,
        gene: str,
        gene_type: Union[Literal["genename", "ensemblid"], str] = "genename"
    ) -> str:
        """
        Provide the geneid with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        gene_type: 'genename' | 'ensemblid' | <database> (default: 'genename')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its geneid.
        """

        if gene_type == "genename":
            gene_name = gene.upper()
            if gene_name in self.__upper_gene_names_mapping:
                return self.__upper_gene_names_mapping[gene_name].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no correspondance for gene '{gene}'", stacklevel=10)
                return None
        elif gene_type == "ensemblid":
            if gene in self.gene_aliases_mapping[gene_type]:
                return self.gene_aliases_mapping[gene_type][gene].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no geneid correspondance for {gene_type} '{gene}'", stacklevel=10)
                return None
        elif gene_type in self.databases:
            if gene in self.gene_aliases_mapping["databases"][gene_type]:
                return self.gene_aliases_mapping["databases"][gene_type][gene].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no geneid correspondance for {gene_type} '{gene}'", stacklevel=10)
                return None
        else:
            raise ValueError(f"invalid argument value for 'gene_type': {gene_type}")
    
    def get_reference_name(
        self,
        gene: str,
        gene_type: Union[GeneType, str] = "genename"
    ) -> str:
        """
        Provide the NCBI reference name with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its NCBI reference name.
        """
    
        geneid = self.get_geneid(gene, gene_type) if gene_type != "geneid" else gene
        if geneid in self.gene_aliases_mapping["geneid"]:
            return self.gene_aliases_mapping["geneid"][geneid].reference_genename
        else:
            if self.show_warnings:
                warnings.warn(f"no reference genename correspondance for {gene_type} '{gene}'", stacklevel=10)
            return None

    def get_ensemblid(
        self,
        gene: str,
        gene_type: Union[Literal["genename", "geneid"], str] = "genename"
    ) -> str:
        """
        Provide the ensemblid with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        gene_type: 'genename' | 'geneid' | <database> (default: 'genename')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its ensemblid.
        """
    
        geneid = self.get_geneid(gene, gene_type) if gene_type != "geneid" else gene
        if geneid in self.gene_aliases_mapping["geneid"]:
            return self.gene_aliases_mapping["geneid"][geneid].ensemblid
        else:
            if self.show_warnings:
                warnings.warn(f"no ensemblid correspondance for {gene_type} '{gene}'", stacklevel=10)
            return None
    
    def get_alias_from_database(
        self,
        gene: str,
        database: str,
        gene_type: GeneType = "genename"
    ) -> str:
        """
        Provide the database-defined gene name with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        database: <database>
            Organism-related database name providing gene identifiers.
            See self.databases for enumerating valid database names.
        gene_type: 'genename' | 'geneid' | 'ensemblid' (default: 'genename')
            Identifier type of the given gene.

        Returns
        -------
        Given a gene identifier, return its alias derived from a database.
        """

        if database not in self.databases:
            raise ValueError(f"invalid argument value for 'database': got '{database}' but expected a value in {self.databases})")
    
        geneid = self.get_geneid(gene, gene_type) if gene_type != "geneid" else gene
        if geneid in self.gene_aliases_mapping["geneid"]:
            if database in self.gene_aliases_mapping["geneid"][geneid].databases:
                return self.gene_aliases_mapping["geneid"][geneid].databases[database]
            else:
                if self.show_warnings:
                    warnings.warn(f"no {database} correspondance for {gene_type} '{gene}'", stacklevel=10)
                return None
        else:
            if self.show_warnings:
                warnings.warn(f"no {database} correspondance for {gene_type} '{gene}'", stacklevel=10)
            return None

    def conversion(
        self,
        gene: str,
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename"
    ) -> str:
        """
        Convert gene identifiers into the user-defined alias type.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Input identifier type of the given gene.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Output identifier type for the given gene.

        Returns
        -------
        Given a gene identifier, return its alias.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        if alias_type == "referencename":
            alias_type = "reference_name"
        if alias_type in ["geneid", "reference_name", "ensemblid"]:
            convert = eval(f"self.get_{alias_type}")
        elif alias_type in self.databases:
            convert = partial(self.get_alias_from_database, database=alias_type)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute 'get_{alias_type}'")
        
        return convert(gene=gene, gene_type=gene_type)
    
    def __conversion_function(
        self,
        alias_type: Union[AliasType, str] = "referencename",
        *args: Sequence[Any],
        **kwargs: Mapping[str, Any]
    ) -> Callable:
        """
        Function converting gene identifiers.

        Parameters
        ----------
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.

        Returns
        -------
        Return Function object converting gene identifiers.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        if alias_type == "referencename":
            alias_type = "reference_name"
        if alias_type in ["geneid", "reference_name", "ensemblid"]:
            return eval(f"self.get_{alias_type}")
        elif alias_type in self.databases:
            return partial(self.get_alias_from_database, database=alias_type)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute 'get_{alias_type}'")

    def convert_sequence(
        self,
        genes: Sequence[str],
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename",
        keep_if_missing: bool = True
    ) -> Sequence[str]:
        """
        Create a copy of the Sequence object, with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        genes: Sequence[str]
            List of gene identifiers.
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Gene identifier input format.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.
        keep_if_missing: bool (default: True)
            If true, keep origin gene identifier instead of None value if origin gene identifier is missing from NCBI database.
        
        Returns
        -------
        Return Sequence object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """
        
        # keep_if_missing = keep_if_missing if (gene_type=="genename" and alias_type=="referencename") else False
        aliases = list()
        alias_conversion = self.__conversion_function(alias_type)
        for gene in genes:
            output_alias = alias_conversion(gene=gene, gene_type=gene_type)
            output_alias = gene if (keep_if_missing and output_alias is None) else output_alias
            aliases.append(output_alias)
        
        aliases = type(genes)(aliases)
        
        return aliases

    def convert_interaction_list(
        self,
        interaction_list: InteractionList,
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename",
        keep_if_missing: bool = True
    ) -> InteractionList:
        """
        Create a copy of the pairwise InteractionList object, with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            List of tuples containing string (source) + string (target) + dict (sign = -1 or 1).
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Gene identifier input format.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.
        keep_if_missing: bool (default: True)
            If true, keep origin gene alias instead of None value if origin gene alias is missing from NCBI database.
        
        Returns
        -------
        Return InteractionList object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        # keep_if_missing = keep_if_missing if (gene_type=="genename" and alias_type=="referencename") else False
        converted_interactions_list = list()
        alias_conversion = self.__conversion_function(alias_type)
        for interaction in interaction_list:
            source = alias_conversion(gene=interaction[0], gene_type=gene_type)
            source = interaction[0] if (keep_if_missing and source is None) else source
            target = alias_conversion(gene=interaction[0], gene_type=gene_type)
            target = interaction[1] if (keep_if_missing and target is None) else target
            converted_interactions_list.append((source, target, interaction[2]))

        return converted_interactions_list

    def convert_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename",
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Replace gene identifiers in DataFrame object with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        df: pd.DataFrame
            DataFrame object where gene aliases are converted into the desired identifiers.
        axis
            whether to rename labels from the index (0 or 'index') or columns (1 or 'columns').
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Gene identifier input format.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.
        copy: bool (default: True)
            Return a copy instead of updating DataFrame object.
        
        Returns
        -------
        Depending on 'copy', update 'df' or return DataFrame object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        df = df.copy() if copy is True else df
        alias_conversion = self.__conversion_function(alias_type)

        genes = list()

        if axis == 0 or axis == "index":
            iterator = iter(df.index)
        elif axis == 1 or axis == "columns":
            iterator = iter(df.columns)
        else:
            raise TypeError(f"unsupported argument type for 'axis': {axis}")

        for gene in iterator:
            output_alias = alias_conversion(gene=gene, gene_type=gene_type)
            output_alias = gene if output_alias is None else output_alias
            genes.append(output_alias)

        if axis == 0 or axis == "index":
            df.index = genes
        elif axis == 1 or axis == "columns":
            df.columns = genes

        if copy is True:
            return df
    
    def convert_graph(
        self,
        graph: Graph,
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename",
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene identifiers in Graph object with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        graph: nx.Graph
            Graph object where gene nodes are converted into the desired identifiers.
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Gene identifier input format.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.
        copy: bool (default: True)
            Return a copy instead of updating Graph object.

        Returns
        -------
        Depending on 'copy', update 'df' or return Graph object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        aliases_mapping = dict()
        alias_conversion = self.__conversion_function(alias_type)
        for gene in graph.nodes:
            output_alias = alias_conversion(gene=gene, gene_type=gene_type)
            output_alias = gene if output_alias is None else output_alias
            aliases_mapping[gene] = output_alias
        if copy is True:
            return nx.relabel_nodes(graph, mapping=aliases_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=aliases_mapping, copy=False)
            return None

    def convert_bn(
        self,
        bn: MPBooleanNetwork, # type: ignore
        gene_type: Union[GeneType, str] = "genename",
        alias_type: Union[AliasType, str] = "referencename",
        copy: bool = False
    ) -> MPBooleanNetwork: # type: ignore
        """
        Replace gene identifiers in MPBooleanNetwork object with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        bn: MPBooleanNetwork
            MPBooleanNetwork object where gene variables are converted into the desired identifiers.
        gene_type: 'genename' | 'geneid' | 'ensemblid' | <database> (default: 'genename')
            Gene identifier input format.
        alias_type: 'referencename' | 'geneid' | 'ensemblid' | <database> (default: 'referencename')
            Gene identifier output format.
        copy: bool (default: True)
            Return a copy instead of updating MPBooleanNetwork object.
        
        Returns
        -------
        Depending on 'copy', update 'bn' or return MPBooleanNetwork object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        bn = bn.copy() if copy else bn

        alias_conversion = self.__conversion_function(alias_type)
        genes = tuple(bn.keys())
        for gene in genes:
            output_alias = alias_conversion(gene=gene, gene_type=gene_type)
            output_alias = gene if output_alias is None else output_alias
            bn.rename(gene, output_alias)
                
        return bn if copy else None

    def standardize_sequence(
        self,
        genes: Sequence[str],
        keep_if_missing: bool = True
    ) -> Sequence[str]:
        """
        Create a copy of the Sequence object, with corresponding NCBI reference names.
        Each gene name is converted into its NCBI reference name.

        Parameters
        ----------
        genes: Sequence[str]
            List of gene names.
        keep_if_missing: bool (default: True)
            If true, keep origin gene name instead of None value if origin gene name is missing from NCBI database.
        
        Returns
        -------
        Return Sequence object.
        """

        return self.convert_sequence(
            genes=genes,
            gene_type="genename",
            alias_type="referencename",
            keep_if_missing=keep_if_missing
        )

    def standardize_interaction_list(
        self,
        interaction_list: InteractionList,
        keep_if_missing: bool = True
    ) -> InteractionList:
        """
        Create a copy of the pairwise InteractionList object, with corresponding NCBI reference names.
        Each gene name is converted into its NCBI reference name.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            List of tuples containing string (source) + string (target) + dict (sign = -1 or 1).
        keep_if_missing: bool (default: True)
            If true, keep origin gene alias instead of None value if origin gene alias is missing from NCBI database.
        
        Returns
        -------
        Return InteractionList object.
        """

        return self.convert_interaction_list(
            interaction_list=interaction_list,
            gene_type="genename",
            alias_type="referencename",
            keep_if_missing=keep_if_missing
        )
    
    def standardize_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Replace gene names in DataFrame object with corresponding NCBI reference names.

        Parameters
        ----------
        df: pd.DataFrame
            DataFrame object where gene identifiers are expected being standardized.
        axis
            whether to rename labels from the index (0 or 'index') or columns (1 or 'columns').
        copy: bool (default: True)
            Return a copy instead of updating DataFrame object.
        
        Returns
        -------
        Depending on 'copy', update 'df' or return DataFrame object.
        """

        return self.convert_df(
            df=df,
            gene_type="genename",
            alias_type="referencename",
            axis=axis,
            copy=copy
        )

    def standardize_graph(
        self,
        graph: Graph,
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene names in Graph object with corresponding NCBI reference names.

        Parameters
        ----------
        graph: nx.Graph
            Graph object where gene nodes are expected being standardized.
        copy: bool (default: True)
            Return a copy instead of updating Graph object.

        Returns
        -------
        Depending on 'copy', update 'graph' or return Graph object.
        """

        return self.convert_graph(
            graph=graph,
            gene_type="genename",
            alias_type="referencename",
            copy=copy
        )

    def standardize_bn(
        self,
        bn: MPBooleanNetwork, # type: ignore
        copy: bool = False
    ) -> MPBooleanNetwork: # type: ignore
        """
        Replace gene names in MPBooleanNetwork object with corresponding NCBI reference names.

        Parameters
        ----------
        bn: MPBooleanNetwork
            MPBooleanNetwork object where gene variables are expected being standardized.
        copy: bool (default: True)
            Return a copy instead of updating MPBooleanNetwork object.
        
        Returns
        -------
        Depending on 'copy', update 'bn' or return MPBooleanNetwork object.
        """

        return self.convert_bn(
            bn=bn,
            gene_type="genename",
            alias_type="referencename",
            copy=copy
        )
