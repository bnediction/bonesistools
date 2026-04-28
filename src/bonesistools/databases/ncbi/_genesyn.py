#!/usr/bin/env python

from typing import (
    Union,
    Optional,
    Mapping,
    Sequence,
    Dict,
    Tuple,
    Callable,
    Any,
)
try:
    from typing import Literal, get_args
except ImportError:
    from typing_extensions import Literal, get_args # type: ignore
try:
    from collections import Sequence as SequenceInstance
except:
    from collections.abc import Sequence as SequenceInstance

import copy
import ctypes
import inspect
from collections import namedtuple
from functools import partial, wraps

from pandas import DataFrame
from pandas._typing import Axis

from networkx import Graph

import subprocess
import warnings
from pathlib import Path

import re

import networkx as nx
from .._typing import MPBooleanNetwork

InputIdentifierType = Literal[
    "name",
    "gene_id",
    "ensembl_id"
]

OutputIdentifierType = Literal[
    "official_name",
    "ncbi_name",
    "gene_id",
    "ensembl_id"
]

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

NCBI_DIR = Path(__file__).resolve().parent / "data" / "gi"

NCBI_GENE_INFO_FILES = {
    "mouse": NCBI_DIR / "mus_musculus_gene_info.tsv",
    "human": NCBI_DIR / "homo_sapiens_gene_info.tsv",
    "escherichia coli": NCBI_DIR / "escherichia_coli_gene_info.tsv",
}

GENE_TYPE_PRIORITY = {
    "protein-coding": 100,
    "ncRNA": 80,
    "rRNA": 70,
    "scRNA": 70,
    "snRNA": 70,
    "snoRNA": 70,
    "tRNA": 70,
    "other": 30,
    "unknown": 20,
    "biological-region": 10,
    "pseudo": 0
}

_GENE_SYNONYMS_DEPRECATED_ARGS = {
    "gene_type": "input_identifier_type",
    "alias_gene": "output_identifier_type",
}

def support_legacy_gene_synonyms_args(func):
    
    valid_parameters = inspect.signature(func).parameters
    
    @wraps(func)
    def wrapper(*args, **kwargs):
        for old_name, new_name in _GENE_SYNONYMS_DEPRECATED_ARGS.items():

            if old_name not in kwargs:
                continue

            if new_name not in valid_parameters:
                continue

            warnings.warn(
                f"'{old_name}' is deprecated and will be removed in a future version; "
                f"use '{new_name}' instead.",
                FutureWarning,
                stacklevel=2,
            )

            if new_name in kwargs:
                raise TypeError(
                    f"Use either '{old_name}' or '{new_name}', not both."
                )

            kwargs[new_name] = kwargs.pop(old_name)

        return func(*args, **kwargs)

    return wrapper

class GeneSynonyms:
    """
    Mapping between gene aliases.

    Parameters
    ----------
    organism: str (default: mouse)
        Common name of the organism of interest.
    force_download: bool (default: False)
        Request to the ncbi ftp protocol for downloading gene_info data.
    show_warnings: bool (default: False)
        Print warning messages.
    """

    def __init__(
        self,
        organism: str = "mouse",
        force_download: bool = False,
        show_warnings: bool = False,
    ) -> None:

        organism = organism.lower().replace("-", " ")

        if organism not in get_args(ORGANISMS):
            raise ValueError(
                f"invalid argument value for 'organism': {organism}"
            )

        if not isinstance(force_download, bool):
            raise TypeError(
                f"unsupported argument type for 'force_download': "
                f"expected {bool} but received {type(force_download)}"
            )

        if not isinstance(show_warnings, bool):
            raise TypeError(
                f"unsupported argument type for 'show_warnings': "
                f"expected {bool} but received {type(show_warnings)}"
            )

        self.organism = organism
        self.ncbi_file = NCBI_GENE_INFO_FILES[self.organism]
        self.force_download = force_download

        self.__download_gene_info(force_download=force_download)
        self.__initialize_mappings(show_warnings=show_warnings)

    def reset(
        self,
        organism: str = None,
        force_download: bool = False,
        show_warnings: bool = False,
    ) -> None:

        if organism is None:
            organism = self.organism
        else:
            organism = organism.lower().replace("-", " ")

        if organism not in get_args(ORGANISMS):
            raise ValueError(
                f"invalid argument value for 'organism': {organism}"
            )

        if not isinstance(force_download, bool):
            raise TypeError(
                f"unsupported argument type for 'force_download': "
                f"expected {bool} but received {type(force_download)}"
            )

        if not isinstance(show_warnings, bool):
            raise TypeError(
                f"unsupported argument type for 'show_warnings': "
                f"expected {bool} but received {type(show_warnings)}"
            )

        self.organism = organism
        self.ncbi_file = NCBI_GENE_INFO_FILES[self.organism]
        self.force_download = force_download

        self.__download_gene_info(force_download=force_download)
        self.__initialize_mappings(show_warnings=show_warnings)

    def __download_gene_info(self, force_download: bool = False) -> None:

        if not force_download:
            return

        cmd = (
            f"wget --quiet --show-progress -cO {self.ncbi_file}.gz "
            f"{FTP_GENE_INFO[self.organism]} && "
            f"gunzip --quiet {self.ncbi_file}.gz"
        )

        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("failed to download NCBI gene_info file") from e

    def __parse_ncbi_gene_info(self, gi_file: Path) -> Dict:
        """
        Parse NCBI gene_info file and create mappings between multiple gene identifiers.

        Parameters
        ----------
        gi_file: Path
            Path to the NCBI gene_info data.

        Returns
        -------
        Dict
            Mapping structure for gene identifiers and aliases.
        """

        gi_file_cut = Path(f"{gi_file}_cut")

        # Output columns:
        # gene_id | official_name | ncbi_name | synonyms | dbXrefs | gene_type
        cmd = (
            "awk -F'\t' 'NR>1 {print $2 \"\t\" $11 \"\t\" $3 \"\t\" $5 \"\t\" $6 \"\t\" $10}' "
            + str(gi_file) + " > " + str(gi_file_cut)
        )
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("failed to parse NCBI gene_info file") from e

        Identifiers = namedtuple(
            "Identifiers",
            ["official_name", "ncbi_name", "ensembl_id", "databases", "gene_type"]
        )

        gene_aliases_mapping = {
            "gene_id": dict(),
            "name": dict(),
            "ensembl_id": dict(),
            "databases": dict(),
        }

        gene_info_rows = []
        protected_official_names = set()

        try:
            with open(gi_file_cut, "r") as file:
                for line in file:
                    fields = line.rstrip("\n").split("\t")

                    gene_id = fields[0]
                    official_name = fields[1]
                    ncbi_name = fields[2]
                    gene_type = fields[5]
                    if official_name == "-":
                        official_name = ncbi_name

                    gene_info_rows.append((
                        gene_id,
                        official_name,
                        ncbi_name,
                        fields[3],
                        fields[4],
                        gene_type,
                    ))

                    if GENE_TYPE_PRIORITY.get(gene_type, -1) > 0:
                        protected_official_names.add(official_name.upper())

            for (
                gene_id,
                official_name,
                ncbi_name,
                synonyms_field,
                db_xrefs_field,
                gene_type,
            ) in gene_info_rows:

                synonyms = [
                    synonym
                    for synonym in (
                        synonyms_field.split("|") + [ncbi_name]
                    )
                    if synonym != "-" and synonym.upper() != official_name.upper()
                ]

                ensembl_id = None
                database_aliases = {}

                for db_entry in db_xrefs_field.split("|"):
                    if db_entry == "-":
                        continue

                    database, database_name = db_entry.split(":", maxsplit=1)

                    if database == "Ensembl":
                        ensembl_match = re.findall(
                            "[A-Z]{7}[0-9]{11}",
                            database_name
                        )
                        ensembl_id = ensembl_match[0] if ensembl_match else None
                    else:
                        database_aliases[database] = database_name

                gene_id_pointer = ctypes.create_string_buffer(gene_id.encode())

                gene_aliases_mapping["gene_id"][gene_id] = Identifiers(
                    official_name=official_name,
                    ncbi_name=ncbi_name,
                    ensembl_id=ensembl_id,
                    databases=database_aliases,
                    gene_type=gene_type,
                )

                priority = GENE_TYPE_PRIORITY.get(gene_type, -1)

                # Official name always wins
                gene_aliases_mapping["name"][official_name.upper()] = (
                    gene_id_pointer,
                    priority,
                )

                for synonym in synonyms:
                    synonym_upper = synonym.upper()
                    # Prevent overwriting another protected official name
                    if synonym_upper in protected_official_names:
                        continue
                    if synonym_upper not in gene_aliases_mapping["name"]:
                        gene_aliases_mapping["name"][synonym_upper] = (gene_id_pointer, priority)
                    else:
                        _, old_priority = gene_aliases_mapping["name"][synonym_upper]
                        if priority > old_priority:
                            gene_aliases_mapping["name"][synonym_upper] = (gene_id_pointer, priority)

                if ensembl_id:
                    gene_aliases_mapping["ensembl_id"][ensembl_id] = gene_id_pointer

                for database, identifier in database_aliases.items():
                    gene_aliases_mapping["databases"].setdefault(database, {})
                    gene_aliases_mapping["databases"][database][identifier] = gene_id_pointer

            # Remove stored priorities after conflict resolution
            gene_aliases_mapping["name"] = {
                alias: pointer
                for alias, (pointer, _) in gene_aliases_mapping["name"].items()
            }

        finally:
            if gi_file_cut.exists():
                gi_file_cut.unlink()

        return gene_aliases_mapping
    
    def __initialize_mappings(self, show_warnings: bool) -> None:
        self.show_warnings = show_warnings
        self.gene_aliases_mapping = self.__parse_ncbi_gene_info(self.ncbi_file)

        self.__upper_gene_names_mapping = {
            key.upper(): value
            for key, value in self.gene_aliases_mapping["name"].items()
        }

        self.databases = set(self.gene_aliases_mapping["databases"].keys())

        self.valid_input_identifier_types = (
            "name",
            "gene_id",
            "ensembl_id",
            *self.databases,
        )

        self.valid_output_identifier_types = (
            "official_name",
            "ncbi_name",
            "gene_id",
            "ensembl_id",
            *self.databases,
        )

    def get_mapping(self):
        return copy.deepcopy(self.gene_aliases_mapping)
    
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

    @support_legacy_gene_synonyms_args
    def get_gene_id(
        self,
        gene: str,
        input_identifier_type: Union[Literal["name", "ensembl_id"], str] = "name"
    ) -> Optional[str]:
        """
        Provide the gene_id with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'ensembl_id' | <database> (default: 'name')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its gene_id.
        """

        if input_identifier_type == "name":
            gene_name = gene.upper()
            if gene_name in self.__upper_gene_names_mapping:
                return self.__upper_gene_names_mapping[gene_name].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no correspondence for gene '{gene}'", stacklevel=10)
                return None
        elif input_identifier_type == "ensembl_id":
            if gene in self.gene_aliases_mapping[input_identifier_type]:
                return self.gene_aliases_mapping[input_identifier_type][gene].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no gene_id correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
                return None
        elif input_identifier_type in self.databases:
            if gene in self.gene_aliases_mapping["databases"][input_identifier_type]:
                return self.gene_aliases_mapping["databases"][input_identifier_type][gene].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(f"no gene_id correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
                return None
        else:
            raise ValueError(f"invalid argument value for 'input_identifier_type': {input_identifier_type}")
    
    @support_legacy_gene_synonyms_args
    def get_ncbi_name(
        self,
        gene: str,
        input_identifier_type: Union[InputIdentifierType, str] = "name"
    ) -> Optional[str]:
        """
        Provide the NCBI reference name with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its NCBI reference name.
        """
    
        gene_id = self.get_gene_id(gene, input_identifier_type) if input_identifier_type != "gene_id" else gene
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].ncbi_name
        else:
            if self.show_warnings:
                warnings.warn(f"no NCBI reference name correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
            return None

    @support_legacy_gene_synonyms_args
    def get_official_name(
        self,
        gene: str,
        input_identifier_type: Union[InputIdentifierType, str] = "name"
    ) -> Optional[str]:
        """
        Provide the official name from nomenclature authority with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its official name from nomenclature authority.
        """
    
        gene_id = self.get_gene_id(gene, input_identifier_type) if input_identifier_type != "gene_id" else gene
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].official_name
        else:
            if self.show_warnings:
                warnings.warn(f"no official name correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
            return None

    @support_legacy_gene_synonyms_args
    def get_ensembl_id(
        self,
        gene: str,
        input_identifier_type: Union[Literal["name", "gene_id"], str] = "name"
    ) -> Optional[str]:
        """
        Provide the ensembl_id with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | <database> (default: 'name')
            Identifier type of the given gene.
            See self.databases for enumerating valid database names.

        Returns
        -------
        Given a gene identifier, return its ensembl_id.
        """
    
        gene_id = self.get_gene_id(gene, input_identifier_type) if input_identifier_type != "gene_id" else gene
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].ensembl_id
        else:
            if self.show_warnings:
                warnings.warn(f"no Ensembl id correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
            return None
    
    @support_legacy_gene_synonyms_args
    def get_alias_from_database(
        self,
        gene: str,
        database: str,
        input_identifier_type: InputIdentifierType = "name"
    ) -> Optional[str]:
        """
        Provide the database-defined gene name with respect to a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        database: <database>
            Organism-related database name providing gene identifiers.
            See self.databases for enumerating valid database names.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' (default: 'name')
            Identifier type of the given gene.

        Returns
        -------
        Given a gene identifier, return its alias derived from a database.
        """

        if database not in self.databases:
            raise ValueError(f"invalid argument value for 'database': got '{database}' but expected a value in {self.databases})")
    
        gene_id = self.get_gene_id(gene, input_identifier_type) if input_identifier_type != "gene_id" else gene
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            if database in self.gene_aliases_mapping["gene_id"][gene_id].databases:
                return self.gene_aliases_mapping["gene_id"][gene_id].databases[database]
            else:
                if self.show_warnings:
                    warnings.warn(f"no {database} correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
                return None
        else:
            if self.show_warnings:
                warnings.warn(f"no {database} correspondence for {input_identifier_type} '{gene}'", stacklevel=10)
            return None

    @support_legacy_gene_synonyms_args
    def conversion(
        self,
        gene: str,
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name"
    ) -> str:
        """
        Convert gene identifiers into the user-defined alias type.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Input identifier type of the given gene.
        output_identifier_type: 'gene_id' | 'official_name' | 'ncbi_name' | 'ensembl_id' | <database> (default: 'official_name')
            Output identifier type for the given gene.

        Returns
        -------
        Given a gene identifier, return its alias.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        if output_identifier_type in ["gene_id", "official_name", "ncbi_name", "ensembl_id"]:
            convert = getattr(self, f"get_{output_identifier_type}")
        elif output_identifier_type in self.databases:
            convert = partial(self.get_alias_from_database, database=output_identifier_type)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute 'get_{output_identifier_type}'")
        
        return convert(gene=gene, input_identifier_type=input_identifier_type)
    
    @support_legacy_gene_synonyms_args
    def __conversion_function(
        self,
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
        *args: Sequence[Any],
        **kwargs: Mapping[str, Any]
    ) -> Callable:
        """
        Function converting gene identifiers.

        Parameters
        ----------
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
            Gene identifier output format.

        Returns
        -------
        Return Function object converting gene identifiers.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        if output_identifier_type in ["gene_id", "official_name", "ncbi_name", "ensembl_id"]:
            return getattr(self, f"get_{output_identifier_type}")
        elif output_identifier_type in self.databases:
            return partial(self.get_alias_from_database, database=output_identifier_type)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute 'get_{output_identifier_type}'")

    @support_legacy_gene_synonyms_args
    def convert_sequence(
        self,
        genes: Sequence[str],
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
        keep_if_missing: bool = True
    ) -> Sequence[str]:
        """
        Create a copy of the Sequence object, with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        genes: Sequence[str]
            List of gene identifiers.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Gene identifier input format.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
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
        
        aliases = list()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for gene in genes:
            output_alias = alias_conversion(gene=gene, input_identifier_type=input_identifier_type)
            output_alias = gene if (keep_if_missing and output_alias is None) else output_alias
            aliases.append(output_alias)
        
        aliases = type(genes)(aliases)
        
        return aliases

    @support_legacy_gene_synonyms_args
    def convert_interaction_list(
        self,
        interaction_list: InteractionList,
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
        keep_if_missing: bool = True
    ) -> InteractionList:
        """
        Create a copy of the pairwise InteractionList object, with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            List of tuples containing string (source) + string (target) + dict (sign = -1 or 1).
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Gene identifier input format.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
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

        converted_interactions_list = list()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for interaction in interaction_list:
            source = alias_conversion(gene=interaction[0], input_identifier_type=input_identifier_type)
            source = interaction[0] if (keep_if_missing and source is None) else source
            target = alias_conversion(gene=interaction[1], input_identifier_type=input_identifier_type)
            target = interaction[1] if (keep_if_missing and target is None) else target
            converted_interactions_list.append((source, target, interaction[2]))

        return converted_interactions_list

    @support_legacy_gene_synonyms_args
    def convert_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Gene identifier input format.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
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
        alias_conversion = self.__conversion_function(output_identifier_type)

        genes = list()

        if axis == 0 or axis == "index":
            iterator = iter(df.index)
        elif axis == 1 or axis == "columns":
            iterator = iter(df.columns)
        else:
            raise TypeError(f"unsupported argument type for 'axis': {axis}")

        for gene in iterator:
            output_alias = alias_conversion(gene=gene, input_identifier_type=input_identifier_type)
            output_alias = gene if output_alias is None else output_alias
            genes.append(output_alias)

        if axis == 0 or axis == "index":
            df.index = genes
        elif axis == 1 or axis == "columns":
            df.columns = genes

        if copy is True:
            return df
    
    @support_legacy_gene_synonyms_args
    def convert_graph(
        self,
        graph: Graph,
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene identifiers in Graph object with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        graph: nx.Graph
            Graph object where gene nodes are converted into the desired identifiers.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Gene identifier input format.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
            Gene identifier output format.
        copy: bool (default: True)
            Return a copy instead of updating Graph object.

        Returns
        -------
        Depending on 'copy', update 'graph' or return Graph object.

        See Also
        --------
        self.get_database() for enumerating valid database names.
        """

        aliases_mapping = dict()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for gene in graph.nodes:
            output_alias = alias_conversion(gene=gene, input_identifier_type=input_identifier_type)
            output_alias = gene if output_alias is None else output_alias
            aliases_mapping[gene] = output_alias
        if copy is True:
            return nx.relabel_nodes(graph, mapping=aliases_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=aliases_mapping, copy=False)
            return None

    @support_legacy_gene_synonyms_args
    def convert_bn(
        self,
        bn: MPBooleanNetwork, # type: ignore
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
        copy: bool = False
    ) -> MPBooleanNetwork: # type: ignore
        """
        Replace gene identifiers in MPBooleanNetwork object with corresponding aliases.
        Each gene identifier is converted into the user-defined alias type.

        Parameters
        ----------
        bn: MPBooleanNetwork
            MPBooleanNetwork object where gene variables are converted into the desired identifiers.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database> (default: 'name')
            Gene identifier input format.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' | 'ensembl_id' | <database> (default: 'official_name')
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

        alias_conversion = self.__conversion_function(output_identifier_type)
        genes = tuple(bn.keys())
        for gene in genes:
            output_alias = alias_conversion(gene=gene, input_identifier_type=input_identifier_type)
            output_alias = gene if output_alias is None else output_alias
            bn.rename(gene, output_alias)
                
        return bn if copy else None

    def standardize_sequence(
        self,
        genes: Sequence[str],
        keep_if_missing: bool = True
    ) -> Sequence[str]:
        """
        Create a copy of the Sequence object, with each gene converted into its official name.

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
            input_identifier_type="name",
            output_identifier_type="official_name",
            keep_if_missing=keep_if_missing
        )

    def standardize_interaction_list(
        self,
        interaction_list: InteractionList,
        keep_if_missing: bool = True
    ) -> InteractionList:
        """
        Create a copy of the pairwise InteractionList object, with each gene converted into its official name.

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
            input_identifier_type="name",
            output_identifier_type="official_name",
            keep_if_missing=keep_if_missing
        )
    
    def standardize_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Replace gene names in DataFrame object with corresponding official names.

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
            input_identifier_type="name",
            output_identifier_type="official_name",
            axis=axis,
            copy=copy
        )

    def standardize_graph(
        self,
        graph: Graph,
        copy: bool = True
    ) -> Union[Graph, None]:
        """
        Replace gene names in Graph object with corresponding official names.

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
            input_identifier_type="name",
            output_identifier_type="official_name",
            copy=copy
        )

    def standardize_bn(
        self,
        bn: MPBooleanNetwork, # type: ignore
        copy: bool = False
    ) -> MPBooleanNetwork: # type: ignore
        """
        Replace gene names in MPBooleanNetwork object with corresponding official names.

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
            input_identifier_type="name",
            output_identifier_type="official_name",
            copy=copy
        )
