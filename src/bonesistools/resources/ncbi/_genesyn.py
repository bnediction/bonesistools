#!/usr/bin/env python

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

if TYPE_CHECKING:
    from ...logic.boolean_network._typing import BooleanNetworkLike

import copy
import ctypes
import gzip
import inspect
import re
import shutil
import tempfile
import urllib.error
import urllib.request
import warnings
from collections import namedtuple
from collections.abc import Mapping as MappingInstance
from collections.abc import Sequence as SequenceInstance
from functools import partial, wraps
from os import PathLike
from pathlib import Path

import networkx as nx
from networkx import Graph
from pandas import DataFrame
from pandas._typing import Axis
from typing_extensions import Literal

from ..._compat import get_args
from ..._validation import _as_boolean, _as_dataframe_axis
from ..._warnings import _warn_deprecated, _warn_deprecated_argument
from ._typing import (
    InputIdentifierType,
)

InteractionList = Sequence[Tuple[str, str, Dict[str, int]]]
GeneInfoVersion = Union[str, PathLike]
_SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING = False

ORGANISMS = Literal[
    "all archaea bacteria",
    "all data",
    "all fungi",
    "all invertebrates",
    "all mammalia",
    "all non mammalian vertebrates",
    "all plants",
    "all protozoa",
    "all viruses",
    "anopheles gambiae",
    "arabidopsis thaliana",
    "archaea",
    "ascomycota",
    "bacteria",
    "bos taurus",
    "caenorhabditis elegans",
    "canis familiaris",
    "chlamydomonas reinhardtii",
    "danio rerio",
    "drosophila melanogaster",
    "escherichia coli",
    "escherichia coli str. k 12 substr. mg1655",
    "gallus gallus",
    "homo sapiens",
    "human",
    "microsporidia",
    "mouse",
    "mus musculus",
    "organelles",
    "oryza sativa",
    "pan troglodytes",
    "penicillium rubens",
    "plasmids",
    "plasmodium falciparum",
    "pseudomonas aeruginosa pao1",
    "rattus norvegicus",
    "retroviridae",
    "saccharomyces cerevisiae",
    "sus scrofa",
    "xenopus laevis",
    "xenopus tropicalis",
    "zea mays",
]

NCBI_GENE_INFO_FTP = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO"

NCBI_GENE_INFO_SOURCES = {
    "all archaea bacteria": (
        "Archaea_Bacteria",
        "All_Archaea_Bacteria.gene_info.gz",
    ),
    "all data": ("", "All_Data.gene_info.gz"),
    "all fungi": ("Fungi", "All_Fungi.gene_info.gz"),
    "all invertebrates": ("Invertebrates", "All_Invertebrates.gene_info.gz"),
    "all mammalia": ("Mammalia", "All_Mammalia.gene_info.gz"),
    "all non mammalian vertebrates": (
        "Non-mammalian_vertebrates",
        "All_Non-mammalian_vertebrates.gene_info.gz",
    ),
    "all plants": ("Plants", "All_Plants.gene_info.gz"),
    "all protozoa": ("Protozoa", "All_Protozoa.gene_info.gz"),
    "all viruses": ("Viruses", "All_Viruses.gene_info.gz"),
    "anopheles gambiae": ("Invertebrates", "Anopheles_gambiae.gene_info.gz"),
    "arabidopsis thaliana": ("Plants", "Arabidopsis_thaliana.gene_info.gz"),
    "archaea": ("Archaea_Bacteria", "Archaea.gene_info.gz"),
    "ascomycota": ("Fungi", "Ascomycota.gene_info.gz"),
    "bacteria": ("Archaea_Bacteria", "Bacteria.gene_info.gz"),
    "bos taurus": ("Mammalia", "Bos_taurus.gene_info.gz"),
    "caenorhabditis elegans": (
        "Invertebrates",
        "Caenorhabditis_elegans.gene_info.gz",
    ),
    "canis familiaris": ("Mammalia", "Canis_familiaris.gene_info.gz"),
    "chlamydomonas reinhardtii": (
        "Plants",
        "Chlamydomonas_reinhardtii.gene_info.gz",
    ),
    "danio rerio": ("Non-mammalian_vertebrates", "Danio_rerio.gene_info.gz"),
    "drosophila melanogaster": (
        "Invertebrates",
        "Drosophila_melanogaster.gene_info.gz",
    ),
    "escherichia coli": (
        "Archaea_Bacteria",
        "Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz",
    ),
    "escherichia coli str. k 12 substr. mg1655": (
        "Archaea_Bacteria",
        "Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz",
    ),
    "gallus gallus": ("Non-mammalian_vertebrates", "Gallus_gallus.gene_info.gz"),
    "homo sapiens": ("Mammalia", "Homo_sapiens.gene_info.gz"),
    "human": ("Mammalia", "Homo_sapiens.gene_info.gz"),
    "microsporidia": ("Fungi", "Microsporidia.gene_info.gz"),
    "mouse": ("Mammalia", "Mus_musculus.gene_info.gz"),
    "mus musculus": ("Mammalia", "Mus_musculus.gene_info.gz"),
    "organelles": ("", "Organelles.gene_info.gz"),
    "oryza sativa": ("Plants", "Oryza_sativa.gene_info.gz"),
    "pan troglodytes": ("Mammalia", "Pan_troglodytes.gene_info.gz"),
    "penicillium rubens": ("Fungi", "Penicillium_rubens.gene_info.gz"),
    "plasmids": ("", "Plasmids.gene_info.gz"),
    "plasmodium falciparum": ("Protozoa", "Plasmodium_falciparum.gene_info.gz"),
    "pseudomonas aeruginosa pao1": (
        "Archaea_Bacteria",
        "Pseudomonas_aeruginosa_PAO1.gene_info.gz",
    ),
    "rattus norvegicus": ("Mammalia", "Rattus_norvegicus.gene_info.gz"),
    "retroviridae": ("Viruses", "Retroviridae.gene_info.gz"),
    "saccharomyces cerevisiae": ("Fungi", "Saccharomyces_cerevisiae.gene_info.gz"),
    "sus scrofa": ("Mammalia", "Sus_scrofa.gene_info.gz"),
    "xenopus laevis": (
        "Non-mammalian_vertebrates",
        "Xenopus_laevis.gene_info.gz",
    ),
    "xenopus tropicalis": (
        "Non-mammalian_vertebrates",
        "Xenopus_tropicalis.gene_info.gz",
    ),
    "zea mays": ("Plants", "Zea_mays.gene_info.gz"),
}

FTP_GENE_INFO = {
    organism: (
        f"{NCBI_GENE_INFO_FTP}/{directory}/{filename}"
        if directory
        else f"{NCBI_GENE_INFO_FTP}/{filename}"
    )
    for organism, (directory, filename) in NCBI_GENE_INFO_SOURCES.items()
}

NCBI_DIR = Path(__file__).resolve().parent / "data" / "gi"

NCBI_GENE_INFO_CACHE_NAMES = {
    "escherichia coli": "escherichia_coli_gene_info.tsv.gz",
    "escherichia coli str. k 12 substr. mg1655": "escherichia_coli_gene_info.tsv.gz",
    "homo sapiens": "homo_sapiens_gene_info.tsv.gz",
    "human": "homo_sapiens_gene_info.tsv.gz",
    "mouse": "mus_musculus_gene_info.tsv.gz",
    "mus musculus": "mus_musculus_gene_info.tsv.gz",
}

NCBI_GENE_INFO_FILES = {
    organism: NCBI_DIR
    / NCBI_GENE_INFO_CACHE_NAMES.get(
        organism,
        f"{re.sub(r'[^a-z0-9]+', '_', organism).strip('_')}_gene_info.tsv.gz",
    )
    for organism in NCBI_GENE_INFO_SOURCES
}

NCBI_GENE_INFO_VERSION_DIR = NCBI_DIR / "versions"

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
    "pseudo": 0,
}

_GENE_SYNONYMS_DEPRECATED_ARGS = {
    "gene_type": "input_identifier_type",
    "alias_gene": "output_identifier_type",
}


def support_legacy_gene_synonyms_args(func):
    """
    Decorate GeneSynonyms methods to accept deprecated argument names.

    """

    valid_parameters = inspect.signature(func).parameters

    @wraps(func)
    def wrapper(*args, **kwargs):
        for old_name, new_name in _GENE_SYNONYMS_DEPRECATED_ARGS.items():

            if old_name not in kwargs:
                continue

            if new_name not in valid_parameters:
                continue

            _warn_deprecated_argument(old_name, new_name, stacklevel=2)

            if new_name in kwargs:
                raise TypeError(
                    f"invalid argument combination: use either '{old_name}' "
                    f"or '{new_name}', not both"
                )

            kwargs[new_name] = kwargs.pop(old_name)

        return func(*args, **kwargs)

    return wrapper


class GeneSynonyms:
    """
    Converter between gene identifiers using NCBI gene_info resources.

    GeneSynonyms loads organism-specific NCBI gene information and builds
    mappings between gene IDs, official symbols, NCBI names, Ensembl IDs and
    database-specific aliases when available.

    The instance is callable and dispatches conversion according to the input
    object type:

    - sequences of gene identifiers are converted with `convert_sequence`,
    - interaction lists are converted with `convert_interaction_list`,
    - pandas DataFrames are converted with `convert_df`,
    - NetworkX graphs are converted with `convert_graph`,
    - BooleanNetworkLike objects are converted with `convert_bn`.

    Examples
    --------
    Convert a list of gene identifiers through the callable interface:

    >>> gene_synonyms = genesyn(organism="mouse")
    >>> gene_synonyms(["Trp53", "Myc"])
    ['Trp53', 'Myc']

    Convert an interaction list while preserving edge attributes:

    >>> interactions = [("Trp53", "Myc", {"sign": 1})]
    >>> gene_synonyms(interactions)
    [('Trp53', 'Myc', {'sign': 1})]

    Parameters
    ----------
    organism: str (default: "mouse")
        Common name of the organism of interest. Supported organisms match the
        NCBI gene_info files listed in `FTP_GENE_INFO`. Common aliases such as
        `"mouse"`, `"human"` and `"escherichia coli"` are accepted.
    version: "bundled", "latest", date string, str path or PathLike (default: "bundled")
        NCBI gene_info version to load.

        If `"bundled"`, use the local organism-specific gene_info file if
        available; if it is missing, download the latest NCBI gene_info file
        before building mappings.

        If `"latest"`, download the latest NCBI gene_info file before building
        mappings.

        If a date such as `"2024-01-01"` or `"20240101"` is provided, load the
        matching local version from `data/gi/versions/<date>/` or
        `data/gi/<organism>_gene_info_<date>.tsv`.

        If a path-like object or path string is provided, load that gene_info
        file directly.
    show_warnings: bool (default: False)
        If `True`, warn when a requested gene identifier has no correspondence.

    Attributes
    ----------
    databases: set
        Database names available for the selected organism, such as `MGI` when
        present in the NCBI gene_info resource.
    valid_input_identifier_types: tuple
        Identifier types accepted as input.
    valid_output_identifier_types: tuple
        Identifier types accepted as output.

    Raises
    ------
    ValueError
        If `organism`, identifier types, database names or axis values are
        unsupported.
    AttributeError
        If an output identifier type has no corresponding conversion method.
    RuntimeError
        If downloading or parsing NCBI gene_info data fails.
    """

    def __init__(
        self,
        organism: str = "mouse",
        version: GeneInfoVersion = "bundled",
        show_warnings: bool = False,
    ) -> None:

        if not _SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING:
            _warn_deprecated(
                "`bt.resources.ncbi.GeneSynonyms()`",
                replacement="`bt.resources.ncbi.genesyn()`",
                stacklevel=2,
            )

        organism = organism.lower().replace("-", " ")

        if organism not in get_args(ORGANISMS):
            raise ValueError(
                f"invalid argument value for 'organism': "
                f"expected one of {get_args(ORGANISMS)} but received {organism!r}"
            )

        show_warnings = _as_boolean(show_warnings, "show_warnings")

        self.organism = organism
        self.version = self.__normalize_gene_info_version(version)
        self.ncbi_file = self.__resolve_gene_info_file(self.version)

        self.__download_gene_info()
        self._initialize_mappings(show_warnings=show_warnings)

    def __call__(
        self,
        data: Union[InteractionList, DataFrame, Graph[Any], "BooleanNetworkLike"],
        *args: Any,
        **kwargs: Any,
    ):
        """
        Convert gene identifiers in a supported object.

        This method dispatches to the appropriate conversion method according
        to `data` type. Additional positional and keyword arguments are passed
        to the selected conversion method.

        Parameters
        ----------
        data: sequence, InteractionList, DataFrame, Graph or BooleanNetworkLike
            Object containing gene identifiers to convert.
        *args: Any
            Positional arguments forwarded to the selected conversion method.
        **kwargs: Any
            Keyword arguments forwarded to the selected conversion method.

        Returns
        -------
        object
            Converted object. The exact type depends on `data` and the selected
            conversion method.

        """

        from ...logic.boolean_network._typing import is_boolean_network_like

        if _is_interaction_list(data):
            return self.convert_interaction_list(
                cast(InteractionList, data),
                *args,
                **kwargs,
            )
        elif (
            isinstance(data, SequenceInstance) and not isinstance(data, str)
        ) or isinstance(data, set):
            return self.convert_sequence(cast(Sequence[str], data), *args, **kwargs)
        elif isinstance(data, DataFrame):
            return self.convert_df(data, *args, **kwargs)
        elif isinstance(data, Graph):
            return self.convert_graph(data, *args, **kwargs)
        elif is_boolean_network_like(data):
            return self.convert_bn(data, *args, **kwargs)
        else:
            raise TypeError(
                f"unsupported argument type for 'data': "
                f"expected sequence, {DataFrame}, {Graph} or Boolean network-like "
                f"object but received {type(data)}"
            )

    @support_legacy_gene_synonyms_args
    def conversion(
        self,
        gene: str,
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
    ) -> Optional[str]:
        """
        Convert gene identifiers into the user-defined alias type.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        str or None
            Converted gene identifier, or None if no match is found.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        if output_identifier_type in [
            "gene_id",
            "official_name",
            "ncbi_name",
            "ensembl_id",
        ]:
            convert = getattr(self, f"get_{output_identifier_type}")
        elif output_identifier_type in self.databases:
            convert = partial(
                self.get_alias_from_database, database=output_identifier_type
            )
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute "
                f"'get_{output_identifier_type}'"
            )

        return convert(
            gene=gene,
            input_identifier_type=cast(InputIdentifierType, input_identifier_type),
        )

    @support_legacy_gene_synonyms_args
    def convert_sequence(
        self,
        genes: Sequence[str],
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        keep_if_missing: bool = True,
    ) -> Sequence[str]:
        """
        Convert a sequence of gene identifiers.

        Each gene identifier is converted into the requested output identifier
        type.

        Parameters
        ----------
        genes: Sequence[str]
            Gene identifiers to convert.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.
        keep_if_missing: bool (default: True)
            If `True`, keep the original gene identifier when no correspondence
            is found.

        Returns
        -------
        Sequence[str]
            Converted sequence, preserving the input sequence type when
            possible.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        aliases = list()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for gene in genes:
            output_alias = alias_conversion(
                gene=gene, input_identifier_type=input_identifier_type
            )
            output_alias = (
                gene if (keep_if_missing and output_alias is None) else output_alias
            )
            aliases.append(output_alias)

        sequence_constructor = cast(
            Callable[[Sequence[str]], Sequence[str]], type(genes)
        )
        aliases = sequence_constructor(aliases)

        return aliases

    @support_legacy_gene_synonyms_args
    def convert_interaction_list(
        self,
        interaction_list: InteractionList,
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        keep_if_missing: bool = True,
    ) -> InteractionList:
        """
        Convert source and target identifiers in an interaction list.

        Edge attributes are preserved.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            Sequence of `(source, target, attributes)` interactions.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.
        keep_if_missing: bool (default: True)
            If `True`, keep the original gene identifier when no correspondence
            is found.

        Returns
        -------
        InteractionList
            Converted interaction list with source and target identifiers
            replaced and edge attributes preserved.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        converted_interactions_list = list()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for interaction in interaction_list:
            source = alias_conversion(
                gene=interaction[0], input_identifier_type=input_identifier_type
            )
            source = interaction[0] if (keep_if_missing and source is None) else source
            target = alias_conversion(
                gene=interaction[1], input_identifier_type=input_identifier_type
            )
            target = interaction[1] if (keep_if_missing and target is None) else target
            converted_interactions_list.append((source, target, interaction[2]))

        return converted_interactions_list

    @support_legacy_gene_synonyms_args
    def convert_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Convert gene identifiers in a DataFrame index or columns.

        Parameters
        ----------
        df: pd.DataFrame
            DataFrame whose gene identifiers are converted.
        axis: {0, 1, "index", "columns"} (default: 0)
            If 0 or `"index"`, convert `df.index`. If 1 or `"columns"`,
            convert `df.columns`.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.
        copy: bool (default: True)
            Return a copy instead of modifying `df`.

        Returns
        -------
        DataFrame or None
            Converted DataFrame if `copy=True`; otherwise None.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        df = df.copy() if copy is True else df
        alias_conversion = self.__conversion_function(output_identifier_type)

        genes = list()
        axis = _as_dataframe_axis(axis)

        if axis == "index":
            iterator = iter(df.index)
        elif axis == "columns":
            iterator = iter(df.columns)

        for gene in iterator:
            output_alias = alias_conversion(
                gene=gene, input_identifier_type=input_identifier_type
            )
            output_alias = gene if output_alias is None else output_alias
            genes.append(output_alias)

        if axis == "index":
            df.index = genes
        elif axis == "columns":
            df.columns = genes

        if copy is True:
            return df

    @support_legacy_gene_synonyms_args
    def convert_graph(
        self,
        graph: Graph[Any],
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        copy: bool = True,
    ) -> Union[Graph[Any], None]:
        """
        Convert gene identifiers in graph node labels.

        Parameters
        ----------
        graph: nx.Graph
            Graph whose node identifiers are converted.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.
        copy: bool (default: True)
            Return a copy instead of modifying `graph`.

        Returns
        -------
        Graph or None
            Converted graph if `copy=True`; otherwise None.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        aliases_mapping = dict()
        alias_conversion = self.__conversion_function(output_identifier_type)
        for gene in graph.nodes:
            output_alias = alias_conversion(
                gene=gene, input_identifier_type=input_identifier_type
            )
            output_alias = gene if output_alias is None else output_alias
            aliases_mapping[gene] = output_alias

        from ...logic.influence_graph import InfluenceGraph

        if type(graph) is InfluenceGraph:
            converted_graph = graph.copy()
            converted_graph.relabel(aliases_mapping)
            if copy is True:
                return converted_graph
            graph._replace_with_graph(converted_graph)
            return None

        if copy is True:
            return nx.relabel_nodes(graph, mapping=aliases_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=aliases_mapping, copy=False)
            return None

    @support_legacy_gene_synonyms_args
    def convert_bn(
        self,
        bn: "BooleanNetworkLike",
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        copy: bool = False,
    ) -> Union["BooleanNetworkLike", None]:
        """
        Convert gene identifiers in Boolean network components and rules.

        Parameters
        ----------
        bn: BooleanNetworkLike
            BooleanNetworkLike object whose component identifiers are converted.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.
        copy: bool (default: False)
            Return a copy instead of modifying `bn`.

        Returns
        -------
        BooleanNetworkLike or None
            Converted BooleanNetworkLike object if `copy=True`; otherwise None.
            If `copy=True` and `bn` does not provide a `rename` method, a
            temporary `BooleanNetwork` is used to perform the rule-safe
            renaming, then the result is converted back to the input type.

        Notes
        -----
        Valid identifier types are available from
        `valid_input_identifier_types` and `valid_output_identifier_types`.
        """

        from ...logic.boolean_network import BooleanNetwork
        from ...logic.boolean_network._typing import is_boolean_network_like

        if not is_boolean_network_like(bn):
            raise TypeError(
                "unsupported argument type for 'bn': "
                "expected Boolean network-like object "
                f"but received {type(bn)}"
            )

        input_type = type(bn)
        bn_any: Any = bn
        rebuild_as_input_type = False
        rename = getattr(bn_any, "rename", None)

        if copy:
            copy_method = getattr(bn_any, "copy", None)
            bn_any = (
                copy_method()
                if callable(copy_method)
                else BooleanNetwork(bn_any, check=False)
            )
            rename = getattr(bn_any, "rename", None)
            if not callable(rename):
                bn_any = BooleanNetwork(bn_any, check=False)
                rename = bn_any.rename
                rebuild_as_input_type = True
        elif not callable(rename):
            raise TypeError(
                "unsupported argument value for 'bn': in-place Boolean network "
                "conversion requires a 'rename' method; use copy=True to return "
                "a converted copy"
            )

        alias_conversion = self.__conversion_function(output_identifier_type)
        genes = tuple(gene for gene, _ in bn_any.items())
        for gene in genes:
            output_alias = alias_conversion(
                gene=gene, input_identifier_type=input_identifier_type
            )
            output_alias = gene if output_alias is None else output_alias
            rename(gene, output_alias)

        if copy and rebuild_as_input_type:
            try:
                return cast(
                    "BooleanNetworkLike",
                    cast(Callable[[Any], Any], input_type)(bn_any.rules),
                )
            except Exception as e:
                raise TypeError(
                    "unable to return converted Boolean network with original "
                    f"type {input_type}"
                ) from e

        return cast("BooleanNetworkLike", bn_any) if copy else None

    def standardize_sequence(
        self, genes: Sequence[str], keep_if_missing: bool = True
    ) -> Sequence[str]:
        """
        Standardize a sequence of gene names to official gene names.

        Parameters
        ----------
        genes: Sequence[str]
            Gene names to standardize.
        keep_if_missing: bool (default: True)
            If `True`, keep the original gene identifier when no correspondence
            is found.

        Returns
        -------
        Sequence[str]
            Sequence of official gene names, preserving the input sequence type
            when possible.
        """

        return self.convert_sequence(
            genes=genes,
            input_identifier_type="name",
            output_identifier_type="official_name",
            keep_if_missing=keep_if_missing,
        )

    def standardize_interaction_list(
        self, interaction_list: InteractionList, keep_if_missing: bool = True
    ) -> InteractionList:
        """
        Standardize source and target identifiers in an interaction list.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            Sequence of `(source, target, attributes)` interactions.
        keep_if_missing: bool (default: True)
            If `True`, keep the original gene identifier when no correspondence
            is found.

        Returns
        -------
        InteractionList
            Interaction list with source and target identifiers standardized to
            official gene names.
        """

        return self.convert_interaction_list(
            interaction_list=interaction_list,
            input_identifier_type="name",
            output_identifier_type="official_name",
            keep_if_missing=keep_if_missing,
        )

    def standardize_df(
        self,
        df: DataFrame,
        axis: Axis = 0,
        copy: bool = True,
    ) -> Union[DataFrame, None]:
        """
        Standardize gene names in a DataFrame index or columns.

        Parameters
        ----------
        df: pd.DataFrame
            DataFrame whose gene identifiers are standardized.
        axis: {0, 1, "index", "columns"} (default: 0)
            If 0 or `"index"`, standardize `df.index`. If 1 or `"columns"`,
            standardize `df.columns`.
        copy: bool (default: True)
            Return a copy instead of modifying `df`.

        Returns
        -------
        DataFrame or None
            Standardized DataFrame if `copy=True`; otherwise None.
        """

        return self.convert_df(
            df=df,
            input_identifier_type="name",
            output_identifier_type="official_name",
            axis=axis,
            copy=copy,
        )

    def standardize_graph(
        self,
        graph: Graph[Any],
        copy: bool = True,
    ) -> Union[Graph[Any], None]:
        """
        Standardize gene names in graph node labels.

        Parameters
        ----------
        graph: nx.Graph
            Graph whose node identifiers are standardized.
        copy: bool (default: True)
            Return a copy instead of modifying `graph`.

        Returns
        -------
        Graph or None
            Standardized graph if `copy=True`; otherwise None.
        """

        return self.convert_graph(
            graph=graph,
            input_identifier_type="name",
            output_identifier_type="official_name",
            copy=copy,
        )

    def standardize_bn(
        self, bn: "BooleanNetworkLike", copy: bool = False
    ) -> Union["BooleanNetworkLike", None]:
        """
        Standardize gene names in Boolean network components and rules.

        Parameters
        ----------
        bn: BooleanNetworkLike
            BooleanNetworkLike object whose component identifiers are standardized.
        copy: bool (default: False)
            Return a copy instead of modifying `bn`.

        Returns
        -------
        BooleanNetworkLike or None
            Standardized BooleanNetworkLike object if `copy=True`; otherwise None.
        """

        return self.convert_bn(
            bn=bn,
            input_identifier_type="name",
            output_identifier_type="official_name",
            copy=copy,
        )

    def reset(
        self,
        organism: Optional[str] = None,
        version: Optional[GeneInfoVersion] = None,
        show_warnings: bool = False,
    ) -> None:
        """
        Reinitialize the converter with a potentially different organism.

        Parameters
        ----------
        organism: str, optional
            Common name of the organism to load. If `None`, reload the current
            organism.
        version: "bundled", "latest", date string, str path or PathLike, optional
            NCBI gene_info version to load. If `None`, keep the current version.
        show_warnings: bool (default: False)
            If `True`, warn when a requested gene identifier has no correspondence.

        Raises
        ------
        ValueError
            If `organism` is unsupported.
        RuntimeError
            If downloading or parsing NCBI gene_info data fails.
        """

        if organism is None:
            organism = self.organism
        else:
            organism = organism.lower().replace("-", " ")

        if organism not in get_args(ORGANISMS):
            raise ValueError(
                f"invalid argument value for 'organism': "
                f"expected one of {get_args(ORGANISMS)} but received {organism!r}"
            )

        show_warnings = _as_boolean(show_warnings, "show_warnings")

        self.organism = organism
        resolved_version: GeneInfoVersion
        if version is None:
            resolved_version = cast(
                GeneInfoVersion,
                getattr(self, "version", "bundled"),
            )
        else:
            resolved_version = version
        self.version = self.__normalize_gene_info_version(resolved_version)
        self.ncbi_file = self.__resolve_gene_info_file(self.version)

        self.__download_gene_info()
        self._initialize_mappings(show_warnings=show_warnings)

    def get_mapping(self):
        """
        Return a deep copy of the internal gene alias mapping.

        Returns
        -------
        Dict
            Copy of the mapping structure used for identifier conversion.
        """

        return copy.deepcopy(self.gene_aliases_mapping)

    def contains(
        self,
        *genes: str,
        identifier_type: Union[Literal["name", "gene_id", "ensembl_id"], str] = "name",
    ) -> List[bool]:
        """
        Test whether gene identifiers are present in the synonym resource.

        Examples
        --------
        >>> gene_synonyms = genesyn(organism="mouse")
        >>> gene_synonyms.contains("Trp53")
        [True]
        >>> gene_synonyms.contains("Trp53", "not-a-gene")
        [True, False]
        >>> gene_synonyms.contains("22059", identifier_type="gene_id")
        [True]

        Parameters
        ----------
        *genes: str
            Gene identifiers to test.
        identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Identifier type used to interpret `genes`.

        Returns
        -------
        list of bool
            Boolean values ordered like `genes`. An empty input returns an
            empty list.

        Raises
        ------
        ValueError
            If `identifier_type` is not supported.
        """

        if identifier_type == "name":
            return [gene.upper() in self.__upper_gene_names_mapping for gene in genes]

        if identifier_type == "gene_id":
            return [gene in self.gene_aliases_mapping["gene_id"] for gene in genes]

        if identifier_type == "ensembl_id":
            return [gene in self.gene_aliases_mapping["ensembl_id"] for gene in genes]

        if identifier_type in self.databases:
            database_mapping = self.gene_aliases_mapping["databases"][identifier_type]
            return [gene in database_mapping for gene in genes]

        raise ValueError(
            f"invalid argument value for 'identifier_type': "
            f"expected one of {self.valid_input_identifier_types} "
            f"but received {identifier_type!r}"
        )

    def find(
        self,
        *genes: str,
        identifier_type: Union[Literal["name", "gene_id", "ensembl_id"], str] = "name",
    ) -> List[str]:
        """
        Return gene identifiers found in the synonym resource.

        Examples
        --------
        >>> gene_synonyms = genesyn(organism="mouse")
        >>> gene_synonyms.find("Trp53", "not-a-gene", "Myc")
        ['Trp53', 'Myc']
        >>> gene_synonyms.find("22059", "bad-id", identifier_type="gene_id")
        ['22059']

        Parameters
        ----------
        *genes: str
            Gene identifiers to test.
        identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Identifier type used to interpret `genes`.

        Returns
        -------
        list of str
            Input gene identifiers present in the synonym resource, preserving
            their input order.

        Raises
        ------
        ValueError
            If `identifier_type` is not supported.
        """

        return [
            gene
            for gene, found in zip(
                genes,
                self.contains(*genes, identifier_type=identifier_type),
            )
            if found
        ]

    @support_legacy_gene_synonyms_args
    def get_gene_id(
        self,
        gene: str,
        input_identifier_type: Union[Literal["name", "ensembl_id"], str] = "name",
    ) -> Optional[str]:
        """
        Return the NCBI gene ID for a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'ensembl_id' | <database> (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        str or None
            NCBI gene ID corresponding to `gene`, or None if no match is
            found.
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
                return self.gene_aliases_mapping[input_identifier_type][
                    gene
                ].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(
                        f"no gene_id correspondence found for {gene!r} "
                        f"using input_identifier_type={input_identifier_type!r}",
                        stacklevel=10,
                    )
                return None
        elif input_identifier_type in self.databases:
            if gene in self.gene_aliases_mapping["databases"][input_identifier_type]:
                return self.gene_aliases_mapping["databases"][input_identifier_type][
                    gene
                ].value.decode()
            else:
                if self.show_warnings:
                    warnings.warn(
                        f"no gene_id correspondence found for {gene!r} "
                        f"using input_identifier_type={input_identifier_type!r}",
                        stacklevel=10,
                    )
                return None
        else:
            raise ValueError(
                f"invalid argument value for 'input_identifier_type': "
                f"expected 'name', 'ensembl_id' or one of {self.databases} "
                f"but received {input_identifier_type!r}"
            )

    @support_legacy_gene_synonyms_args
    def get_ncbi_name(
        self,
        gene: str,
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
    ) -> Optional[str]:
        """
        Return the NCBI reference name for a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        str or None
            NCBI reference name corresponding to `gene`, or None if no match is
            found.
        """

        gene_id = (
            self.get_gene_id(gene, input_identifier_type)
            if input_identifier_type != "gene_id"
            else gene
        )
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].ncbi_name
        else:
            if self.show_warnings:
                warnings.warn(
                    f"no NCBI reference name correspondence found for {gene!r} "
                    f"using input_identifier_type={input_identifier_type!r}",
                    stacklevel=10,
                )
            return None

    @support_legacy_gene_synonyms_args
    def get_official_name(
        self,
        gene: str,
        input_identifier_type: Union[
            Literal["name", "gene_id", "ensembl_id"], str
        ] = "name",
    ) -> Optional[str]:
        """
        Return the official nomenclature name for a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' | <database>
            (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        str or None
            Official gene name corresponding to `gene`, or None if no match is
            found.
        """

        gene_id = (
            self.get_gene_id(gene, input_identifier_type)
            if input_identifier_type != "gene_id"
            else gene
        )
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].official_name
        else:
            if self.show_warnings:
                warnings.warn(
                    f"no official name correspondence found for {gene!r} "
                    f"using input_identifier_type={input_identifier_type!r}",
                    stacklevel=10,
                )
            return None

    @support_legacy_gene_synonyms_args
    def get_ensembl_id(
        self,
        gene: str,
        input_identifier_type: Union[Literal["name", "gene_id"], str] = "name",
    ) -> Optional[str]:
        """
        Return the Ensembl ID for a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        input_identifier_type: 'name' | 'gene_id' | <database> (default: 'name')
            Input gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        str or None
            Ensembl ID corresponding to `gene`, or None if no match is found.
        """

        gene_id = (
            self.get_gene_id(gene, input_identifier_type)
            if input_identifier_type != "gene_id"
            else gene
        )
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            return self.gene_aliases_mapping["gene_id"][gene_id].ensembl_id
        else:
            if self.show_warnings:
                warnings.warn(
                    f"no Ensembl id correspondence found for {gene!r} "
                    f"using input_identifier_type={input_identifier_type!r}",
                    stacklevel=10,
                )
            return None

    @support_legacy_gene_synonyms_args
    def get_alias_from_database(
        self,
        gene: str,
        database: str,
        input_identifier_type: InputIdentifierType = "name",
    ) -> Optional[str]:
        """
        Return a database-defined gene alias for a gene identifier.

        Parameters
        ----------
        gene: str
            Identifier of the gene of interest.
        database: <database>
            Organism-related database name providing gene identifiers.
            See `databases` for valid database names.
        input_identifier_type: 'name' | 'gene_id' | 'ensembl_id' (default: 'name')
            Input gene identifier type.

        Returns
        -------
        str or None
            Database-specific alias corresponding to `gene`, or None if no
            match is found.
        """

        if database not in self.databases:
            raise ValueError(
                f"invalid argument value for 'database': "
                f"expected one of {self.databases} but received {database!r}"
            )

        gene_id = (
            self.get_gene_id(gene, input_identifier_type)
            if input_identifier_type != "gene_id"
            else gene
        )
        if gene_id in self.gene_aliases_mapping["gene_id"]:
            if database in self.gene_aliases_mapping["gene_id"][gene_id].databases:
                return self.gene_aliases_mapping["gene_id"][gene_id].databases[database]
            else:
                if self.show_warnings:
                    warnings.warn(
                        f"no {database} correspondence found for {gene!r} "
                        f"using input_identifier_type={input_identifier_type!r}",
                        stacklevel=10,
                    )
                return None
        else:
            if self.show_warnings:
                warnings.warn(
                    f"no {database} correspondence found for {gene!r} "
                    f"using input_identifier_type={input_identifier_type!r}",
                    stacklevel=10,
                )
            return None

    def __normalize_gene_info_version(self, version: GeneInfoVersion) -> str:

        if not isinstance(version, (str, PathLike)):
            raise TypeError(
                f"unsupported argument type for 'version': "
                f"expected {str} or {PathLike} but received {type(version)}"
            )

        if isinstance(version, PathLike):
            return str(version)

        version = version.strip()

        if version in ["bundled", ""]:
            return "bundled"

        if version == "latest":
            return "latest"

        date_match = re.fullmatch(r"([0-9]{4})-?([0-9]{2})-?([0-9]{2})", version)
        if date_match:
            return "".join(date_match.groups())

        return version

    def __resolve_gene_info_file(self, version: str) -> Path:

        bundled_file = NCBI_GENE_INFO_FILES[self.organism]

        if version == "bundled":
            return bundled_file

        if version == "latest":
            latest_file = tempfile.NamedTemporaryFile(
                prefix=f"{self.__gene_info_base_name(bundled_file)}_",
                suffix=".tsv",
                delete=False,
            )
            latest_path = Path(latest_file.name)
            latest_file.close()
            self.__unlink_if_exists(latest_path)
            return latest_path

        version_path = Path(version)
        if version_path.exists():
            return version_path

        if re.fullmatch(r"[0-9]{8}", version):
            base_name = self.__gene_info_base_name(bundled_file)
            dated_candidates = [
                NCBI_GENE_INFO_VERSION_DIR / version / bundled_file.name,
                NCBI_GENE_INFO_VERSION_DIR / version / f"{base_name}.tsv",
                NCBI_DIR / f"{base_name}_{version}.tsv.gz",
                NCBI_DIR / f"{base_name}_{version}.tsv",
            ]
            for candidate in dated_candidates:
                if candidate.exists():
                    return candidate

            candidates = ", ".join(str(candidate) for candidate in dated_candidates)
            raise FileNotFoundError(
                f"NCBI gene_info version not found for {self.organism!r} "
                f"at {version}: expected one of {candidates}"
            )

        raise FileNotFoundError(f"NCBI gene_info version file not found: {version}")

    def __download_gene_info(self) -> None:

        if self.version == "bundled" and self.ncbi_file.exists():
            return

        if self.version == "bundled":
            full_file = self.__make_temporary_gene_info_file()
            try:
                self.__download_full_gene_info(full_file)
                self.__write_bundled_gene_info_file(full_file, self.ncbi_file)
            finally:
                self.__unlink_if_exists(full_file)
                self.__unlink_if_exists(Path(f"{full_file}.gz"))
            return

        if self.version == "latest":
            try:
                self.__download_full_gene_info(self.ncbi_file)
            except RuntimeError:
                self.__unlink_if_exists(self.ncbi_file)
                self.__unlink_if_exists(Path(f"{self.ncbi_file}.gz"))
                raise
            return

    def __download_full_gene_info(self, outfile: Path) -> None:

        url = FTP_GENE_INFO[self.organism]
        outfile.parent.mkdir(parents=True, exist_ok=True)
        temporary_gzip = Path(f"{outfile}.gz.tmp")
        temporary_outfile = Path(f"{outfile}.tmp")

        try:
            with urllib.request.urlopen(url) as response:
                with open(temporary_gzip, "wb") as file:
                    shutil.copyfileobj(response, file)

            with gzip.open(temporary_gzip, "rb") as compressed_file:
                with open(temporary_outfile, "wb") as file:
                    shutil.copyfileobj(compressed_file, file)

            temporary_outfile.replace(outfile)
        except (OSError, urllib.error.URLError) as e:
            raise RuntimeError("failed to download NCBI gene_info file") from e
        finally:
            self.__unlink_if_exists(temporary_gzip)
            self.__unlink_if_exists(temporary_outfile)

    def __write_bundled_gene_info_file(self, infile: Path, outfile: Path) -> None:

        outfile.parent.mkdir(parents=True, exist_ok=True)
        temporary_outfile = Path(f"{outfile}.tmp")

        try:
            with self.__open_gene_info_output(outfile, temporary_outfile) as file:
                file.write(
                    "\t".join(
                        [
                            "gene_id",
                            "official_name",
                            "ncbi_name",
                            "synonyms",
                            "dbXrefs",
                            "chromosome",
                            "gene_type",
                        ]
                    )
                    + "\n"
                )
                for fields in self.__iter_gene_info_rows(infile, reduced=False):
                    file.write("\t".join(fields) + "\n")

            temporary_outfile.replace(outfile)
        except (IndexError, OSError) as e:
            raise RuntimeError("failed to reduce NCBI gene_info file") from e
        finally:
            self.__unlink_if_exists(temporary_outfile)

    def __parse_ncbi_gene_info(self, gi_file: Path) -> Dict[str, Any]:
        """
        Parse NCBI gene_info into mappings between gene identifiers.

        Parameters
        ----------
        gi_file: Path
            Path to the NCBI gene_info data.

        Returns
        -------
        Dict
            Mapping structure for gene identifiers and aliases.
        """

        Identifiers = namedtuple(
            "Identifiers",
            [
                "official_name",
                "ncbi_name",
                "ensembl_id",
                "databases",
                "chromosome",
                "gene_type",
            ],
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
            for fields in self.__iter_gene_info_rows(
                gi_file,
                reduced=self.version == "bundled",
            ):
                gene_id = fields[0]
                official_name = fields[1]
                ncbi_name = fields[2]
                chromosome = fields[5]
                gene_type = fields[6]
                if official_name == "-":
                    official_name = ncbi_name

                gene_info_rows.append(
                    (
                        gene_id,
                        official_name,
                        ncbi_name,
                        fields[3],
                        fields[4],
                        chromosome,
                        gene_type,
                    )
                )

                if GENE_TYPE_PRIORITY.get(gene_type, -1) > 0:
                    protected_official_names.add(official_name.upper())

            for (
                gene_id,
                official_name,
                ncbi_name,
                synonyms_field,
                db_xrefs_field,
                chromosome,
                gene_type,
            ) in gene_info_rows:

                synonyms = [
                    synonym
                    for synonym in (synonyms_field.split("|") + [ncbi_name])
                    if synonym != "-" and synonym.upper() != official_name.upper()
                ]

                ensembl_id = None
                database_aliases = {}

                for db_entry in db_xrefs_field.split("|"):
                    if db_entry == "-":
                        continue

                    database, database_name = db_entry.split(":", maxsplit=1)

                    if database == "Ensembl":
                        ensembl_match = re.findall("[A-Z]{7}[0-9]{11}", database_name)
                        ensembl_id = ensembl_match[0] if ensembl_match else None
                    else:
                        database_aliases[database] = database_name

                gene_id_pointer = ctypes.create_string_buffer(gene_id.encode())

                gene_aliases_mapping["gene_id"][gene_id] = Identifiers(
                    official_name=official_name,
                    ncbi_name=ncbi_name,
                    ensembl_id=ensembl_id,
                    databases=database_aliases,
                    chromosome=chromosome,
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
                        gene_aliases_mapping["name"][synonym_upper] = (
                            gene_id_pointer,
                            priority,
                        )
                    else:
                        _, old_priority = gene_aliases_mapping["name"][synonym_upper]
                        if priority > old_priority:
                            gene_aliases_mapping["name"][synonym_upper] = (
                                gene_id_pointer,
                                priority,
                            )

                if ensembl_id:
                    gene_aliases_mapping["ensembl_id"][ensembl_id] = gene_id_pointer

                for database, identifier in database_aliases.items():
                    gene_aliases_mapping["databases"].setdefault(database, {})
                    gene_aliases_mapping["databases"][database][
                        identifier
                    ] = gene_id_pointer

            # Remove stored priorities after conflict resolution
            gene_aliases_mapping["name"] = {
                alias: pointer
                for alias, (pointer, _) in gene_aliases_mapping["name"].items()
            }
        except (IndexError, OSError) as e:
            raise RuntimeError("failed to parse NCBI gene_info file") from e

        return gene_aliases_mapping

    def _initialize_mappings(self, show_warnings: bool) -> None:

        self.show_warnings = show_warnings
        try:
            self.gene_aliases_mapping = self.__parse_ncbi_gene_info(self.ncbi_file)
        finally:
            if self.version == "latest":
                self.__unlink_if_exists(self.ncbi_file)
                self.__unlink_if_exists(Path(f"{self.ncbi_file}.gz"))

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

    def __initialize_mappings(self, show_warnings: bool) -> None:

        self._initialize_mappings(show_warnings=show_warnings)

    @support_legacy_gene_synonyms_args
    def __conversion_function(
        self,
        output_identifier_type: Union[
            Literal["official_name", "ncbi_name", "gene_id", "ensembl_id"], str
        ] = "official_name",
        *args: Any,
        **kwargs: Any,
    ) -> Callable[..., Optional[str]]:
        """
        Function converting gene identifiers.

        Parameters
        ----------
        output_identifier_type: 'official_name' | 'ncbi_name' | 'gene_id' |
            'ensembl_id' | <database> (default: 'official_name')
            Output gene identifier type. Valid database-specific values are
            listed in `databases`.

        Returns
        -------
        Callable
            Function converting one gene identifier to the requested output
            identifier type.

        Notes
        -----
        Valid output identifier types are available from
        `valid_output_identifier_types`.
        """

        if output_identifier_type in [
            "gene_id",
            "official_name",
            "ncbi_name",
            "ensembl_id",
        ]:
            return getattr(self, f"get_{output_identifier_type}")
        elif output_identifier_type in self.databases:
            return partial(
                self.get_alias_from_database, database=output_identifier_type
            )
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute "
                f"'get_{output_identifier_type}'"
            )

    @staticmethod
    def __gene_info_base_name(gene_info_file: Path) -> str:

        name = gene_info_file.name
        if name.endswith(".tsv.gz"):
            return name[:-7]
        if name.endswith(".tsv"):
            return name[:-4]
        return gene_info_file.stem

    @staticmethod
    def __make_temporary_gene_info_file() -> Path:

        temporary_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)
        temporary_file.close()
        return Path(temporary_file.name)

    @staticmethod
    def __open_gene_info_input(gene_info_file: Path) -> Any:

        if str(gene_info_file).endswith(".gz"):
            return gzip.open(gene_info_file, "rt", encoding="utf-8")
        return open(gene_info_file, "r", encoding="utf-8")

    @staticmethod
    def __open_gene_info_output(outfile: Path, temporary_outfile: Path) -> Any:

        if str(outfile).endswith(".gz"):
            return gzip.open(temporary_outfile, "wt", encoding="utf-8")
        return open(temporary_outfile, "w", encoding="utf-8")

    @classmethod
    def __iter_gene_info_rows(
        cls,
        gene_info_file: Path,
        reduced: bool,
    ) -> Iterator[Tuple[str, str, str, str, str, str, str]]:

        with cls.__open_gene_info_input(gene_info_file) as file:
            next(file, None)
            for line in file:
                fields = line.rstrip("\n").split("\t")
                if reduced:
                    chromosome = fields[5] if len(fields) > 6 else "-"
                    gene_type = fields[6] if len(fields) > 6 else fields[5]
                    yield (
                        fields[0],
                        fields[1],
                        fields[2],
                        fields[3],
                        fields[4],
                        chromosome,
                        gene_type,
                    )
                else:
                    yield (
                        fields[1],
                        fields[10],
                        fields[2],
                        fields[4],
                        fields[5],
                        fields[6],
                        fields[9],
                    )

    @staticmethod
    def __unlink_if_exists(path: Path) -> None:

        try:
            path.unlink()
        except FileNotFoundError:
            pass


def _is_interaction_list(data: Any) -> bool:

    if not (
        (isinstance(data, SequenceInstance) and not isinstance(data, str))
        or isinstance(data, set)
    ):
        return False

    if len(data) == 0:
        return False

    return all(_is_interaction(item) for item in data)


def _is_interaction(item: Any) -> bool:

    return (
        isinstance(item, SequenceInstance)
        and not isinstance(item, str)
        and len(item) == 3
        and isinstance(item[0], str)
        and isinstance(item[1], str)
        and isinstance(item[2], MappingInstance)
    )


def genesyn(
    organism: str = "mouse",
    version: GeneInfoVersion = "bundled",
    show_warnings: bool = False,
) -> GeneSynonyms:
    """
    Create a GeneSynonyms converter.

    Parameters
    ----------
    organism: str (default: "mouse")
        Common name of the organism of interest.
    version: "bundled", "latest", date string, str path or PathLike (default: "bundled")
        NCBI gene_info version to load.
    show_warnings: bool (default: False)
        If `True`, warn when a requested gene identifier has no correspondence.

    Returns
    -------
    GeneSynonyms
        Converter for gene identifier resolution, conversion and
        standardisation.
    """

    global _SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING

    previous = _SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING
    _SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING = True
    try:
        return GeneSynonyms(
            organism=organism,
            version=version,
            show_warnings=show_warnings,
        )
    finally:
        _SUPPRESS_GENE_SYNONYMS_CONSTRUCTOR_WARNING = previous
