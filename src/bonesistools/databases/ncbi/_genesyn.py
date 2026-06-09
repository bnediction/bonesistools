#!/usr/bin/env python

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

if TYPE_CHECKING:
    from ...boolpy.boolean_network._typing import BooleanNetworkLike

import copy
import ctypes
import inspect
import re
import subprocess
import warnings
from collections import namedtuple
from collections.abc import Mapping as MappingInstance
from collections.abc import Sequence as SequenceInstance
from functools import partial, wraps
from pathlib import Path

import networkx as nx
from networkx import Graph
from pandas import DataFrame
from pandas._typing import Axis

from ..._compat import Literal, get_args
from ._typing import (
    InputIdentifierType,
    OutputIdentifierType,
)

InteractionList = Sequence[Tuple[str, str, Dict[str, int]]]

ORGANISMS = Literal["mouse", "human", "escherichia coli"]

FTP_GENE_INFO = {
    "mouse": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    "human": "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "escherichia coli": (
        "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/"
        "Escherichia_coli_str._K-12_substr._MG1655.gene_info.gz"
    ),
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
    "pseudo": 0,
}

_GENE_SYNONYMS_DEPRECATED_ARGS = {
    "gene_type": "input_identifier_type",
    "alias_gene": "output_identifier_type",
}


def _is_interaction(item: Any) -> bool:
    return (
        isinstance(item, SequenceInstance)
        and not isinstance(item, str)
        and len(item) == 3
        and isinstance(item[0], str)
        and isinstance(item[1], str)
        and isinstance(item[2], MappingInstance)
    )


def _is_interaction_list(data: Any) -> bool:
    if not (
        (isinstance(data, SequenceInstance) and not isinstance(data, str))
        or isinstance(data, set)
    ):
        return False

    if len(data) == 0:
        return False

    return all(_is_interaction(item) for item in data)


def support_legacy_gene_synonyms_args(func):
    """
    Decorate GeneSynonyms methods to accept deprecated argument names.

    Raises
    ------
    TypeError
        If both a deprecated argument name and its replacement are provided.
    """

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

    >>> genesyn = GeneSynonyms(organism="mouse")
    >>> genesyn(["Trp53", "Myc"])
    ['Trp53', 'Myc']

    Convert an interaction list while preserving edge attributes:

    >>> interactions = [("Trp53", "Myc", {"sign": 1})]
    >>> genesyn(interactions)
    [('Trp53', 'Myc', {'sign': 1})]

    Parameters
    ----------
    organism: str (default: "mouse")
        Common name of the organism of interest. Supported organisms are
        `"mouse"`, `"human"` and `"escherichia coli"`.
    force_download: bool (default: False)
        If True, download the NCBI gene_info resource before building mappings.
    show_warnings: bool (default: False)
        If True, warn when a requested gene identifier has no correspondence.

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
    TypeError
        If `force_download`, `show_warnings` or converted data has an
        unsupported type.
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
        force_download: bool = False,
        show_warnings: bool = False,
    ) -> None:

        organism = organism.lower().replace("-", " ")

        if organism not in get_args(ORGANISMS):
            raise ValueError(
                f"invalid argument value for 'organism': "
                f"expected one of {get_args(ORGANISMS)} but received {organism!r}"
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

        Raises
        ------
        TypeError
            If `data` has an unsupported type.
        """

        from ...boolpy.boolean_network._typing import is_boolean_network_like

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
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
            If True, keep the original gene identifier when no correspondence
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
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
            If True, keep the original gene identifier when no correspondence
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
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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

        if axis == 0 or axis == "index":
            iterator = iter(df.index)
        elif axis == 1 or axis == "columns":
            iterator = iter(df.columns)
        else:
            raise ValueError(
                f"invalid argument value for 'axis': "
                f"expected 0, 1, 'index' or 'columns' but received {axis!r}"
            )

        for gene in iterator:
            output_alias = alias_conversion(
                gene=gene, input_identifier_type=input_identifier_type
            )
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
        graph: Graph[Any],
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
        if copy is True:
            return nx.relabel_nodes(graph, mapping=aliases_mapping, copy=True)
        else:
            nx.relabel_nodes(graph, mapping=aliases_mapping, copy=False)
            return None

    @support_legacy_gene_synonyms_args
    def convert_bn(
        self,
        bn: "BooleanNetworkLike",
        input_identifier_type: Union[InputIdentifierType, str] = "name",
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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

        from ...boolpy.boolean_network import BooleanNetwork
        from ...boolpy.boolean_network._typing import is_boolean_network_like

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
            If True, keep the original gene identifier when no correspondence
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
            If True, keep the original gene identifier when no correspondence
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
        force_download: bool = False,
        show_warnings: bool = False,
    ) -> None:
        """
        Reinitialize the converter with a potentially different organism.

        Parameters
        ----------
        organism: str, optional
            Common name of the organism to load. If None, reload the current
            organism.
        force_download: bool (default: False)
            If True, download the NCBI gene_info resource before rebuilding
            mappings.
        show_warnings: bool (default: False)
            If True, warn when a requested gene identifier has no
            correspondence.

        Raises
        ------
        TypeError
            If `force_download` or `show_warnings` is not Boolean.
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
        identifier_type: Union[InputIdentifierType, str] = "name",
    ) -> List[bool]:
        """
        Test whether gene identifiers are present in the synonym resource.

        Examples
        --------
        >>> genesyn = GeneSynonyms(organism="mouse")
        >>> genesyn.contains("Trp53")
        [True]
        >>> genesyn.contains("Trp53", "not-a-gene")
        [True, False]
        >>> genesyn.contains("22059", identifier_type="gene_id")
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
        identifier_type: Union[InputIdentifierType, str] = "name",
    ) -> List[str]:
        """
        Return gene identifiers found in the synonym resource.

        Examples
        --------
        >>> genesyn = GeneSynonyms(organism="mouse")
        >>> genesyn.find("Trp53", "not-a-gene", "Myc")
        ['Trp53', 'Myc']
        >>> genesyn.find("22059", "bad-id", identifier_type="gene_id")
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
        self, gene: str, input_identifier_type: Union[InputIdentifierType, str] = "name"
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
        self, gene: str, input_identifier_type: Union[InputIdentifierType, str] = "name"
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

        gi_file_cut = Path(f"{gi_file}_cut")

        # Output columns:
        # gene_id | official_name | ncbi_name | synonyms | dbXrefs | gene_type
        awk_script = 'NR>1 {print $2 "\t" $11 "\t" $3 "\t" $5 "\t" $6 "\t" $10}'
        cmd = f"awk -F'\t' '{awk_script}' {gi_file} > {gi_file_cut}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("failed to parse NCBI gene_info file") from e

        Identifiers = namedtuple(
            "Identifiers",
            ["official_name", "ncbi_name", "ensembl_id", "databases", "gene_type"],
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

                    gene_info_rows.append(
                        (
                            gene_id,
                            official_name,
                            ncbi_name,
                            fields[3],
                            fields[4],
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

    @support_legacy_gene_synonyms_args
    def __conversion_function(
        self,
        output_identifier_type: Union[OutputIdentifierType, str] = "official_name",
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
