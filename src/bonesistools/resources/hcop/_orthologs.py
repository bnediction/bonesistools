#!/usr/bin/env python

from __future__ import annotations

import copy as copylib
from collections.abc import Mapping as MappingInstance
from collections.abc import Sequence as SequenceInstance
from itertools import product
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
    cast,
    overload,
)

import networkx as nx
import pandas as pd
from networkx import Graph
from typing_extensions import Protocol

from ..._compat import Literal
from ..._validation import _as_positive_integer
from ..._warnings import _deprecated, _warn_deprecated

if TYPE_CHECKING:
    from ...logic.boolean_algebra import Hypercube
    from ...logic.boolean_network._typing import BooleanNetworkLike

InteractionList = Sequence[Tuple[str, str, Dict[str, int]]]
HcopVersion = Union[Literal["bundled", "latest"], str, Path]
HCOP_DIR = Path(__file__).resolve().parent / "data"


class Orthologs:
    """
    Converter between genes from HCOP-supported organisms.

    HCOP tables use human genes as their reference. Human-to-organism and
    organism-to-human conversions therefore use one table, while conversions
    between two non-human organisms join their HCOP tables through shared human
    orthologs. The converter exposes helpers for translating common biological
    objects: sequences, interaction lists, DataFrames, hypercubes, NetworkX
    graphs and Boolean networks.

    Single-gene translation may return several orthologs. Object-level
    translations use the first ortholog in the deterministic ranking built from
    the HCOP table: strongest evidence, then most frequent target, then
    alphabetical order.

    Parameters
    ----------
    input_organism: str (default: "human")
        Organism associated with the input gene symbols. Supported values are
        listed by `bt.resources.hcop.organisms()`.
    output_organism: str (default: "mouse")
        Target organism. Supported values are listed by
        `bt.resources.hcop.organisms()`.
    min_evidence: int (default: 3)
        Minimum number of HCOP support sources required for each retained
        mapping.
    version: "bundled", "latest", str path or Path (default: "bundled")
        HCOP version to load. `"bundled"` loads the HCOP snapshot distributed
        with bonesistools. `"latest"` loads the current public HCOP file.
        File paths resolve directly to user-provided HCOP tables.
    table: pandas.DataFrame, optional
        Preloaded direct mapping with columns `source_symbol`, `target_symbol`,
        `support` and `evidence`. For backward compatibility, `human_symbol`
        can replace `source_symbol` when `input_organism="human"`. If `None`,
        the required HCOP table or tables are loaded automatically.
    target_organism: str, optional
        Alias for `output_organism`.

    Attributes
    ----------
    input_organism: str
        Source organism used by the converter.
    output_organism: str
        Target organism used by the converter.
    min_evidence: int
        Minimum HCOP support threshold used to build mappings.
    version: str
        Normalized HCOP version label used by the converter.

    Raises
    ------
    ValueError
        If organisms are unsupported, `min_evidence` is not positive, or the
        HCOP table is invalid.
    """

    download_url = "https://storage.googleapis.com/public-download-files/hcop"
    supported_organisms = [
        "human",
        "anole_lizard",
        "c.elegans",
        "cat",
        "cattle",
        "chicken",
        "chimpanzee",
        "dog",
        "fruitfly",
        "horse",
        "macaque",
        "mouse",
        "opossum",
        "pig",
        "platypus",
        "rat",
        "s.cerevisiae",
        "s.pombe",
        "xenopus",
        "zebrafish",
    ]
    input_organisms = supported_organisms
    _table: pd.DataFrame
    _mapping: Dict[str, List[str]]
    _best_mapping: Dict[str, str]
    _one_to_many_mapping_cache: Optional[Dict[object, List[object]]]

    def __init__(
        self,
        input_organism: str = "human",
        output_organism: str = "mouse",
        min_evidence: int = 3,
        version: HcopVersion = "bundled",
        table: Optional[pd.DataFrame] = None,
        target_organism: Optional[str] = None,
    ) -> None:
        if target_organism is not None:
            output_organism = target_organism
        self.reset(
            input_organism=input_organism,
            output_organism=output_organism,
            min_evidence=min_evidence,
            version=version,
            table=table,
        )

    @overload
    def __call__(
        self,
        data: str,
        *args: Any,
        **kwargs: Any,
    ) -> List[str]: ...

    @overload
    def __call__(
        self,
        data: List[str],
        *args: Any,
        **kwargs: Any,
    ) -> List[str]: ...

    @overload
    def __call__(
        self,
        data: Tuple[str, ...],
        *args: Any,
        **kwargs: Any,
    ) -> Tuple[str, ...]: ...

    @overload
    def __call__(
        self,
        data: Set[str],
        *args: Any,
        **kwargs: Any,
    ) -> Set[str]: ...

    @overload
    def __call__(
        self,
        data: Sequence[str],
        *args: Any,
        **kwargs: Any,
    ) -> Sequence[str]: ...

    @overload
    def __call__(
        self,
        data: InteractionList,
        *args: Any,
        **kwargs: Any,
    ) -> InteractionList: ...

    @overload
    def __call__(
        self,
        data: pd.DataFrame,
        *args: Any,
        **kwargs: Any,
    ) -> pd.DataFrame: ...

    @overload
    def __call__(
        self,
        data: Graph[Any],
        *args: Any,
        **kwargs: Any,
    ) -> Graph[Any]: ...

    @overload
    def __call__(
        self,
        data: "Hypercube",
        *args: Any,
        **kwargs: Any,
    ) -> "Hypercube": ...

    @overload
    def __call__(
        self,
        data: "BooleanNetworkLike",
        *args: Any,
        **kwargs: Any,
    ) -> "BooleanNetworkLike": ...

    def __call__(
        self,
        data: Union[
            str,
            Sequence[str],
            Set[str],
            InteractionList,
            pd.DataFrame,
            Graph[Any],
            "Hypercube",
            "BooleanNetworkLike",
        ],
        *args: Any,
        **kwargs: Any,
    ) -> Any:
        """
        Translate genes in a supported object.

        This method dispatches to the appropriate translation method according
        to `data` type. Additional positional and keyword arguments are passed
        to the selected translation method.

        Parameters
        ----------
        data: str, sequence, set, InteractionList, DataFrame, Graph, Hypercube
            or BooleanNetworkLike
            Object containing symbols from `input_organism` to translate.
        *args: Any
            Positional arguments forwarded to the selected translation method.
        **kwargs: Any
            Keyword arguments forwarded to the selected translation method.

        Returns
        -------
        object
            Translated object. The exact type depends on `data` and the selected
            translation method.

        """

        if isinstance(data, str):
            return self.translate(data, *args, **kwargs)
        if self._is_interaction_list(data):
            return self.translate_interaction_list(
                cast(InteractionList, data),
                *args,
                **kwargs,
            )
        if (
            isinstance(data, SequenceInstance) and not isinstance(data, str)
        ) or isinstance(data, set):
            return self.translate_sequence(
                cast(Union[Sequence[str], Set[str]], data),
                *args,
                **kwargs,
            )
        if isinstance(data, pd.DataFrame):
            return self.translate_dataframe(data, *args, **kwargs)
        if isinstance(data, Graph):
            return self.translate_graph(data, *args, **kwargs)

        from ...logic.boolean_algebra import Hypercube
        from ...logic.boolean_network._typing import is_boolean_network_like

        if isinstance(data, Hypercube):
            return self.translate_hypercube(data, *args, **kwargs)
        if is_boolean_network_like(data):
            return self.translate_boolean_network(data, *args, **kwargs)
        raise TypeError(
            f"unsupported argument type for 'data': "
            f"expected str, sequence, set, interaction list, {pd.DataFrame}, "
            f"{Graph}, Hypercube or Boolean network-like object "
            f"but received {type(data)}"
        )

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return a deep copy of the loaded orthology table.

        Human-to-organism mappings retain the historical `human_symbol`
        column. Other conversion routes use `source_symbol`. For conversions
        between two non-human organisms, each row also records the intermediate
        `human_symbol`, both branch-specific evidence values, the conservative
        path evidence, the best evidence for the target and the number of
        distinct human paths supporting it.

        Returns
        -------
        pandas.DataFrame
            Copy of the filtered HCOP table used by the converter.
        """

        table = copylib.deepcopy(self._table)
        if self.input_organism == "human":
            table = cast(
                pd.DataFrame,
                table.rename(columns={"source_symbol": "human_symbol"}),
            )
        return table

    def to_dict(self) -> Dict[str, List[str]]:
        """
        Return orthology mappings as source symbol -> target symbols.

        Returns
        -------
        dict
            Copy of the mapping used for translation. Values are ordered by the
            deterministic HCOP ranking.
        """

        return {source: list(targets) for source, targets in self._mapping.items()}

    def translate(self, gene: str, keep_if_missing: bool = False) -> List[str]:
        """
        Translate one gene symbol to the selected output organism.

        A single source gene can have several target orthologs. For one-to-one
        object translations, the first value of this list is used.

        Parameters
        ----------
        gene: str
            Source gene symbol to translate. Complex names separated by `_` are
            translated subunit by subunit.
        keep_if_missing: bool (default: False)
            If `True`, keep the original gene symbol when no ortholog is found.

        Returns
        -------
        list of str
            Target ortholog symbols. An empty list is returned when no mapping
            is found and `keep_if_missing=False`.
        """

        if "_" in gene:
            translated_subunits: List[List[str]] = []
            for subunit in gene.split("_"):
                translated = self.translate(subunit, keep_if_missing=keep_if_missing)
                if len(translated) == 0:
                    return []
                translated_subunits.append(translated)
            return [
                "_".join(translated_complex)
                for translated_complex in product(*translated_subunits)
            ]

        targets = self._mapping.get(gene)
        if targets is not None:
            return list(targets)
        return [gene] if keep_if_missing else []

    @overload
    def translate_sequence(
        self,
        genes: List[str],
        keep_if_missing: bool = True,
    ) -> List[str]: ...

    @overload
    def translate_sequence(
        self,
        genes: Tuple[str, ...],
        keep_if_missing: bool = True,
    ) -> Tuple[str, ...]: ...

    @overload
    def translate_sequence(
        self,
        genes: Set[str],
        keep_if_missing: bool = True,
    ) -> Set[str]: ...

    @overload
    def translate_sequence(
        self,
        genes: Sequence[str],
        keep_if_missing: bool = True,
    ) -> Sequence[str]: ...

    def translate_sequence(
        self,
        genes: Union[Sequence[str], Set[str]],
        keep_if_missing: bool = True,
    ) -> Union[Sequence[str], Set[str]]:
        """
        Translate a sequence or set of source gene symbols.

        Each gene is translated to one target symbol using the deterministic
        HCOP ranking. Missing genes are either kept or removed depending on
        `keep_if_missing`.

        Parameters
        ----------
        genes: sequence or set of str
            Source gene symbols to translate.
        keep_if_missing: bool (default: True)
            If `True`, keep original gene symbols with no ortholog.

        Returns
        -------
        sequence or set of str
            Translated values, preserving the input collection type when
            possible.
        """

        best_mapping = self._best_mapping
        if keep_if_missing:
            translated_genes = [
                cast(str, self._translate_best(gene, keep_if_missing=True))
                if "_" in gene
                else best_mapping.get(gene, gene)
                for gene in genes
            ]
        else:
            translated_genes = []
            for gene in genes:
                translated = (
                    self._translate_best(gene, keep_if_missing=False)
                    if "_" in gene
                    else best_mapping.get(gene)
                )
                if translated is not None:
                    translated_genes.append(translated)

        if isinstance(genes, list):
            return translated_genes

        sequence_constructor = cast(Any, type(genes))
        try:
            return sequence_constructor(translated_genes)
        except TypeError:
            return translated_genes

    def translate_hypercube(
        self,
        hypercube: "Hypercube",
        keep_if_missing: bool = True,
        copy: bool = True,
    ) -> Union["Hypercube", None]:
        """
        Translate hypercube component names.

        Components are translated to one target symbol using the deterministic
        HCOP ranking. Missing components are kept by default because removing
        an explicitly specified component would broaden the represented state
        space.

        Parameters
        ----------
        hypercube: Hypercube
            Hypercube whose explicitly specified components are translated.
        keep_if_missing: bool (default: True)
            If `True`, keep component names with no ortholog. If `False`, raise
            an error when a component has no ortholog.
        copy: bool (default: True)
            Return a translated copy instead of modifying `hypercube`.

        Returns
        -------
        Hypercube or None
            Translated hypercube if `copy=True`; otherwise None.

        Raises
        ------
        TypeError
            If `hypercube` is not a Hypercube.
        ValueError
            If `keep_if_missing=False` and a component has no ortholog, or if
            translation would merge explicitly specified components.
        """

        from ...logic.boolean_algebra import Hypercube

        if not isinstance(hypercube, Hypercube):
            raise TypeError(
                f"unsupported argument type for 'hypercube': "
                f"expected {Hypercube} but received {type(hypercube)}"
            )

        mapping = {}

        for component in hypercube:
            translated = self._translate_best(
                component,
                keep_if_missing=keep_if_missing,
            )

            if translated is None:
                raise ValueError(
                    "hypercube translation cannot remove an unmapped "
                    f"component: {component!r}"
                )

            mapping[component] = translated

        translated_hypercube = hypercube.copy() if copy else hypercube
        translated_hypercube.relabel(mapping)

        return translated_hypercube if copy else None

    def translate_interaction_list(
        self,
        interaction_list: InteractionList,
        keep_if_missing: bool = True,
    ) -> InteractionList:
        """
        Translate source and target genes in an interaction list.

        Edge attributes are preserved. Each source and target is translated to
        one target symbol using the deterministic HCOP ranking.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            Sequence of `(source, target, attributes)` interactions.
        keep_if_missing: bool (default: True)
            If `True`, keep original source or target symbols with no ortholog.
            If `False`, interactions with unmapped endpoints are removed.

        Returns
        -------
        InteractionList
            Translated interaction list with preserved edge attributes.
        """

        translate = self._translate_best
        translated_interactions: List[Tuple[str, str, Dict[str, int]]] = []
        for source, target, attributes in interaction_list:
            translated_source = translate(
                source,
                keep_if_missing,
            )
            translated_target = translate(
                target,
                keep_if_missing,
            )
            if translated_source is None or translated_target is None:
                continue
            translated_interactions.append(
                (translated_source, translated_target, attributes)
            )

        return translated_interactions

    @overload
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: Literal[True] = True,
        one_to_many: Optional[int] = None,
    ) -> pd.DataFrame: ...

    @overload
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        *,
        copy: Literal[False],
        one_to_many: Optional[int] = None,
    ) -> None: ...

    @overload
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]: ...

    def translate_dataframe(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]:
        """
        Translate gene symbols in selected DataFrame columns.

        Each selected value is translated to one target symbol using the
        deterministic HCOP ranking by default. If `one_to_many` is provided,
        use decoupler-like expansion: genes with up to `one_to_many` orthologs
        expand to all orthologs, while genes with more orthologs or no mapping
        are removed.

        Parameters
        ----------
        df: pandas.DataFrame
            DataFrame containing symbols from `input_organism`.
        columns: str or sequence of str, optional
            Columns to translate. If `None`, translate any existing column among
            `source`, `target` and `genesymbol`.
        keep_if_missing: bool (default: True)
            If `True`, keep original symbols with no ortholog. If `False`, rows
            with missing translations in selected columns are removed.
        copy: bool (default: True)
            Return a translated copy instead of modifying `df`.
        one_to_many: int, optional
            Maximum number of orthologs allowed per gene. If provided, use the
            same expansion semantics as `decoupler.op.translate`.

        Returns
        -------
        pandas.DataFrame or None
            Translated DataFrame if `copy=True`; otherwise None.

        Raises
        ------
        ValueError
            If none of the requested columns are present in `df`.
        """

        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"unsupported argument type for 'df': "
                f"expected {pd.DataFrame} but received {type(df)}"
            )

        if one_to_many is not None:
            one_to_many = self._validate_one_to_many(one_to_many)
            if not copy:
                raise TypeError(
                    "invalid argument combination: 'one_to_many' can change "
                    "the number of rows and requires copy=True"
                )

        selected_columns = self._select_columns(df, columns)
        translated = df.copy() if copy else df
        if one_to_many is not None:
            for column in selected_columns:
                translated = self._translate_df_one_to_many(
                    translated,
                    column=column,
                    one_to_many=one_to_many,
                )
            translated = translated.reset_index(drop=True)
            return translated

        for column in selected_columns:
            translated[column] = self._translate_series(
                cast(pd.Series, translated[column]),
                keep_if_missing=keep_if_missing,
            )
            if not keep_if_missing:
                translated = cast(pd.DataFrame, translated.dropna(subset=[column]))

        translated = translated.reset_index(drop=True)
        if copy:
            return translated

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: Literal[True] = True,
        one_to_many: Optional[int] = None,
    ) -> pd.DataFrame: ...

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        *,
        copy: Literal[False],
        one_to_many: Optional[int] = None,
    ) -> None: ...

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]: ...

    @_deprecated(replacement="`Orthologs.translate_dataframe()`")
    def translate_df(
        self,
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]] = None,
        keep_if_missing: bool = True,
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]:
        """Deprecated alias for `translate_dataframe()`."""

        return self.translate_dataframe(
            df,
            columns=columns,
            keep_if_missing=keep_if_missing,
            copy=copy,
            one_to_many=one_to_many,
        )

    def translate_graph(
        self,
        graph: Graph[Any],
        keep_if_missing: bool = True,
        copy: bool = True,
    ) -> Union[Graph[Any], None]:
        """
        Translate graph node labels.

        Nodes are translated to one target symbol using the deterministic HCOP
        ranking. Edge and node attributes are preserved by NetworkX relabeling.

        Parameters
        ----------
        graph: networkx.Graph
            Graph whose node identifiers are symbols from `input_organism`.
        keep_if_missing: bool (default: True)
            If `True`, keep original node labels with no ortholog. If `False`,
            unmapped nodes and incident edges are removed.
        copy: bool (default: True)
            Return a translated copy instead of modifying `graph`.

        Returns
        -------
        networkx.Graph or None
            Translated graph if `copy=True`; otherwise None.

        """

        if not isinstance(graph, Graph):
            raise TypeError(
                f"unsupported argument type for 'graph': "
                f"expected {Graph} but received {type(graph)}"
            )

        translated = graph.copy() if copy else graph
        mapping = {}
        missing = []
        for node in translated.nodes:
            target = self._translate_best(
                node,
                keep_if_missing=keep_if_missing,
            )
            if target is None:
                missing.append(node)
            else:
                mapping[node] = target

        translated.remove_nodes_from(missing)
        nx.relabel_nodes(translated, mapping=mapping, copy=False)

        if copy:
            return translated

    def translate_boolean_network(
        self,
        bn: "BooleanNetworkLike",
        keep_if_missing: bool = True,
        copy: bool = False,
    ) -> Union["BooleanNetworkLike", None]:
        """
        Translate Boolean network component names.

        Components are translated to one target symbol using the deterministic
        HCOP ranking. Component references in Boolean rules are updated
        accordingly.

        Parameters
        ----------
        bn: BooleanNetworkLike
            Boolean network-like object whose component identifiers are symbols
            from `input_organism`.
        keep_if_missing: bool (default: True)
            If `True`, keep original component names with no ortholog. If `False`,
            raise an error for unmapped components.
        copy: bool (default: False)
            Return a translated copy instead of modifying `bn`.

        Returns
        -------
        BooleanNetworkLike or None
            Translated Boolean network if `copy=True`; otherwise None.

        Raises
        ------
        ValueError
            If `keep_if_missing=False` and a component has no ortholog.
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
        relabel = getattr(bn_any, "relabel", None)

        if copy:
            copy_method = getattr(bn_any, "copy", None)
            bn_any = (
                copy_method()
                if callable(copy_method)
                else BooleanNetwork(bn_any, check=False)
            )
            rename = getattr(bn_any, "rename", None)
            relabel = getattr(bn_any, "relabel", None)
            if not callable(rename) and not callable(relabel):
                bn_any = BooleanNetwork(bn_any, check=False)
                rename = bn_any.rename
                relabel = bn_any.relabel
                rebuild_as_input_type = True
        elif not callable(rename) and not callable(relabel):
            raise TypeError(
                "unsupported argument value for 'bn': in-place Boolean network "
                "translation requires a 'relabel' or 'rename' method; use "
                "copy=True to return a translated copy"
            )

        genes = tuple(gene for gene, _ in bn_any.items())
        mapping = {}
        for gene in genes:
            translated = self._translate_best(gene, keep_if_missing=keep_if_missing)
            if translated is None:
                raise ValueError(
                    "Boolean network translation cannot remove an unmapped "
                    f"component: {gene!r}"
                )
            mapping[gene] = translated

        if callable(relabel):
            relabel(mapping)
        elif callable(rename):
            for gene, translated in mapping.items():
                rename(gene, translated)

        if copy and rebuild_as_input_type:
            try:
                return cast(
                    "BooleanNetworkLike",
                    cast(Any, input_type)(bn_any.rules),
                )
            except Exception as error:
                raise TypeError(
                    "unable to return translated Boolean network with original "
                    f"type {input_type}"
                ) from error

        return cast("BooleanNetworkLike", bn_any) if copy else None

    @_deprecated(replacement="`Orthologs.translate_boolean_network()`")
    def translate_bn(
        self,
        bn: "BooleanNetworkLike",
        keep_if_missing: bool = True,
        copy: bool = False,
    ) -> Union["BooleanNetworkLike", None]:
        """Deprecated alias for `translate_boolean_network()`."""

        return self.translate_boolean_network(
            bn,
            keep_if_missing=keep_if_missing,
            copy=copy,
        )

    def reset(
        self,
        input_organism: Optional[str] = None,
        output_organism: Optional[str] = None,
        min_evidence: Optional[int] = None,
        version: Optional[HcopVersion] = None,
        table: Optional[pd.DataFrame] = None,
        target_organism: Optional[str] = None,
    ) -> None:
        """
        Reinitialize the converter with a potentially different organism.

        Parameters
        ----------
        input_organism: str, optional
            Source organism. If `None`, keep the current source organism.
        output_organism: str, optional
            Target organism. If `None`, keep the current target organism.
        min_evidence: int, optional
            Minimum HCOP support threshold. If `None`, keep the current threshold.
        version: "bundled", "latest", str path or Path, optional
            HCOP version to load. If `None`, keep the current version.
        table: pandas.DataFrame, optional
            Preloaded direct source-to-target table. If `None`, reload the
            required mapping from HCOP.
        target_organism: str, optional
            Alias for `output_organism`.

        Raises
        ------
        ValueError
            If organisms are unsupported, `min_evidence` is not positive, or
            the HCOP table is invalid.
        """

        if target_organism is not None:
            output_organism = target_organism
        if input_organism is None:
            resolved_input_organism = cast(
                str,
                getattr(self, "input_organism", "human"),
            )
        else:
            resolved_input_organism = input_organism
        if output_organism is None:
            resolved_output_organism = cast(
                str,
                getattr(self, "output_organism", "mouse"),
            )
        else:
            resolved_output_organism = output_organism
        if min_evidence is None:
            resolved_min_evidence = cast(
                int,
                getattr(self, "min_evidence", 3),
            )
        else:
            resolved_min_evidence = min_evidence
        if version is None:
            resolved_version = cast(
                HcopVersion,
                getattr(self, "version", "bundled"),
            )
        else:
            resolved_version = version
        normalized_version = self._normalize_version(resolved_version)

        self._validate_input_organism(resolved_input_organism)
        self._validate_output_organism(resolved_output_organism)
        if resolved_input_organism == resolved_output_organism:
            raise ValueError(
                "invalid organism pair: 'input_organism' and "
                "'output_organism' must differ"
            )
        resolved_min_evidence = self._validate_min_evidence(resolved_min_evidence)

        if table is None:
            table = self._load_orthologs_table(
                input_organism=resolved_input_organism,
                output_organism=resolved_output_organism,
                min_evidence=resolved_min_evidence,
                version=normalized_version,
            )
        else:
            table = copylib.deepcopy(table)
            table = self._normalize_direct_table(
                table,
                input_organism=resolved_input_organism,
            )

        self.input_organism = resolved_input_organism
        self.output_organism = resolved_output_organism
        self.target_organism = resolved_output_organism
        self.min_evidence = resolved_min_evidence
        self.version = normalized_version
        self._table = cast(pd.DataFrame, table)
        self._mapping = self._build_mapping(self._table)
        self._best_mapping = {
            source: targets[0] for source, targets in self._mapping.items()
        }
        self._one_to_many_mapping_cache = None

    def contains(self, *genes: str) -> List[bool]:
        """
        Test whether genes are present in the loaded HCOP resource.

        Parameters
        ----------
        *genes: str
            Human gene symbols to test.

        Returns
        -------
        list of bool
            Boolean values ordered like `genes`. An empty input returns an empty
            list.
        """

        mapping = self._mapping
        return [gene in mapping for gene in genes]

    def find(self, *genes: str) -> List[str]:
        """
        Return genes found in the loaded HCOP resource.

        Parameters
        ----------
        *genes: str
            Human gene symbols to test.

        Returns
        -------
        list of str
            Input gene symbols present in the HCOP resource, preserving their
            input order.
        """

        mapping = self._mapping
        return [gene for gene in genes if gene in mapping]

    def _translate_best(self, gene: str, keep_if_missing: bool) -> Optional[str]:

        if "_" in gene:
            translated_subunits: List[str] = []
            for subunit in gene.split("_"):
                translated = self._best_mapping.get(subunit)
                if translated is None:
                    if not keep_if_missing:
                        return None
                    translated = subunit
                translated_subunits.append(translated)
            return "_".join(translated_subunits)

        translated = self._best_mapping.get(gene)
        if translated is not None:
            return translated
        return gene if keep_if_missing else None

    def _translate_series(
        self,
        values: pd.Series,
        keep_if_missing: bool,
    ) -> pd.Series:
        """Translate a Series while preserving complex-name semantics."""

        source = values.astype(str)
        best_mapping = self._best_mapping
        translate = self._translate_best
        if keep_if_missing:
            translated = [
                translate(gene, keep_if_missing=True)
                if "_" in gene
                else best_mapping.get(gene, gene)
                for gene in source
            ]
        else:
            translated = [
                translate(gene, keep_if_missing=False)
                if "_" in gene
                else best_mapping.get(gene)
                for gene in source
            ]
        return pd.Series(translated, index=values.index, name=values.name)

    def _translate_df_one_to_many(
        self,
        df: pd.DataFrame,
        column: str,
        one_to_many: int,
    ) -> pd.DataFrame:

        expansions = self._ortholog_expansions(
            values=cast(pd.Series, df[column]),
            one_to_many=one_to_many,
        )
        translated = df
        translated[column] = cast(pd.Series, translated[column]).map(expansions)
        translated = cast(pd.DataFrame, translated.explode(column))
        return cast(pd.DataFrame, translated.dropna(subset=[column]))

    @staticmethod
    def _load_orthologs_table(
        input_organism: str = "human",
        output_organism: str = "mouse",
        min_evidence: int = 3,
        version: HcopVersion = "bundled",
        target_organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Load HCOP orthology mappings as a filtered DataFrame.

        Human-to-organism mappings are loaded directly, organism-to-human
        mappings are reversed, and mappings between two non-human organisms are
        joined through shared human orthologs.

        Parameters
        ----------
        input_organism: str (default: "human")
            Source organism.
        output_organism: str (default: "mouse")
            Target organism. Supported values are listed by
            `bt.resources.hcop.organisms()`.
        min_evidence: int (default: 3)
            Minimum number of HCOP support sources required for a mapping.
            The support count is computed from the comma-separated `support`
            column.
        version: "bundled", "latest", str path or Path (default: "bundled")
            HCOP version to load.

            If `"bundled"`, load the HCOP snapshot distributed with
            bonesistools.

            If `"latest"`, load the current public HCOP file.

            If a path is provided, load that file directly. A single custom
            HCOP file can only describe a route with human at one endpoint.
        target_organism: str, optional
            Alias for `output_organism`.

        Returns
        -------
        pandas.DataFrame
            Filtered source-to-target mapping. The internal source column is
            named `source_symbol`.

        Raises
        ------
        ValueError
            If an organism is unsupported, `min_evidence` is not positive,
            `version` is invalid, or the downloaded HCOP table is missing
            required columns.
        FileNotFoundError
            If the bundled or file-based HCOP version cannot be resolved.
        """

        if target_organism is not None:
            output_organism = target_organism
        Orthologs._validate_input_organism(input_organism)
        Orthologs._validate_output_organism(output_organism)
        if input_organism == output_organism:
            raise ValueError(
                "invalid organism pair: 'input_organism' and "
                "'output_organism' must differ"
            )
        min_evidence = Orthologs._validate_min_evidence(min_evidence)
        version_label = Orthologs._normalize_version(version)

        if (
            input_organism != "human"
            and output_organism != "human"
            and version_label not in {"bundled", "latest"}
        ):
            raise ValueError(
                "non-human to non-human conversion requires two HCOP tables; "
                "use version='bundled', version='latest', or provide a direct "
                "table with source_symbol and target_symbol columns"
            )

        if input_organism == "human":
            direct = Orthologs._load_human_orthologs_table(
                output_organism,
                min_evidence=min_evidence,
                version=version_label,
            )
            return cast(
                pd.DataFrame,
                direct.rename(columns={"human_symbol": "source_symbol"}),
            )

        if output_organism == "human":
            direct = Orthologs._load_human_orthologs_table(
                input_organism,
                min_evidence=min_evidence,
                version=version_label,
            )
            valid = Orthologs._without_placeholder_symbols(direct)
            return pd.DataFrame(
                {
                    "source_symbol": valid["target_symbol"],
                    "target_symbol": valid["human_symbol"],
                    "support": valid["support"],
                    "evidence": valid["evidence"],
                }
            ).reset_index(drop=True)

        source = Orthologs._load_human_orthologs_table(
            input_organism,
            min_evidence=min_evidence,
            version=version_label,
        )
        target = Orthologs._load_human_orthologs_table(
            output_organism,
            min_evidence=min_evidence,
            version=version_label,
        )
        return Orthologs._bridge_hcop_tables(source, target)

    @staticmethod
    def _load_human_orthologs_table(
        output_organism: str,
        min_evidence: int,
        version: HcopVersion,
    ) -> pd.DataFrame:
        """Load one filtered human-to-organism HCOP table."""

        hcop = Orthologs._read_hcop_table(output_organism, version=version)
        target_column = Orthologs._target_symbol_column(output_organism)
        if target_column not in hcop and "target_symbol" in hcop:
            target_column = "target_symbol"
        required_columns = ["human_symbol", target_column, "support"]
        missing_columns = [column for column in required_columns if column not in hcop]
        if missing_columns:
            raise ValueError(
                "invalid HCOP table: missing required columns "
                f"{', '.join(missing_columns)}"
            )

        table = cast(pd.DataFrame, hcop[required_columns].dropna()).copy()
        support = cast(pd.Series, table["support"]).astype(str)
        table["evidence"] = support.str.count(",") + 1
        table = cast(pd.DataFrame, table[table["evidence"] >= min_evidence])
        table = cast(
            pd.DataFrame,
            table.rename(columns={target_column: "target_symbol"}).reset_index(
                drop=True
            ),
        )
        return table

    @staticmethod
    def _without_placeholder_symbols(table: pd.DataFrame) -> pd.DataFrame:
        """Remove HCOP rows whose source or target symbol is a placeholder."""

        valid = (table["human_symbol"] != "-") & (table["target_symbol"] != "-")
        return cast(pd.DataFrame, table.loc[valid]).copy()

    @staticmethod
    def _bridge_hcop_tables(
        source: pd.DataFrame,
        target: pd.DataFrame,
    ) -> pd.DataFrame:
        """Project two HCOP mappings through their shared human symbols."""

        source_edges = Orthologs._without_placeholder_symbols(source)
        source_edges = cast(
            pd.DataFrame,
            source_edges.sort_values(
                "evidence",
                ascending=False,
                kind="mergesort",
            ).drop_duplicates(["human_symbol", "target_symbol"]),
        )
        source_edges = cast(
            pd.DataFrame,
            source_edges.rename(
                columns={
                    "target_symbol": "source_symbol",
                    "support": "source_support",
                    "evidence": "source_evidence",
                }
            ),
        )

        target_edges = Orthologs._without_placeholder_symbols(target)
        target_edges = cast(
            pd.DataFrame,
            target_edges.sort_values(
                "evidence",
                ascending=False,
                kind="mergesort",
            ).drop_duplicates(["human_symbol", "target_symbol"]),
        )
        target_edges = cast(
            pd.DataFrame,
            target_edges.rename(
                columns={
                    "support": "target_support",
                    "evidence": "target_evidence",
                }
            ),
        )

        paths = cast(
            pd.DataFrame,
            source_edges.merge(target_edges, on="human_symbol", how="inner"),
        )
        paths["evidence"] = paths[["source_evidence", "target_evidence"]].min(axis=1)
        paths["support"] = (
            paths["source_support"].astype(str)
            + " -> "
            + paths["target_support"].astype(str)
        )

        grouped_paths = paths.groupby(
            ["source_symbol", "target_symbol"],
            sort=False,
        )
        paths["best_evidence"] = grouped_paths["evidence"].transform("max")
        paths["paths"] = grouped_paths["human_symbol"].transform("nunique")
        paths = cast(
            pd.DataFrame,
            paths.sort_values(
                ["source_symbol", "target_symbol", "human_symbol"],
                kind="mergesort",
            ).reset_index(drop=True),
        )
        return paths[
            [
                "source_symbol",
                "target_symbol",
                "human_symbol",
                "source_support",
                "target_support",
                "source_evidence",
                "target_evidence",
                "evidence",
                "best_evidence",
                "paths",
                "support",
            ]
        ]

    @staticmethod
    def _is_interaction(item: Any) -> bool:

        return (
            isinstance(item, SequenceInstance)
            and not isinstance(item, str)
            and len(item) == 3
            and isinstance(item[0], str)
            and isinstance(item[1], str)
            and isinstance(item[2], MappingInstance)
        )

    @staticmethod
    def _is_interaction_list(data: Any) -> bool:

        if not (
            (isinstance(data, SequenceInstance) and not isinstance(data, str))
            or isinstance(data, set)
        ):
            return False
        if len(data) == 0:
            return False
        return Orthologs._is_interaction(next(iter(data)))

    @staticmethod
    def _select_columns(
        df: pd.DataFrame,
        columns: Optional[Union[str, Sequence[str]]],
    ) -> List[str]:

        if columns is None:
            requested_columns = ["source", "target", "genesymbol"]
        elif isinstance(columns, str):
            requested_columns = [columns]
        elif isinstance(columns, Sequence):
            requested_columns = list(columns)
        else:
            raise TypeError(
                f"unsupported argument type for 'columns': "
                f"expected {str} or sequence of {str} but received {type(columns)}"
            )

        if not all(isinstance(column, str) for column in requested_columns):
            raise TypeError("all values in 'columns' must be strings")

        selected_columns = [column for column in requested_columns if column in df]
        if not selected_columns:
            raise ValueError(
                "invalid argument value for 'columns': none of the requested columns "
                f"are present in df ({list(df.columns)})"
            )

        return selected_columns

    def _ortholog_expansions(
        self,
        values: pd.Series,
        one_to_many: int,
    ) -> Dict[str, List[str]]:
        """Return row-ordered target expansions for distinct source values."""

        map_dict = self._one_to_many_mapping()
        expansions: Dict[str, List[str]] = {}
        for value in values.drop_duplicates():
            source = str(value)
            subunit_lists: List[List[str]] = []
            for subunit in source.split("_"):
                targets = map_dict.get(subunit, [])
                if len(targets) == 0 or len(targets) > one_to_many:
                    subunit_lists.append([])
                else:
                    subunit_lists.append([str(target) for target in targets])
            if any(len(subunit_list) == 0 for subunit_list in subunit_lists):
                continue
            expansions[source] = [
                "_".join(translated_complex)
                for translated_complex in product(*subunit_lists)
            ]

        return expansions

    def _one_to_many_mapping(self) -> Dict[object, List[object]]:
        """Return the cached row-ordered mapping used for expansion."""

        mapping = self._one_to_many_mapping_cache
        if mapping is None:
            if self.input_organism == "human":
                row_mapping: Dict[object, List[object]] = {}
                sources = self._table["source_symbol"]
                targets = self._table["target_symbol"]
                for source, target in zip(sources, targets):
                    if pd.isna(source):
                        continue
                    row_mapping.setdefault(source, []).append(target)
                mapping = row_mapping
            else:
                mapping = cast(
                    Dict[object, List[object]],
                    {
                        source: list(targets)
                        for source, targets in self._mapping.items()
                    },
                )
            self._one_to_many_mapping_cache = mapping
        return mapping

    @staticmethod
    def _target_symbol_column(output_organism: str) -> str:

        if output_organism == "anole_lizard":
            return "anole lizard_symbol"
        if output_organism == "fruitfly":
            return "fruit fly_symbol"
        return f"{output_organism}_symbol"

    @staticmethod
    def _validate_input_organism(input_organism: str) -> None:

        if not isinstance(input_organism, str):
            raise TypeError(
                f"unsupported argument type for 'input_organism': "
                f"expected {str} but received {type(input_organism)}"
            )
        if input_organism not in Orthologs.input_organisms:
            raise ValueError(
                "invalid argument value for 'input_organism': expected one of "
                f"{Orthologs.input_organisms}, but received {input_organism!r}"
            )

    @staticmethod
    def _validate_output_organism(output_organism: str) -> None:

        if not isinstance(output_organism, str):
            raise TypeError(
                f"unsupported argument type for 'output_organism': "
                f"expected {str} but received {type(output_organism)}"
            )
        if output_organism not in Orthologs.supported_organisms:
            raise ValueError(
                "invalid argument value for 'output_organism': expected one of "
                f"{Orthologs.supported_organisms}, but received {output_organism!r}"
            )

    @staticmethod
    def _validate_min_evidence(min_evidence: int) -> int:

        return _as_positive_integer(min_evidence, "min_evidence")

    @staticmethod
    def _validate_one_to_many(one_to_many: int) -> int:

        return _as_positive_integer(one_to_many, "one_to_many")

    @staticmethod
    def _normalize_version(version: HcopVersion) -> str:

        if isinstance(version, Path):
            return str(version)

        if not isinstance(version, str):
            raise TypeError(
                f"unsupported argument type for 'version': "
                f"expected {str} or {Path} but received {type(version)}"
            )

        version_label = version.strip()
        if version_label == "":
            return "bundled"
        if version_label.lower() in {"bundled", "latest"}:
            return version_label.lower()

        return version_label

    @staticmethod
    def _normalize_direct_table(
        table: pd.DataFrame,
        input_organism: str,
    ) -> pd.DataFrame:
        """Normalize a user-provided direct source-to-target mapping."""

        normalized = table
        if "source_symbol" not in normalized and input_organism == "human":
            if "human_symbol" in normalized:
                normalized = cast(
                    pd.DataFrame,
                    normalized.rename(columns={"human_symbol": "source_symbol"}),
                )
        Orthologs._validate_orthologs_table(normalized)
        return normalized

    @staticmethod
    def _validate_orthologs_table(table: pd.DataFrame) -> None:

        required_columns = ["source_symbol", "target_symbol", "support", "evidence"]
        missing_columns = [column for column in required_columns if column not in table]
        if missing_columns:
            raise ValueError(
                "invalid HCOP table: missing required columns "
                f"{', '.join(missing_columns)}"
            )

    @staticmethod
    def _build_mapping(table: pd.DataFrame) -> Dict[str, List[str]]:

        scores: Dict[Tuple[str, str], Tuple[int, int]] = {}
        bridge = "best_evidence" in table and "paths" in table
        valid = table["source_symbol"].notna() & table["target_symbol"].notna()
        source_values = table.loc[valid, "source_symbol"].tolist()
        target_values = table.loc[valid, "target_symbol"].tolist()
        evidence_column = "best_evidence" if bridge else "evidence"
        evidence_values = table.loc[valid, evidence_column].tolist()

        if bridge:
            count_values = table.loc[valid, "paths"].tolist()
            for source, target, evidence_value, count_value in zip(
                source_values,
                target_values,
                evidence_values,
                count_values,
            ):
                key = (str(source), str(target))
                evidence = int(evidence_value)
                count = int(count_value)
                previous = scores.get(key)
                scores[key] = (
                    evidence if previous is None else max(previous[0], evidence),
                    count if previous is None else max(previous[1], count),
                )
        else:
            for source, target, evidence_value in zip(
                source_values,
                target_values,
                evidence_values,
            ):
                key = (str(source), str(target))
                evidence = int(evidence_value)
                previous = scores.get(key)
                scores[key] = (
                    evidence if previous is None else max(previous[0], evidence),
                    1 if previous is None else previous[1] + 1,
                )

        grouped: Dict[str, List[Tuple[str, int, int]]] = {}
        for (source, target), (evidence, count) in scores.items():
            grouped.setdefault(source, []).append((target, evidence, count))

        return {
            source: [
                target
                for target, _, _ in sorted(
                    grouped[source],
                    key=lambda value: (-value[1], -value[2], value[0]),
                )
            ]
            for source in sorted(grouped)
        }

    @staticmethod
    def _resolve_hcop_table_source(output_organism: str, version: str) -> str:

        if version == "bundled":
            bundled_files = [
                HCOP_DIR / f"human_{output_organism}_hcop.tsv.gz",
                HCOP_DIR / f"human_{output_organism}_hcop.tsv",
                HCOP_DIR / f"human_{output_organism}_hcop_fifteen_column.txt.gz",
                HCOP_DIR / f"human_{output_organism}_hcop_fifteen_column.txt",
            ]
            for bundled_file in bundled_files:
                if bundled_file.exists():
                    return str(bundled_file)
            raise FileNotFoundError(
                "bundled HCOP file not found for "
                f"{output_organism!r}: "
                f"{', '.join(str(file) for file in bundled_files)}"
            )

        if version == "latest":
            return (
                f"{Orthologs.download_url}/"
                f"human_{output_organism}_hcop_fifteen_column.txt.gz"
            )

        version_path = Path(version)
        if version_path.exists():
            return str(version_path)

        raise FileNotFoundError(
            f"HCOP version file not found: {version}. "
            f"Available versions: {', '.join(versions())}"
        )

    @staticmethod
    def _read_hcop_table(
        output_organism: str, version: HcopVersion = "bundled"
    ) -> pd.DataFrame:

        source = Orthologs._resolve_hcop_table_source(
            output_organism=output_organism,
            version=Orthologs._normalize_version(version),
        )
        return pd.read_csv(
            source,
            sep="\t",
            compression="infer",
            low_memory=False,
        )


def organisms() -> List[str]:
    """
    List organisms supported by the HCOP translation table.
    """

    return list(Orthologs.supported_organisms)


def versions() -> List[str]:
    """
    List version labels accepted by the HCOP resource loader.

    Returns
    -------
    list of str
        `"bundled"` denotes the HCOP snapshot distributed with bonesistools.
        `"latest"` denotes the current public HCOP files.
    """

    return ["bundled", "latest"]


class _OrthologsCallable(Protocol):
    def __call__(
        self,
        output_organism: str = "mouse",
        min_evidence: int = 3,
        version: HcopVersion = "bundled",
        target_organism: Optional[str] = None,
        *,
        input_organism: str = "human",
    ) -> Orthologs: ...


def _orthologs(
    output_organism: str = "mouse",
    min_evidence: int = 3,
    version: HcopVersion = "bundled",
    target_organism: Optional[str] = None,
    *,
    input_organism: str = "human",
) -> Orthologs:
    """
    Create an HCOP orthology converter.

    HCOP uses human genes as its reference. Conversions involving human use one
    HCOP table; conversions between two non-human organisms join their tables
    through shared human orthologs. Single-gene translation can return several
    orthologs, while object-level translations use a deterministic one-to-one
    choice: strongest evidence, then most supporting paths, then alphabetical
    order.

    Parameters
    ----------
    output_organism: str (default: "mouse")
        Target organism. Supported values are listed by
        `bt.resources.hcop.organisms()`.
    min_evidence: int (default: 3)
        Minimum number of HCOP support sources required for each retained
        mapping.
    version: "bundled", "latest", str path or Path (default: "bundled")
        HCOP version to load. `"bundled"` loads the snapshot distributed with
        bonesistools, `"latest"` loads the current public HCOP file, and paths
        resolve to custom HCOP-like files.
    target_organism: str, optional
        Alias for `output_organism`.
    input_organism: str (default: "human")
        Organism associated with the input gene symbols. Supported values are
        listed by `bt.resources.hcop.organisms()`.

    Returns
    -------
    Orthologs
        Converter initialized with mappings from `input_organism` to
        `output_organism`.

    Raises
    ------
    ValueError
        If `output_organism` is unsupported, `min_evidence` is not positive, or
        the HCOP table is invalid.
    """

    if target_organism is not None:
        output_organism = target_organism
    return Orthologs(
        input_organism=input_organism,
        output_organism=output_organism,
        min_evidence=min_evidence,
        version=version,
    )


def _orthologs_versions(output_organism: Optional[str] = None) -> List[str]:
    """Deprecated compatibility alias for `bt.resources.hcop.versions()`."""

    _warn_deprecated(
        "`bt.resources.hcop.orthologs.versions()`",
        replacement="`bt.resources.hcop.versions()`",
        stacklevel=3,
    )
    if output_organism is not None:
        Orthologs._validate_output_organism(output_organism)
    return versions()


setattr(_orthologs, "versions", _orthologs_versions)
orthologs: _OrthologsCallable = cast(_OrthologsCallable, _orthologs)
