#!/usr/bin/env python

from __future__ import annotations

import copy as copylib
from collections.abc import Mapping as MappingInstance
from collections.abc import MutableMapping as MutableMappingInstance
from collections.abc import Sequence as SequenceInstance
from functools import wraps
from itertools import product
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    TypeVar,
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
from ..._warnings import (
    _deprecated,
    _warn_deprecated,
    _warn_deprecated_argument,
)

if TYPE_CHECKING:
    from ...logic.boolean_algebra import Hypercube
    from ...logic.boolean_algebra._typing import HypercubeLike
    from ...logic.boolean_network._typing import BooleanNetworkLike

InteractionList = Sequence[Tuple[str, str, Dict[str, int]]]
HcopVersion = Union[Literal["bundled", "latest"], str, Path]
HCOP_DIR = Path(__file__).resolve().parent / "data"
_F = TypeVar("_F", bound=Callable[..., Any])
_HypercubeLikeT = TypeVar("_HypercubeLikeT", bound="HypercubeLike")


def _support_legacy_keep_if_missing(
    *,
    positional_parameters: Tuple[str, ...],
    false_policy: Literal["raise", "drop"],
) -> Callable[[_F], _F]:
    """Accept the deprecated missing-symbol boolean as an unmapped policy."""

    def decorator(function: _F) -> _F:
        @wraps(function)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            positional = args[:2]
            legacy_values = args[2:]
            if len(legacy_values) > len(positional_parameters):
                raise TypeError(
                    f"expected at most {len(positional_parameters) + 1} "
                    "positional arguments"
                )
            for name, value in zip(positional_parameters, legacy_values):
                if name in kwargs:
                    raise TypeError(f"multiple values for argument {name!r}")
                kwargs[name] = value

            if "keep_if_missing" in kwargs:
                if "on_unmapped" in kwargs:
                    raise TypeError(
                        "cannot specify both 'keep_if_missing' and 'on_unmapped'"
                    )
                keep_if_missing = kwargs.pop("keep_if_missing")
                if not isinstance(keep_if_missing, bool):
                    raise TypeError("'keep_if_missing' must be a boolean")
                _warn_deprecated_argument(
                    "keep_if_missing",
                    "on_unmapped",
                    stacklevel=4,
                )
                kwargs["on_unmapped"] = "keep" if keep_if_missing else false_policy

            return function(*positional, **kwargs)

        return cast(_F, wrapped)

    return decorator


def _validate_on_unmapped(
    on_unmapped: str,
    *,
    allowed: Tuple[str, ...],
) -> str:
    """Validate a missing-symbol policy against the supported values."""

    if on_unmapped not in allowed:
        expected = ", ".join(repr(policy) for policy in allowed)
        raise ValueError(
            "invalid argument value for 'on_unmapped': "
            f"expected one of {expected} but received {on_unmapped!r}"
        )
    return on_unmapped


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
        **kwargs: Any,
    ) -> List[str]: ...

    @overload
    def __call__(
        self,
        data: List[str],
        **kwargs: Any,
    ) -> List[str]: ...

    @overload
    def __call__(
        self,
        data: Tuple[str, ...],
        **kwargs: Any,
    ) -> Tuple[str, ...]: ...

    @overload
    def __call__(
        self,
        data: Set[str],
        **kwargs: Any,
    ) -> Set[str]: ...

    @overload
    def __call__(
        self,
        data: Sequence[str],
        **kwargs: Any,
    ) -> Sequence[str]: ...

    @overload
    def __call__(
        self,
        data: InteractionList,
        **kwargs: Any,
    ) -> InteractionList: ...

    @overload
    def __call__(
        self,
        data: pd.DataFrame,
        **kwargs: Any,
    ) -> pd.DataFrame: ...

    @overload
    def __call__(
        self,
        data: Graph[Any],
        **kwargs: Any,
    ) -> Graph[Any]: ...

    @overload
    def __call__(
        self,
        data: "Hypercube",
        **kwargs: Any,
    ) -> "Hypercube": ...

    @overload
    def __call__(
        self,
        data: "BooleanNetworkLike",
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
        **kwargs: Any,
    ) -> Any:
        """
        Translate genes in a supported object.

        This method dispatches to the appropriate translation method according
        to `data` type. Additional keyword arguments are passed to the selected
        translation method.

        Parameters
        ----------
        data: str, sequence, set, InteractionList, DataFrame, Graph, Hypercube
            or BooleanNetworkLike
            Object containing symbols from `input_organism` to translate.
        **kwargs: Any
            Keyword arguments forwarded to the selected translation method.

        Returns
        -------
        object
            Translated object. The exact type depends on `data` and the selected
            translation method.

        """

        if isinstance(data, str):
            return self.translate(data, **kwargs)
        if self._is_interaction_list(data):
            return self.translate_interaction_list(
                cast(InteractionList, data),
                **kwargs,
            )
        if (
            isinstance(data, SequenceInstance) and not isinstance(data, str)
        ) or isinstance(data, set):
            return self.translate_sequence(
                cast(Union[Sequence[str], Set[str]], data),
                **kwargs,
            )
        if isinstance(data, pd.DataFrame):
            return self.translate_dataframe(data, **kwargs)
        if isinstance(data, Graph):
            return self.translate_graph(data, **kwargs)

        from ...logic.boolean_algebra import Hypercube
        from ...logic.boolean_network._typing import is_boolean_network_like

        if isinstance(data, Hypercube):
            return self.translate_hypercube(data, **kwargs)
        if is_boolean_network_like(data):
            return self.translate_boolean_network(data, **kwargs)
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

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing",),
        false_policy="drop",
    )
    def translate(
        self,
        gene: str,
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> List[str]:
        """
        Translate one gene symbol to the selected output organism.

        A single source gene can have several target orthologs. For one-to-one
        object translations, the first value of this list is used.

        Parameters
        ----------
        gene: str
            Source gene symbol to translate. Complex names separated by `_` are
            translated subunit by subunit.
        on_unmapped: {"raise", "keep", "drop"} (default: "raise")
            Strategy for a gene symbol with no target ortholog. `"raise"`
            raises a `ValueError`, `"keep"` returns the original symbol and
            `"drop"` returns an empty list.

        Returns
        -------
        list of str
            Target ortholog symbols.

        Raises
        ------
        ValueError
            If no target ortholog exists and `on_unmapped="raise"`.
        """

        on_unmapped = cast(
            Literal["raise", "keep", "drop"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep", "drop"),
            ),
        )

        if "_" in gene:
            translated_subunits: List[List[str]] = []
            for subunit in gene.split("_"):
                translated = self.translate(subunit, on_unmapped=on_unmapped)
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
        if on_unmapped == "raise":
            raise ValueError(f"no target ortholog found for gene symbol {gene!r}")
        return [gene] if on_unmapped == "keep" else []

    @overload
    def translate_sequence(
        self,
        genes: List[str],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> List[str]: ...

    @overload
    def translate_sequence(
        self,
        genes: Tuple[str, ...],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> Tuple[str, ...]: ...

    @overload
    def translate_sequence(
        self,
        genes: Set[str],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> Set[str]: ...

    @overload
    def translate_sequence(
        self,
        genes: Sequence[str],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> Sequence[str]: ...

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing",),
        false_policy="drop",
    )
    def translate_sequence(
        self,
        genes: Union[Sequence[str], Set[str]],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> Union[Sequence[str], Set[str]]:
        """
        Translate a sequence or set of source gene symbols.

        Each gene is translated to one target symbol using the deterministic
        HCOP ranking.

        Parameters
        ----------
        genes: sequence or set of str
            Source gene symbols to translate.
        on_unmapped: {"raise", "keep", "drop"} (default: "raise")
            Strategy for gene symbols with no target ortholog. `"raise"`
            raises a `ValueError`, `"keep"` preserves the original symbols and
            `"drop"` removes them from the result.

        Returns
        -------
        sequence or set of str
            Translated values, preserving the input collection type when
            possible.
        """

        on_unmapped = cast(
            Literal["raise", "keep", "drop"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep", "drop"),
            ),
        )
        best_mapping = self._best_mapping
        if on_unmapped == "keep":
            translated_genes = [
                cast(str, self._translate_best(gene, on_unmapped="keep"))
                if "_" in gene
                else best_mapping.get(gene, gene)
                for gene in genes
            ]
        elif on_unmapped == "drop":
            translated_genes = []
            for gene in genes:
                translated = (
                    self._translate_best(gene, on_unmapped="drop")
                    if "_" in gene
                    else best_mapping.get(gene)
                )
                if translated is not None:
                    translated_genes.append(translated)
        else:
            translated_genes = [
                cast(str, self._translate_best(gene, on_unmapped="raise"))
                for gene in genes
            ]

        if isinstance(genes, list):
            return translated_genes

        sequence_constructor = cast(Any, type(genes))
        try:
            return sequence_constructor(translated_genes)
        except TypeError:
            return translated_genes

    @overload
    def translate_hypercube(
        self,
        hypercube: _HypercubeLikeT,
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
        copy: Literal[True] = True,
    ) -> _HypercubeLikeT: ...

    @overload
    def translate_hypercube(
        self,
        hypercube: "HypercubeLike",
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
        copy: Literal[False],
    ) -> None: ...

    @overload
    def translate_hypercube(
        self,
        hypercube: _HypercubeLikeT,
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
        copy: bool = True,
    ) -> Union[_HypercubeLikeT, None]: ...

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing", "copy"),
        false_policy="raise",
    )
    def translate_hypercube(
        self,
        hypercube: _HypercubeLikeT,
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
        copy: bool = True,
    ) -> Union[_HypercubeLikeT, None]:
        """
        Translate hypercube component names.

        Components are translated to one target symbol using the deterministic
        HCOP ranking. An explicitly specified component cannot be dropped
        because doing so would broaden the represented state space.

        Parameters
        ----------
        hypercube: HypercubeLike
            Hypercube or compatible component mapping to translate.
        on_unmapped: {"raise", "keep"} (default: "raise")
            Strategy for component symbols with no target ortholog. `"raise"`
            raises a `ValueError` and `"keep"` preserves the original symbol.
        copy: bool (default: True)
            Return a translated object of the same type instead of modifying
            `hypercube`. In-place translation requires a mutable mapping.

        Returns
        -------
        HypercubeLike or None
            Translated object with the same type as `hypercube` if `copy=True`;
            otherwise None.

        Raises
        ------
        TypeError
            If `hypercube` is not hypercube-like, or if `copy=False` is used
            with an immutable mapping.
        ValueError
            If a component has no target ortholog and `on_unmapped="raise"`, or
            if translation would merge explicitly specified components.
        """

        from ...logic.boolean_algebra import Hypercube
        from ...logic.boolean_algebra._typing import is_hypercube_like

        if not isinstance(hypercube, Hypercube) and not is_hypercube_like(hypercube):
            raise TypeError(
                f"unsupported argument type for 'hypercube': "
                "expected a hypercube-like mapping "
                f"but received {type(hypercube)}"
            )
        if not copy and not isinstance(hypercube, MutableMappingInstance):
            raise TypeError(
                "invalid argument combination: copy=False requires a "
                "mutable hypercube-like mapping"
            )

        on_unmapped = cast(
            Literal["raise", "keep"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep"),
            ),
        )
        mapping = {}

        for component in hypercube:
            translated = self._translate_best(
                component,
                on_unmapped="keep" if on_unmapped == "keep" else "drop",
            )
            if translated is None:
                raise ValueError(
                    "hypercube translation cannot remove an unmapped "
                    f"component: {component!r}"
                )
            mapping[component] = translated

        translated_values: Dict[str, Any] = {}
        for component, value in hypercube.items():
            translated_component = mapping[component]
            if translated_component in translated_values:
                raise ValueError("relabeling would merge hypercube components")
            translated_values[translated_component] = value

        if not copy:
            translated_hypercube = cast(Any, hypercube)
            translated_hypercube.clear()
            translated_hypercube.update(translated_values)
            return None

        hypercube_type = cast(Any, type(hypercube))
        try:
            translated_hypercube = hypercube_type(translated_values)
        except Exception as constructor_error:
            if not isinstance(hypercube, MutableMappingInstance):
                raise TypeError(
                    "cannot create a translated copy while preserving "
                    f"hypercube type {type(hypercube)}"
                ) from constructor_error

            translated_hypercube = copylib.copy(hypercube)
            translated_hypercube.clear()
            translated_hypercube.update(translated_values)

        if type(translated_hypercube) is not type(hypercube):
            raise TypeError(
                "cannot create a translated copy while preserving "
                f"hypercube type {type(hypercube)}"
            )

        return cast(_HypercubeLikeT, translated_hypercube)

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing",),
        false_policy="drop",
    )
    def translate_interaction_list(
        self,
        interaction_list: InteractionList,
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
    ) -> InteractionList:
        """
        Translate source and target genes in an interaction list.

        Edge attributes are preserved. Each source and target is translated to
        one target symbol using the deterministic HCOP ranking.

        Parameters
        ----------
        interaction_list: Sequence[Tuple[str, str, Dict[str, int]]]
            Sequence of `(source, target, attributes)` interactions.
        on_unmapped: {"raise", "keep", "drop"} (default: "raise")
            Strategy for endpoint symbols with no target ortholog. `"raise"`
            raises a `ValueError`, `"keep"` preserves the original endpoint and
            `"drop"` removes interactions containing an unmapped endpoint.

        Returns
        -------
        InteractionList
            Translated interaction list with preserved edge attributes.
        """

        on_unmapped = cast(
            Literal["raise", "keep", "drop"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep", "drop"),
            ),
        )
        translate = self._translate_best
        translated_interactions: List[Tuple[str, str, Dict[str, int]]] = []
        for source, target, attributes in interaction_list:
            translated_source = translate(
                source,
                on_unmapped,
            )
            translated_target = translate(
                target,
                on_unmapped,
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
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: Literal[True] = True,
        one_to_many: Optional[int] = None,
    ) -> pd.DataFrame: ...

    @overload
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: Literal[False],
        one_to_many: Optional[int] = None,
    ) -> None: ...

    @overload
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]: ...

    @_support_legacy_keep_if_missing(
        positional_parameters=(
            "columns",
            "keep_if_missing",
            "copy",
            "one_to_many",
        ),
        false_policy="drop",
    )
    def translate_dataframe(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]:
        """
        Translate gene symbols in selected DataFrame columns.

        Each selected value is translated to one target symbol using the
        deterministic HCOP ranking by default. If `one_to_many` is provided,
        values with up to that many orthologs expand to one row per ortholog.
        Values with more orthologs are removed.

        Examples
        --------
        >>> import pandas as pd
        >>> from bonesistools.resources.hcop import orthologs
        >>> converter = orthologs(
        ...     output_organism="mouse",
        ...     input_organism="human",
        ...     version="bundled",
        ... )
        >>> df = pd.DataFrame(
        ...     {"gene": ["EZH2", "ABCB1"], "lambda": [0.2, 0.8]}
        ... )
        >>> converter.translate_dataframe(
        ...     df,
        ...     columns="gene",
        ...     one_to_many=2,
        ... )
             gene  lambda
        0    Ezh2     0.2
        1  Abcb1a     0.8
        2  Abcb1b     0.8

        Parameters
        ----------
        df: pandas.DataFrame
            DataFrame containing symbols from `input_organism`.
        columns: str or sequence of str, optional
            Columns to translate. If `None`, translate any existing column among
            `source`, `target` and `genesymbol`.
        on_unmapped: {"raise", "keep", "drop", "nan"} (default: "raise")
            Strategy for gene symbols with no target ortholog.

            - `"raise"`: raise a `ValueError` if any symbol is unmapped.
            - `"keep"`: preserve the original symbol.
            - `"drop"`: remove rows containing an unmapped symbol.
            - `"nan"`: replace unmapped symbols with `pandas.NA`.
        copy: bool (default: True)
            Return a translated copy instead of modifying `df`.
        one_to_many: int, optional
            Maximum number of orthologs allowed per gene. Each accepted
            ortholog produces one output row.

        Returns
        -------
        pandas.DataFrame or None
            Translated DataFrame if `copy=True`; otherwise None.

        Raises
        ------
        ValueError
            If none of the requested columns are present in `df`, if
            `on_unmapped` is invalid, or if an unmapped symbol is encountered
            with `on_unmapped="raise"`.
        """

        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"unsupported argument type for 'df': "
                f"expected {pd.DataFrame} but received {type(df)}"
            )

        on_unmapped = cast(
            Literal["raise", "keep", "drop", "nan"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep", "drop", "nan"),
            ),
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
                    on_unmapped=on_unmapped,
                )
            translated = translated.reset_index(drop=True)
            return translated

        for column in selected_columns:
            translated[column] = self._translate_series(
                cast(pd.Series, translated[column]),
                on_unmapped=on_unmapped,
            )
            if on_unmapped == "drop":
                translated.dropna(subset=[column], inplace=True)

        translated.reset_index(drop=True, inplace=True)
        if copy:
            return translated

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: Literal[True] = True,
        one_to_many: Optional[int] = None,
    ) -> pd.DataFrame: ...

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: Literal[False],
        one_to_many: Optional[int] = None,
    ) -> None: ...

    @overload
    def translate_df(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]: ...

    @_deprecated(replacement="`Orthologs.translate_dataframe()`")
    @_support_legacy_keep_if_missing(
        positional_parameters=(
            "columns",
            "keep_if_missing",
            "copy",
            "one_to_many",
        ),
        false_policy="drop",
    )
    def translate_df(
        self,
        df: pd.DataFrame,
        *,
        columns: Optional[Union[str, Sequence[str]]] = None,
        on_unmapped: Literal["raise", "keep", "drop", "nan"] = "raise",
        copy: bool = True,
        one_to_many: Optional[int] = None,
    ) -> Union[pd.DataFrame, None]:
        """Deprecated alias for `translate_dataframe()`."""

        return self.translate_dataframe(
            df,
            columns=columns,
            on_unmapped=on_unmapped,
            copy=copy,
            one_to_many=one_to_many,
        )

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing", "copy"),
        false_policy="drop",
    )
    def translate_graph(
        self,
        graph: Graph[Any],
        *,
        on_unmapped: Literal["raise", "keep", "drop"] = "raise",
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
        on_unmapped: {"raise", "keep", "drop"} (default: "raise")
            Strategy for node symbols with no target ortholog. `"raise"` raises
            a `ValueError`, `"keep"` preserves the original node label and
            `"drop"` removes unmapped nodes and their incident edges.
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

        on_unmapped = cast(
            Literal["raise", "keep", "drop"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep", "drop"),
            ),
        )
        translated = graph.copy() if copy else graph
        mapping = {}
        missing = []
        for node in translated.nodes:
            target = self._translate_best(
                node,
                on_unmapped=on_unmapped,
            )
            if target is None:
                missing.append(node)
            else:
                mapping[node] = target

        translated.remove_nodes_from(missing)
        nx.relabel_nodes(translated, mapping=mapping, copy=False)

        if copy:
            return translated

    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing", "copy"),
        false_policy="raise",
    )
    def translate_boolean_network(
        self,
        bn: "BooleanNetworkLike",
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
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
        on_unmapped: {"raise", "keep"} (default: "raise")
            Strategy for component symbols with no target ortholog. `"raise"`
            raises a `ValueError` and `"keep"` preserves the original symbol.
        copy: bool (default: False)
            Return a translated copy instead of modifying `bn`.

        Returns
        -------
        BooleanNetworkLike or None
            Translated Boolean network if `copy=True`; otherwise None.

        Raises
        ------
        ValueError
            If a component has no target ortholog and `on_unmapped="raise"`.
        """

        from ...logic.boolean_network import BooleanNetwork
        from ...logic.boolean_network._typing import is_boolean_network_like

        if not is_boolean_network_like(bn):
            raise TypeError(
                "unsupported argument type for 'bn': "
                "expected Boolean network-like object "
                f"but received {type(bn)}"
            )

        on_unmapped = cast(
            Literal["raise", "keep"],
            _validate_on_unmapped(
                on_unmapped,
                allowed=("raise", "keep"),
            ),
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
            translated = self._translate_best(
                gene,
                on_unmapped="keep" if on_unmapped == "keep" else "drop",
            )
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
    @_support_legacy_keep_if_missing(
        positional_parameters=("keep_if_missing", "copy"),
        false_policy="raise",
    )
    def translate_bn(
        self,
        bn: "BooleanNetworkLike",
        *,
        on_unmapped: Literal["raise", "keep"] = "raise",
        copy: bool = False,
    ) -> Union["BooleanNetworkLike", None]:
        """Deprecated alias for `translate_boolean_network()`."""

        return self.translate_boolean_network(
            bn,
            on_unmapped=on_unmapped,
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

    def _translate_best(
        self,
        gene: str,
        on_unmapped: Literal["raise", "keep", "drop"],
    ) -> Optional[str]:

        if "_" in gene:
            translated_subunits: List[str] = []
            for subunit in gene.split("_"):
                translated = self._best_mapping.get(subunit)
                if translated is None:
                    if on_unmapped == "raise":
                        raise ValueError(
                            f"no target ortholog found for gene symbol {subunit!r}"
                        )
                    if on_unmapped == "drop":
                        return None
                    translated = subunit
                translated_subunits.append(translated)
            return "_".join(translated_subunits)

        translated = self._best_mapping.get(gene)
        if translated is not None:
            return translated
        if on_unmapped == "raise":
            raise ValueError(f"no target ortholog found for gene symbol {gene!r}")
        return gene if on_unmapped == "keep" else None

    def _translate_series(
        self,
        values: pd.Series,
        on_unmapped: Literal["raise", "keep", "drop", "nan"],
    ) -> pd.Series:
        """Translate a Series while preserving complex-name semantics."""

        source = values.astype(str)
        best_mapping = self._best_mapping
        translate = self._translate_best
        if on_unmapped == "keep":
            translated = [
                translate(gene, on_unmapped="keep")
                if "_" in gene
                else best_mapping.get(gene, gene)
                for gene in source
            ]
        elif on_unmapped in {"drop", "nan"}:
            translated = [
                translate(gene, on_unmapped="drop")
                if "_" in gene
                else best_mapping.get(gene)
                for gene in source
            ]
            translated = [pd.NA if gene is None else gene for gene in translated]
        else:
            translated = [translate(gene, on_unmapped="raise") for gene in source]
        return pd.Series(translated, index=values.index, name=values.name)

    def _translate_df_one_to_many(
        self,
        df: pd.DataFrame,
        column: str,
        one_to_many: int,
        on_unmapped: Literal["raise", "keep", "drop", "nan"],
    ) -> pd.DataFrame:

        expansions = self._ortholog_expansions(
            values=cast(pd.Series, df[column]),
            one_to_many=one_to_many,
            on_unmapped=on_unmapped,
        )
        expanded = cast(pd.Series, df[column]).astype(str).map(expansions)
        retained = expanded.map(bool)
        translated = cast(pd.DataFrame, df.loc[retained].copy())
        translated[column] = expanded.loc[retained]
        return cast(pd.DataFrame, translated.explode(column))

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
        return cast(
            pd.DataFrame,
            paths[
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
            ],
        )

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
        on_unmapped: Literal["raise", "keep", "drop", "nan"],
    ) -> Dict[str, List[Any]]:
        """Return row-ordered target expansions for distinct source values."""

        map_dict = self._one_to_many_mapping()
        expansions: Dict[str, List[Any]] = {}
        for value in values.drop_duplicates():
            source = str(value)
            subunit_lists: List[List[str]] = []
            has_excess = False
            has_unmapped = False
            for subunit in source.split("_"):
                targets = map_dict.get(subunit, [])
                if len(targets) == 0:
                    if on_unmapped == "raise":
                        raise ValueError(
                            f"no target ortholog found for gene symbol {subunit!r}"
                        )
                    has_unmapped = True
                    if on_unmapped == "keep":
                        subunit_lists.append([subunit])
                    else:
                        subunit_lists.append([])
                elif len(targets) > one_to_many:
                    has_excess = True
                    subunit_lists.append([])
                else:
                    subunit_lists.append([str(target) for target in targets])

            if has_excess or (has_unmapped and on_unmapped == "drop"):
                expansions[source] = []
                continue
            if has_unmapped and on_unmapped == "nan":
                expansions[source] = [pd.NA]
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
