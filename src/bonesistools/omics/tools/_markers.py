#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Collection as CollectionInstance
from collections.abc import Mapping as MappingInstance
from collections.abc import Sequence as Seq
from inspect import signature
from typing import (
    Any,
    Callable,
    Collection,
    Iterable,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
    cast,
    overload,
)

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy import log2
from scipy.sparse import issparse
from scipy.stats import hypergeom, ks_2samp

from ..._compat import Literal
from ..._validation import _as_boolean, _as_callable, _as_literal, _as_probability
from ..._warnings import _warn_deprecated
from .._typing import Matrix, VarSubset, anndata_checker
from ._conversion import to_dataframe
from ._stats import (
    CORRECTION_METHODS,
    CorrectionMethod,
    _adjust_pvalues,
    _resolve_background_groups,
    _resolve_groups,
    welch_tests,
    wilcoxon_tests,
)
from ._utils import _as_dense_matrix_chunk, _get_expression_with_gene_names

_CorrMethod = Literal["benjamini-hochberg", "bonferroni"]
_Alternatives = Literal["two-sided", "less", "greater"]
DEAMethod = Literal["welch", "welch_overestimate", "wilcoxon"]
SignatureCollection = Union[
    Mapping[str, Collection[str]],
    Sequence[Tuple[str, Collection[str]]],
]
DEA_METHODS: Tuple[DEAMethod, ...] = (
    "welch",
    "welch_overestimate",
    "wilcoxon",
)
_SMIRNOV_CHUNK_SIZE = 64


class _SmirnovInputs(NamedTuple):
    expression_mtx: Matrix
    gene_names: np.ndarray
    groups: Tuple[Any, ...]
    group_indices: Tuple[np.ndarray, ...]
    reference_indices: Tuple[np.ndarray, ...]
    unique_gene_names: bool


class _SmirnovResults(NamedTuple):
    statistics: np.ndarray
    locations: np.ndarray
    signs: np.ndarray
    pvals: np.ndarray


def _ks_supports_axis() -> bool:
    """Return whether the installed SciPy supports vectorized KS tests."""

    try:
        return "axis" in signature(ks_2samp).parameters
    except (TypeError, ValueError):
        return False


def _ks_returns_details() -> bool:
    """Return whether KS results include statistic location and sign."""

    result = ks_2samp(np.asarray([0.0]), np.asarray([1.0]))
    return hasattr(result, "statistic_location") and hasattr(
        result,
        "statistic_sign",
    )


_KS_SUPPORTS_AXIS = _ks_supports_axis()
_KS_RETURNS_DETAILS = _ks_returns_details()


@anndata_checker
def logfoldchanges(
    adata: AnnData,
    groupby: str,
    layer: Optional[str] = None,
    column_name: str = "logfoldchanges",
    is_log: bool = False,
    cluster_rebalancing: bool = False,
    filter_logfoldchanges: Optional[Callable[[np.ndarray], Any]] = None,
) -> pd.DataFrame:
    """
    Compute log2 fold changes between each group and the remaining observations.

    The computation follows the fold-change convention used by Seurat's
    FindAllMarkers rather than Scanpy's `rank_gene_groups` output.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation key in `adata.obs` defining groups.
    layer: str, optional
        Layer to use instead of `adata.X`.
    column_name: str (default: "logfoldchanges")
        Column name in output dataframe storing log2 fold-change values.
    is_log: bool (default: False)
        Whether values are already log1p-transformed.
    cluster_rebalancing: bool (default: False)
        Whether to average groups before computing fold changes, so groups
        contribute equally regardless of size.
    filter_logfoldchanges: Callable, optional
        Callable used to filter rows from the resulting log-fold-change table.

    Returns
    -------
    DataFrame
        Table with `group`, `names` and `<column_name>` columns.

    See Also
    --------
    Get more information about the difference between log2 fold-changes derived
    with Seurat and Scanpy here:
    <https://www.biostars.org/p/453129/>

    """

    def compute_logfc(mean_in, mean_out, cluster):
        with np.errstate(divide="ignore", invalid="ignore"):
            logfc = log2(mean_in) - log2(mean_out)

        __df = pd.DataFrame(
            logfc,
            columns=cast(Any, [column_name]),
        )
        __df.reset_index(names="names", inplace=True)
        __df.insert(0, "group", cluster)
        return __df

    logfoldchanges = []
    counts_df = cast(
        pd.DataFrame,
        to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log),
    )

    if cluster_rebalancing:
        mean_counts_df = cast(
            pd.DataFrame,
            counts_df.groupby(
                by=groupby,
                sort=True,
                observed=False,
            ).mean(),
        )

        for cluster in sorted(adata.obs[groupby].dropna().unique()):
            _mean_in = mean_counts_df.loc[cluster]
            _mean_out = mean_counts_df.drop(index=cluster, inplace=False).mean()

            logfoldchanges.append(compute_logfc(_mean_in, _mean_out, cluster))

    else:
        for cluster in sorted(adata.obs[groupby].dropna().unique()):
            _mean_in = counts_df.loc[
                counts_df[groupby] == cluster,
                counts_df.columns != groupby,
            ].mean()

            _mean_out = counts_df.loc[
                counts_df[groupby] != cluster,
                counts_df.columns != groupby,
            ].mean()

            logfoldchanges.append(compute_logfc(_mean_in, _mean_out, cluster))

    if len(logfoldchanges) == 0:
        return pd.DataFrame(columns=cast(Any, ["group", "names", column_name]))

    logfoldchanges_df = pd.concat(logfoldchanges, ignore_index=True)

    filter_logfoldchanges = _as_callable(
        filter_logfoldchanges,
        "filter_logfoldchanges",
        allow_none=True,
    )
    if filter_logfoldchanges is not None:
        logfoldchanges_df = logfoldchanges_df.loc[
            filter_logfoldchanges(
                cast(np.ndarray, logfoldchanges_df[column_name].values)
            )
        ]

    return logfoldchanges_df.reset_index(drop=True)


def calculate_logfoldchanges(*args: Any, **kwargs: Any) -> pd.DataFrame:
    """
    Deprecated alias for `logfoldchanges`.
    """

    _warn_deprecated(
        "`bt.omics.tl.calculate_logfoldchanges`",
        replacement="`bt.omics.tl.logfoldchanges`",
        stacklevel=2,
    )
    return logfoldchanges(*args, **kwargs)


@anndata_checker
def dea(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Sequence[Any]] = "all",
    background: Union[Literal["rest"], Sequence[Any]] = "rest",
    method: DEAMethod = "welch",
    expression: Optional[str] = None,
    is_log: bool = False,
    var_subset: VarSubset = None,
    correction: Optional[CorrectionMethod] = "benjamini-hochberg",
    alpha: Optional[float] = 0.05,
    filter_logfoldchanges: Optional[Callable[[np.ndarray], Any]] = None,
    max_memory: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    """
    Run differential expression analysis.

    Differential expression analysis combines statistical testing and log2
    fold-change estimation.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation column in `adata.obs` defining groups.
    groups: 'all' or sequence (default: "all")
        Groups to test. If `"all"`, test all groups in `adata.obs[groupby]`;
        when `background` is a sequence, background groups themselves are not
        tested.
    background: 'rest' or sequence (default: "rest")
        Background population. If `"rest"`, each target group is compared
        against all observations outside that group. If a sequence is provided,
        each target group is compared against the union of observations
        belonging to those background labels. To compare against an actual group
        named `"rest"`, pass `background=["rest"]`.
    method: {'welch', 'welch_overestimate', 'wilcoxon'} (default: 'welch')
        Statistical test used before log2 fold-change estimation.
    expression: str, optional
        Expression source. If `None` or `"X"`, use `adata.X`. If `"raw.X"`, use
        `adata.raw.X`. Otherwise, interpret as a layer key in `adata.layers`.
    is_log: bool (default: False)
        Whether selected expression values are already log1p-transformed. This
        only affects log2 fold-change estimation, not statistical testing.
    var_subset: str or collection of str, optional
        Variables to test. If a string is provided, it is interpreted as a
        boolean column in `adata.var`. If a collection is provided, it is
        interpreted as variable names. If `None`, test all variables in the
        selected expression matrix.
    correction: {'benjamini-hochberg', 'bonferroni'} or None
        Multiple-testing correction applied independently for each group. If
        None, p-values are not adjusted.
    alpha: float, optional
        Significance threshold applied to adjusted p-values (`pvals_adj`). If
        None, p-value filtering is disabled.
    filter_logfoldchanges: Callable, optional
        Callable used to filter rows from the resulting table after p-value
        filtering. The callable receives the log2 fold-change vector and must
        return a boolean mask.
    max_memory: int or str, optional
        Approximate maximum memory allocated to dense working arrays. If `None`,
        process all variables in a single chunk. Integers are interpreted as
        bytes. Human-readable strings such as `"512MB"`, `"2GB"` or `"1GiB"`
        are accepted.

    Returns
    -------
    DataFrame
        Differential expression table. The table contains:

        - `feature`: variable name;
        - `group`: tested group;
        - `statistics`: primary test statistic;
        - `pvals`: two-sided p-values;
        - `pvals_adj`: adjusted p-values;
        - `logfoldchanges`: log2 fold changes.

        The `statistics` column contains the primary test statistic: Welch t
        statistic for `method="welch"` and `method="welch_overestimate"`, and z
        statistic for `method="wilcoxon"`.

    See Also
    --------
    welch_tests
    wilcoxon_tests
    logfoldchanges
    """

    method = _as_literal(method, choices=DEA_METHODS, name="method")
    alpha = None if alpha is None else _as_probability(alpha, "alpha")
    filter_logfoldchanges = _as_callable(
        filter_logfoldchanges,
        "filter_logfoldchanges",
        allow_none=True,
    )

    if method == "welch":
        test_results = welch_tests(
            adata,
            groupby=groupby,
            groups=groups,
            background=background,
            expression=expression,
            var_subset=var_subset,
            correction=correction,
            overestimate_variance=False,
            max_memory=max_memory,
        )
    elif method == "welch_overestimate":
        test_results = welch_tests(
            adata,
            groupby=groupby,
            groups=groups,
            background=background,
            expression=expression,
            var_subset=var_subset,
            correction=correction,
            overestimate_variance=True,
            max_memory=max_memory,
        )
    else:
        test_results = wilcoxon_tests(
            adata,
            groupby=groupby,
            groups=groups,
            background=background,
            expression=expression,
            var_subset=var_subset,
            correction=correction,
            max_memory=max_memory,
        )

    if test_results.empty:
        return pd.DataFrame(
            columns=cast(
                Any,
                [
                    "feature",
                    "group",
                    "statistics",
                    "pvals",
                    "pvals_adj",
                    "logfoldchanges",
                ],
            )
        )

    test_results_df = test_results.reset_index().rename(columns={"names": "feature"})
    logfoldchanges_df = _dea_logfoldchanges(
        adata,
        groupby,
        groups,
        background,
        expression,
        is_log,
        var_subset,
    ).rename(columns={"names": "feature"})
    dea_df = test_results_df.merge(
        logfoldchanges_df,
        on=["group", "feature"],
        how="left",
        sort=False,
    )

    first_columns = [
        "feature",
        "group",
        "statistics",
        "pvals",
        "pvals_adj",
        "logfoldchanges",
    ]
    dea_df = dea_df[first_columns]

    if alpha is not None:
        dea_df = cast(pd.DataFrame, dea_df.loc[dea_df["pvals_adj"] <= alpha])

    if filter_logfoldchanges is not None:
        dea_df = cast(
            pd.DataFrame,
            dea_df.loc[
                filter_logfoldchanges(cast(np.ndarray, dea_df["logfoldchanges"].values))
            ],
        )

    return cast(pd.DataFrame, dea_df.reset_index(drop=True))


def ora(
    query_set: Collection[str],
    signatures: SignatureCollection,
    background: Collection[str],
    correction: CorrectionMethod = "benjamini-hochberg",
    include_overlap: bool = False,
) -> pd.DataFrame:
    """
    Run over-representation analysis (ORA).

    Over-representation analysis (ORA) is a statistical method used to identify
    biological signatures that are significantly enriched in a query gene set
    relative to a background gene universe.

    For each tested signature, enrichment significance is evaluated using a
    one-sided hypergeometric test.

    Genes absent from `background` are discarded from both `query_set` and each
    tested signature before enrichment testing.

    Parameters
    ----------
    query_set: collection of str
        Genes of interest.
    signatures: mapping or sequence of tuple
        Named gene signatures. Provide either a mapping from signature names to
        genes, or a sequence of `(name, genes)` pairs.
    background: collection of str
        Tested gene universe.
    correction: {'benjamini-hochberg', 'bonferroni'}
        Multiple-testing correction used to compute `pvals_adj`.
    include_overlap: bool (default: False)
        Whether to include the sorted tuple of overlapping genes for each
        signature.

    Returns
    -------
    DataFrame
        Table with one row per tested signature. The table contains:

        - `pvals`: one-sided hypergeometric p-value;
        - `pvals_adj`: adjusted p-value;
        - `observed_overlap`: number of query genes in the signature;
        - `expected_overlap`: expected overlap under the null hypothesis;
        - `fold_enrichment`: observed over expected overlap.
        - `signature_size`: number of tested signature genes;

        If `include_overlap=True`, the table also contains:

        - `overlap`: sorted tuple of overlapping genes.

        The numbers of tested query and background genes are stored in
        `DataFrame.attrs["query_size"]` and
        `DataFrame.attrs["background_size"]`.

    """

    correction = _as_literal(
        correction,
        choices=CORRECTION_METHODS,
        name="correction",
    )
    include_overlap = _as_boolean(include_overlap, "include_overlap")

    background_set = _as_gene_set(background, "background")
    if not background_set:
        raise ValueError(
            "invalid argument value for 'background': expected at least one gene"
        )

    query = _as_gene_set(query_set, "query_set") & background_set
    if not query:
        raise ValueError(
            "invalid argument value for 'query_set': "
            "no query genes remain after intersection with background"
        )

    signature_items = _as_signature_items(signatures)

    rows = []
    background_size = len(background_set)
    query_size = len(query)

    for signature_name, signature_genes in signature_items:
        signature = _as_gene_set(signature_genes, f"signatures[{signature_name!r}]")
        tested_signature = signature & background_set

        if not tested_signature:
            continue

        overlap = query & tested_signature
        overlap_size = len(overlap)
        signature_size = len(tested_signature)

        pval = _hypergeometric_pvalue(
            overlap_size=overlap_size,
            signature_size=signature_size,
            query_size=query_size,
            background_size=background_size,
        )

        expected_overlap = query_size * signature_size / background_size
        fold_enrichment = (
            overlap_size / expected_overlap if expected_overlap > 0 else np.nan
        )

        rows.append(
            {
                "signature": signature_name,
                "pvals": pval,
                "observed_overlap": overlap_size,
                "expected_overlap": expected_overlap,
                "fold_enrichment": fold_enrichment,
                "signature_size": signature_size,
                **({"overlap": tuple(sorted(overlap))} if include_overlap else {}),
            }
        )

    columns = [
        "pvals",
        "pvals_adj",
        "observed_overlap",
        "expected_overlap",
        "fold_enrichment",
        "signature_size",
    ]
    if include_overlap:
        columns.append("overlap")
    df = pd.DataFrame(rows)
    if df.empty:
        empty_df = pd.DataFrame(columns=cast(Any, columns))
        empty_df.index.name = "signature"
        empty_df.attrs["query_size"] = query_size
        empty_df.attrs["background_size"] = background_size
        return empty_df

    df["pvals_adj"] = _adjust_pvalues(
        df["pvals"].to_numpy(dtype=float),
        correction,
    )
    df = cast(pd.DataFrame, df.set_index("signature"))
    df = cast(pd.DataFrame, df[columns])
    df.attrs["query_size"] = query_size
    df.attrs["background_size"] = background_size

    return cast(
        pd.DataFrame,
        df.sort_values(
            by=["pvals", "observed_overlap"],
            ascending=[True, False],
            kind="mergesort",
        ),
    )


def _dea_logfoldchanges(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Sequence[Any]],
    background: Union[Literal["rest"], Sequence[Any]],
    expression: Optional[str],
    is_log: bool,
    var_subset: VarSubset,
) -> pd.DataFrame:

    expression_mtx, gene_names = _get_expression_with_gene_names(
        adata,
        expression,
        var_subset,
    )
    labels = cast(pd.Series, adata.obs[groupby])
    background_groups = _resolve_background_groups(labels, background)
    target_groups = _resolve_groups(labels, groups, background_groups)

    logfoldchange_tables = []
    for group in target_groups:
        group_mask = np.asarray(labels == group, dtype=bool)
        if background_groups is None:
            background_mask = ~group_mask
        else:
            background_mask = np.asarray(labels.isin(background_groups), dtype=bool)

        combined_mask = group_mask | background_mask
        temporary_obs_index = pd.Index(labels.index.to_numpy()[combined_mask])
        temporary_adata = AnnData(
            X=cast(Any, expression_mtx)[combined_mask, :],
            obs=pd.DataFrame(
                {"__group": group_mask[combined_mask]},
                index=temporary_obs_index,
            ),
            var=pd.DataFrame(index=gene_names),
        )
        group_logfoldchanges = logfoldchanges(
            temporary_adata,
            groupby="__group",
            is_log=is_log,
        )
        group_logfoldchanges = group_logfoldchanges.loc[
            cast(pd.Series, group_logfoldchanges["group"]).astype(bool),
            ["names", "logfoldchanges"],
        ].copy()
        group_logfoldchanges.insert(0, "group", group)
        logfoldchange_tables.append(group_logfoldchanges)

    if not logfoldchange_tables:
        return pd.DataFrame(columns=cast(Any, ["group", "names", "logfoldchanges"]))

    return pd.concat(logfoldchange_tables, ignore_index=True)


def _hypergeometric_pvalue(
    overlap_size: int,
    signature_size: int,
    query_size: int,
    background_size: int,
) -> float:
    """
    Compute `P[X >= overlap_size]` for over-representation analysis.

    `X` follows a hypergeometric distribution with population size
    `background_size`, number of success states `signature_size`, and number of
    draws `query_size`.
    """

    return float(
        hypergeom.sf(
            overlap_size - 1,
            background_size,
            signature_size,
            query_size,
        )
    )


def _as_gene_set(genes: Collection[str], name: str) -> Set[str]:

    if isinstance(genes, str) or not isinstance(genes, CollectionInstance):
        raise TypeError(
            f"unsupported argument type for {name!r}: "
            f"expected a collection of {str} but received {type(genes)}"
        )

    gene_set = set(genes)
    invalid = [gene for gene in gene_set if not isinstance(gene, str)]
    if invalid:
        gene = invalid[0]
        raise TypeError(
            f"unsupported element type in {name!r}: "
            f"expected {str} but received {type(gene)}"
        )

    return gene_set


def _as_signature_items(
    signatures: SignatureCollection,
) -> List[Tuple[str, Collection[str]]]:

    if isinstance(signatures, MappingInstance):
        items = list(signatures.items())
    elif isinstance(signatures, Seq) and not isinstance(signatures, str):
        items = list(signatures)
    else:
        raise TypeError(
            "unsupported argument type for 'signatures': "
            "expected a mapping or a sequence of (name, genes) pairs"
        )

    if not items:
        raise ValueError(
            "invalid argument value for 'signatures': expected at least one signature"
        )

    resolved = []
    for item in items:
        if (
            isinstance(item, str)
            or not isinstance(item, Seq)
            or len(item) != 2
            or not isinstance(item[0], str)
        ):
            raise TypeError(
                "unsupported element type in 'signatures': "
                "expected (str, Collection[str]) pairs"
            )

        name, genes = item
        resolved.append((name, genes))

    names = [name for name, _ in resolved]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(
            "invalid argument value for 'signatures': duplicate signature name(s) "
            + ", ".join(repr(name) for name in duplicates)
        )

    return resolved


def _prepare_smirnov_inputs(
    adata: AnnData,
    groupby: str,
    groups: Sequence[Any],
    reference: str,
    layer: Optional[str],
) -> _SmirnovInputs:
    """Prepare expression data and observation indices for KS tests."""

    expression_mtx = adata.layers[layer] if layer is not None else adata.X
    if issparse(expression_mtx):
        expression_mtx = cast(
            Matrix,
            cast(Any, expression_mtx).tocsc(copy=True),
        )
        cast(Any, expression_mtx).eliminate_zeros()
    else:
        expression_mtx = cast(Matrix, np.asarray(expression_mtx))

    labels = adata.obs[groupby].to_numpy()
    resolved_groups = tuple(groups)
    group_indices = tuple(np.flatnonzero(labels == group) for group in resolved_groups)
    if reference == "rest":
        reference_indices = tuple(
            np.flatnonzero(labels != group) for group in resolved_groups
        )
    else:
        fixed_reference = np.flatnonzero(labels == reference)
        reference_indices = (fixed_reference,) * len(resolved_groups)

    return _SmirnovInputs(
        expression_mtx=expression_mtx,
        gene_names=np.asarray(adata.var_names),
        groups=resolved_groups,
        group_indices=group_indices,
        reference_indices=reference_indices,
        unique_gene_names=adata.var_names.is_unique,
    )


def _empty_smirnov_results(inputs: _SmirnovInputs) -> _SmirnovResults:
    """Allocate arrays for one KS result per gene and group."""

    shape = (inputs.gene_names.size, len(inputs.groups))
    matrix_dtype = np.dtype(cast(Any, inputs.expression_mtx).dtype)
    has_empty_sample = any(
        indices.size == 0 for indices in inputs.group_indices + inputs.reference_indices
    )
    probe_dtype = (
        np.dtype(float)
        if has_empty_sample and not np.issubdtype(matrix_dtype, np.floating)
        else matrix_dtype
    )
    probe = cast(
        Any,
        ks_2samp(
            np.asarray([0], dtype=probe_dtype),
            np.asarray([1], dtype=probe_dtype),
        ),
    )
    location = getattr(
        probe,
        "statistic_location",
        np.asarray(0, dtype=probe_dtype),
    )
    return _SmirnovResults(
        statistics=np.empty(shape, dtype=np.asarray(probe.statistic).dtype),
        locations=np.empty(shape, dtype=np.asarray(location).dtype),
        signs=np.empty(shape, dtype=float),
        pvals=np.empty(shape, dtype=np.asarray(probe.pvalue).dtype),
    )


def _finalize_smirnov_results(results: _SmirnovResults) -> _SmirnovResults:
    """Use SciPy's compact integer dtype when all KS signs are defined."""

    signs = results.signs
    if signs.size and not np.isnan(signs).any():
        signs = signs.astype(np.int8)
    return _SmirnovResults(
        statistics=results.statistics,
        locations=results.locations,
        signs=signs,
        pvals=results.pvals,
    )


def _vectorized_smirnov_results(
    inputs: _SmirnovInputs,
    alternative: _Alternatives,
) -> _SmirnovResults:
    """Compute KS tests by vectorizing over contiguous gene chunks."""

    results = _empty_smirnov_results(inputs)
    n_genes = inputs.gene_names.size
    for start in range(0, n_genes, _SMIRNOV_CHUNK_SIZE):
        end = min(start + _SMIRNOV_CHUNK_SIZE, n_genes)
        values = _as_dense_matrix_chunk(inputs.expression_mtx, start, end)
        for group_index, (sample_indices, reference_indices) in enumerate(
            zip(inputs.group_indices, inputs.reference_indices)
        ):
            ks = cast(
                Any,
                ks_2samp(
                    values[sample_indices],
                    values[reference_indices],
                    alternative=alternative,
                    axis=0,
                ),
            )
            results.statistics[start:end, group_index] = ks.statistic
            results.locations[start:end, group_index] = ks.statistic_location
            results.signs[start:end, group_index] = ks.statistic_sign
            results.pvals[start:end, group_index] = ks.pvalue

    return _finalize_smirnov_results(results)


def _smirnov_gene_indices(inputs: _SmirnovInputs) -> Tuple[np.ndarray, ...]:
    """Return source columns represented by each output gene name."""

    if inputs.unique_gene_names:
        return tuple(
            np.asarray([gene_index], dtype=int)
            for gene_index in range(inputs.gene_names.size)
        )

    indices_by_name = {}
    for gene_index, gene_name in enumerate(inputs.gene_names):
        indices_by_name.setdefault(gene_name, []).append(gene_index)
    return tuple(
        np.asarray(indices_by_name[gene_name], dtype=int)
        for gene_name in inputs.gene_names
    )


def _smirnov_gene_values(
    inputs: _SmirnovInputs,
    gene_indices: np.ndarray,
) -> np.ndarray:
    """Extract all observation values associated with one gene name."""

    values = cast(Any, inputs.expression_mtx)[:, gene_indices]
    if issparse(values):
        values = values.toarray()
    return np.asarray(values)


def _scalar_smirnov_results(
    inputs: _SmirnovInputs,
    alternative: _Alternatives,
) -> _SmirnovResults:
    """Compute detailed KS results one gene and group at a time."""

    results = _empty_smirnov_results(inputs)
    for gene_index, source_indices in enumerate(_smirnov_gene_indices(inputs)):
        values = _smirnov_gene_values(inputs, source_indices)
        for group_index, (sample_indices, reference_indices) in enumerate(
            zip(inputs.group_indices, inputs.reference_indices)
        ):
            ks = cast(
                Any,
                ks_2samp(
                    values[sample_indices].reshape(-1),
                    values[reference_indices].reshape(-1),
                    alternative=alternative,
                ),
            )
            results.statistics[gene_index, group_index] = ks.statistic
            results.locations[gene_index, group_index] = ks.statistic_location
            results.signs[gene_index, group_index] = ks.statistic_sign
            results.pvals[gene_index, group_index] = ks.pvalue

    return _finalize_smirnov_results(results)


def _legacy_smirnov_details(
    sample: np.ndarray,
    reference: np.ndarray,
    statistic: float,
    alternative: _Alternatives,
) -> Tuple[Any, Union[int, float]]:
    """Recover the KS statistic location and sign for older SciPy releases."""

    if sample.size == 0 or reference.size == 0 or np.isnan(statistic):
        return np.nan, np.nan

    sample = np.sort(sample)
    reference = np.sort(reference)
    pooled = np.concatenate((sample, reference))
    values = np.sort(pooled)
    ranks = np.searchsorted(values, pooled, side="left") + 1
    total_size = pooled.size
    sample_counts = np.diff(
        np.concatenate(([1], ranks[: sample.size], [total_size + 1]))
    )
    reference_counts = np.diff(
        np.concatenate(([1], ranks[sample.size :], [total_size + 1]))
    )
    probability_dtype = np.result_type(sample.dtype, reference.dtype, np.float32)
    sample_cdf = np.repeat(
        np.linspace(0, 1, sample.size + 1, dtype=probability_dtype),
        sample_counts,
    )
    reference_cdf = np.repeat(
        np.linspace(0, 1, reference.size + 1, dtype=probability_dtype),
        reference_counts,
    )
    differences = sample_cdf - reference_cdf
    minimum_index = int(np.argmin(differences))
    maximum_index = int(np.argmax(differences))
    minimum = float(np.clip(-differences[minimum_index], 0, 1))
    maximum = float(differences[maximum_index])

    if alternative == "less" or (alternative == "two-sided" and minimum > maximum):
        return values[minimum_index], -1
    return values[maximum_index], 1


def _legacy_smirnov_results(
    inputs: _SmirnovInputs,
    alternative: _Alternatives,
) -> _SmirnovResults:
    """Compute KS tests and reconstruct details missing from older SciPy."""

    results = _empty_smirnov_results(inputs)
    for gene_index, source_indices in enumerate(_smirnov_gene_indices(inputs)):
        values = _smirnov_gene_values(inputs, source_indices)
        for group_index, (sample_indices, reference_indices) in enumerate(
            zip(inputs.group_indices, inputs.reference_indices)
        ):
            sample = values[sample_indices].reshape(-1)
            reference = values[reference_indices].reshape(-1)
            ks = cast(
                Any,
                ks_2samp(
                    sample,
                    reference,
                    alternative=alternative,
                ),
            )
            location, sign = _legacy_smirnov_details(
                sample,
                reference,
                float(ks.statistic),
                alternative,
            )
            results.statistics[gene_index, group_index] = ks.statistic
            results.locations[gene_index, group_index] = location
            results.signs[gene_index, group_index] = sign
            results.pvals[gene_index, group_index] = ks.pvalue

    return _finalize_smirnov_results(results)


def _format_smirnov_results(
    inputs: _SmirnovInputs,
    results: _SmirnovResults,
    corr_method: _CorrMethod,
    pval_cutoff: Optional[float],
) -> pd.DataFrame:
    """Build, sort and adjust the public Smirnov result table."""

    n_genes = inputs.gene_names.size
    if n_genes == 0 or not inputs.groups:
        dataframe = pd.DataFrame(
            columns=cast(
                Any,
                [
                    "group",
                    "names",
                    "statistics",
                    "locations",
                    "signs",
                    "pvals",
                ],
            )
        )
    else:
        flattened = (
            results.statistics.reshape(-1),
            results.locations.reshape(-1),
            results.signs.reshape(-1),
            results.pvals.reshape(-1),
        )
        has_nan = any(
            np.issubdtype(values.dtype, np.floating) and np.isnan(values).any()
            for values in flattened
        )
        numeric_values = (
            tuple(pd.Series(values.tolist(), dtype=object) for values in flattened)
            if has_nan
            else flattened
        )
        dataframe = pd.DataFrame(
            {
                "group": np.tile(np.asarray(inputs.groups), n_genes),
                "names": np.repeat(inputs.gene_names, len(inputs.groups)),
                "statistics": numeric_values[0],
                "locations": numeric_values[1],
                "signs": numeric_values[2],
                "pvals": numeric_values[3],
            }
        )

    dataframe.sort_values(
        by=["statistics"],
        ascending=False,
        inplace=True,
        ignore_index=True,
    )

    if corr_method == "benjamini-hochberg":
        pvals = dataframe["pvals"].to_numpy(dtype=float)
        if pvals.size:
            order = np.argsort(pvals)
            ranks = np.arange(1, pvals.size + 1)
            adjusted = pvals[order] * pvals.size / ranks
            adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
            pvals_adj = np.empty_like(adjusted)
            pvals_adj[order] = np.minimum(adjusted, 1.0)
            dataframe["pvals_adj"] = pvals_adj
        else:
            dataframe["pvals_adj"] = []
    elif corr_method == "bonferroni":
        dataframe["pvals_adj"] = np.minimum(
            dataframe["pvals"] * n_genes,
            1.0,
        )

    if pval_cutoff is not None:
        dataframe = cast(
            pd.DataFrame,
            dataframe[dataframe["pvals_adj"] < pval_cutoff],
        )
    return dataframe


@overload
def smirnov_tests(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    layer: Optional[str] = None,
    alternative: _Alternatives = "two-sided",
    corr_method: _CorrMethod = "benjamini-hochberg",
    pval_cutoff: Optional[float] = None,
    key_added: Optional[str] = None,
    copy: Literal[False] = False,
) -> None: ...


@overload
def smirnov_tests(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    layer: Optional[str] = None,
    alternative: _Alternatives = "two-sided",
    corr_method: _CorrMethod = "benjamini-hochberg",
    pval_cutoff: Optional[float] = None,
    key_added: Optional[str] = None,
    *,
    copy: Literal[True],
) -> pd.DataFrame: ...


@overload
def smirnov_tests(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    layer: Optional[str] = None,
    alternative: _Alternatives = "two-sided",
    corr_method: _CorrMethod = "benjamini-hochberg",
    pval_cutoff: Optional[float] = None,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[pd.DataFrame]: ...


def smirnov_tests(
    adata: AnnData,
    groupby: str,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    layer: Optional[str] = None,
    alternative: _Alternatives = "two-sided",
    corr_method: _CorrMethod = "benjamini-hochberg",
    pval_cutoff: Optional[float] = None,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[pd.DataFrame]:
    """
    Compare whether a subsample and a reference sample have the same distribution
    using a two-sample Kolmogorov-Smirnov test.

    Expects logarithmized data.

    Parameters
    ----------
    adata: AnnData
        Unimodal annotated data matrix.
    groupby: str
        Observation key in `adata.obs` defining groups.
    groups: 'all' or sequence of str (default: "all")
        Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparisons
        shall be restricted, or 'all' for performing comparisons for all groups.
    reference: 'rest' | str
        Name of the group used as reference.
        If 'rest', compare each group to all observations not in the group.
    layer: str, optional
        Layer to use instead of `adata.X`.
    alternative: 'two-sided' | 'less' | 'greater' (default: two-sided)
        Defines the null and alternative hypotheses.
    corr_method: 'benjamini-hochberg' | 'bonferroni'
        Correction method for computing adjusted p-values.
    pval_cutoff: float (optional, default: None)
        Return only adjusted p-values below the cutoff.
    key_added: str, optional
        Key in `adata.uns` where information is saved.
    copy: bool (default: False)
        Return a copy instead of modifying `adata`.

    Returns
    -------
    DataFrame or None
        If `copy=True`, returns the marker ranking results. Otherwise, updates
        `adata` in place and returns None.

        Marker ranking results are stored in:

        - `adata.uns[key_added]`: marker ranking table and metadata.

    Notes
    -----
    Results are stored in `adata.uns[key_added]` with the following values:
    `group`, `names`, `statistics`, `locations`, `signs`, `pvals` and
    `pvals_adj`.

    """

    adata = adata.copy() if copy else adata

    if groups == "all":
        groups = sorted(list(adata.obs[groupby].cat.categories))
    elif isinstance(groups, Seq) and not isinstance(groups, str):
        groups = list(groups)
    else:
        raise TypeError(
            f"unsupported argument type for 'groups': "
            f"expected sequence but received {type(groups)}"
        )

    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        cats = sorted(list(adata.obs[groupby].cat.categories))
        raise ValueError(
            f"invalid argument value for 'reference': "
            f"expected 'rest' or one of {cats} but received {reference!r}"
        )

    if key_added is None:
        key_added = "smirnov_tests"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {
        "groupby": groupby,
        "reference": reference,
        "layer": layer,
        "corr_method": corr_method,
    }

    inputs = _prepare_smirnov_inputs(
        adata,
        groupby,
        groups,
        reference,
        layer,
    )
    if _KS_SUPPORTS_AXIS and _KS_RETURNS_DETAILS and inputs.unique_gene_names:
        results = _vectorized_smirnov_results(inputs, alternative)
    elif _KS_RETURNS_DETAILS:
        results = _scalar_smirnov_results(inputs, alternative)
    else:
        results = _legacy_smirnov_results(inputs, alternative)
    dataframe = _format_smirnov_results(
        inputs,
        results,
        corr_method,
        pval_cutoff,
    )

    if copy:
        return dataframe
    else:
        adata.uns[key_added]["results"] = dataframe
