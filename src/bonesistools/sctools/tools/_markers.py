#!/usr/bin/env python

from __future__ import annotations

from collections.abc import Collection as CollectionInstance
from collections.abc import Mapping as MappingInstance
from collections.abc import Sequence as Seq
from typing import (
    Any,
    Callable,
    Collection,
    Iterable,
    List,
    Mapping,
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
from scipy.stats import hypergeom

from ..._compat import Literal
from ..._validation import _as_boolean, _as_callable, _as_literal, _as_probability
from ..._warnings import _warn_deprecated
from .._typing import VarSubset, anndata_checker
from ._conversion import anndata_to_dataframe
from ._stats import (
    CORRECTION_METHODS,
    CorrectionMethod,
    _adjust_pvalues,
    _resolve_background_groups,
    _resolve_groups,
    welch_tests,
    wilcoxon_tests,
)
from ._utils import _get_expression_with_gene_names

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
    counts_df = anndata_to_dataframe(adata, obs=groupby, layer=layer, is_log=is_log)

    if cluster_rebalancing:
        mean_counts_df = counts_df.groupby(
            by=groupby,
            sort=True,
            observed=False,
        ).mean()

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
        "`bt.sct.tl.calculate_logfoldchanges`",
        replacement="`bt.sct.tl.logfoldchanges`",
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
        Expression source. If None or `"X"`, use `adata.X`. If `"raw.X"`, use
        `adata.raw.X`. Otherwise, interpret as a layer key in `adata.layers`.
    var_subset: str or collection of str, optional
        Variables to test. If a string is provided, it is interpreted as a
        boolean column in `adata.var`. If a collection is provided, it is
        interpreted as variable names. If None, test all variables in the
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
        Approximate maximum memory allocated to dense working arrays. If None,
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
        var_subset,
    ).rename(
        columns={"names": "feature"}
    )
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
                filter_logfoldchanges(
                    cast(np.ndarray, dea_df["logfoldchanges"].values)
                )
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
        )
        group_logfoldchanges = group_logfoldchanges.loc[
            cast(pd.Series, group_logfoldchanges["group"]).astype(bool),
            ["names", "logfoldchanges"],
        ].copy()
        group_logfoldchanges.insert(0, "group", group)
        logfoldchange_tables.append(group_logfoldchanges)

    if not logfoldchange_tables:
        return pd.DataFrame(
            columns=cast(Any, ["group", "names", "logfoldchanges"])
        )

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

    from scipy.sparse import issparse
    from scipy.stats import kstest

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

    if layer is not None:
        expression_mtx = adata.layers[layer]
    else:
        expression_mtx = adata.X

    if issparse(expression_mtx):
        cast(Any, expression_mtx).eliminate_zeros()

    def get_gene_values(obs_mask, gene_name):
        obs_indices = np.where(np.asarray(obs_mask))[0]
        var_indices = np.where(np.asarray(adata.var.index == gene_name))[0]
        if issparse(expression_mtx):
            values = cast(Any, expression_mtx)[obs_indices][:, var_indices].toarray()
        else:
            values = cast(np.ndarray, expression_mtx)[np.ix_(obs_indices, var_indices)]
        return np.asarray(values).reshape(-1)

    df = pd.DataFrame(
        columns=cast(
            Any,
            ["group", "names", "statistics", "locations", "signs", "pvals"],
        )
    )
    index = 0

    for name in adata.var_names:
        reference_sample = (
            get_gene_values(adata.obs[groupby] == reference, name)
            if reference != "rest"
            else None
        )

        for group in groups:
            if reference == "rest":
                ref_sample = get_gene_values(adata.obs[groupby] != group, name)
            else:
                assert reference_sample is not None
                ref_sample = reference_sample
            sample = get_gene_values(adata.obs[groupby] == group, name)
            ks = kstest(sample, ref_sample, alternative=alternative)
            df.loc[index] = [
                group,
                name,
                ks.statistic,
                ks.statistic_location,
                ks.statistic_sign,
                ks.pvalue,
            ]
            index += 1

    df.sort_values(by=["statistics"], ascending=False, inplace=True, ignore_index=True)

    if corr_method == "benjamini-hochberg":
        pvals = df["pvals"].to_numpy(dtype=float)
        if pvals.size:
            order = np.argsort(pvals)
            ranks = np.arange(1, pvals.size + 1)
            adjusted = pvals[order] * pvals.size / ranks
            adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
            pvals_adj = np.empty_like(adjusted)
            pvals_adj[order] = np.minimum(adjusted, 1.0)
            df["pvals_adj"] = pvals_adj
        else:
            df["pvals_adj"] = []
    elif corr_method == "bonferroni":
        n_genes = adata.var.shape[0]
        df["pvals_adj"] = np.minimum(df["pvals"] * n_genes, 1.0)

    if pval_cutoff is not None:
        df = cast(pd.DataFrame, df[df["pvals_adj"] < pval_cutoff])

    if copy:
        return df
    else:
        adata.uns[key_added]["results"] = df
